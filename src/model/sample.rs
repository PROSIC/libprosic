use std::str;
use std::collections::{HashMap, VecDeque, vec_deque};
use std::cmp;
use std::fmt;
use std::error::Error;
use std::f64::consts;
use std::ascii::AsciiExt;
use std::cell::RefCell;

use rgsl::randist::gaussian::{gaussian_pdf, ugaussian_P};
use rgsl::error::erfc;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::{Cigar, CigarStringView};
use bio::stats::{LogProb, PHREDProb, Prob};

use model;
use model::Variant;
use pairhmm::PairHMM;

lazy_static! {
    static ref PROB_CONFUSION: LogProb = LogProb::from(Prob(0.3333));
}


/// Calculate probability of read_base given ref_base.
pub fn prob_read_base(read_base: u8, ref_base: u8, base_qual: u8) -> LogProb {
    let prob_miscall = prob_read_base_miscall(base_qual);

    if read_base.to_ascii_uppercase() == ref_base.to_ascii_uppercase() {
        prob_miscall.ln_one_minus_exp()
    } else {
        // TODO replace the second term with technology specific confusion matrix
        prob_miscall + *PROB_CONFUSION
    }
}


/// Calculate probability of read_base given ref_base.
pub fn prob_read_base_miscall(base_qual: u8) -> LogProb {
    LogProb::from(PHREDProb::from((base_qual) as f64))
}


/// For an SNV, calculate likelihood of ref (first) or alt (second) for a given read, i.e.
/// Pr(Z_i | ref) and Pr(Z_i | alt). This follows Samtools and GATK.
pub fn prob_read_snv(
    record: &bam::Record,
    cigar: &CigarStringView,
    start: u32,
    variant: &Variant,
    ref_seq: &[u8]
) -> Result<(LogProb, LogProb), Box<Error>> {
    if let &Variant::SNV(base) = variant {
        if let Some(qpos) = cigar.read_pos(start, false, false)? {
            let read_base = record.seq()[qpos as usize];
            let base_qual = record.qual()[qpos as usize];
            let prob_alt = prob_read_base(read_base, base, base_qual);
            let prob_ref = prob_read_base(read_base, ref_seq[start as usize], base_qual);
            Ok((prob_ref, prob_alt))
        } else {
            // such a read should be filtered out before
            panic!("bug: read does not enclose SNV");
        }
    } else {
        panic!("bug: unsupported variant");
    }
}


/// Convert MAPQ (from read mapper) to LogProb for the event that the read maps correctly.
pub fn prob_mapping(mapq: u8) -> LogProb {
    LogProb::from(PHREDProb(mapq as f64)).ln_one_minus_exp()
}


/// We assume that the mapping quality of split alignments provided by the aligner is conditional on being no artifact.
/// Artifact alignments can be caused by aligning short ends as splits due to e.g. repeats.
/// We can calculate the probability of having no artifact by investigating if there is at least
/// one full fragment supporting the alternative allele (via enlarged or reduced insert size).
/// Then, the final mapping quality can be obtained by multiplying this probability.
pub fn adjust_mapq(observations: &mut [Observation]) {
    // calculate probability of at least one alt fragment observation
    let prob_no_alt_fragment: LogProb = observations.iter().filter_map(|obs| {
        if !obs.is_alignment_evidence() {
            // Calculate posterior probability of having the alternative allele in that read.
            if obs.prob_alt == LogProb::ln_zero() && obs.prob_ref == LogProb::ln_zero() {
                // Both events are zero. Hence, the read is should be ignored, it is likely another artifact.
                let prob_not_alt = LogProb::ln_one();
                Some(prob_not_alt)
            } else {
                // This assumes that alt and ref are the only possible events, which is reasonable in case of indels.
                let prob_alt = obs.prob_alt - (obs.prob_alt.ln_add_exp(obs.prob_ref));
                let prob_not_alt = (obs.prob_mapping + prob_alt).ln_one_minus_exp();
                Some(prob_not_alt)
            }
        } else {
            None
        }
    }).fold(LogProb::ln_one(), |s, e| s + e);

    let prob_no_artifact = prob_no_alt_fragment.ln_one_minus_exp();
    for obs in observations.iter_mut() {
        if obs.is_alignment_evidence() && obs.prob_alt > obs.prob_ref {
            // adjust as Pr(mapping) = Pr(no artifact) * Pr(mapping|no artifact) + Pr(artifact) * Pr(mapping|artifact)
            // with Pr(mapping|artifact) = 0
            obs.prob_mapping = prob_no_artifact + obs.prob_mapping;
        }
    }
}


quick_error! {
    #[derive(Debug)]
    pub enum RecordBufferError {
        UnknownSequence(chrom: String) {
            description("unknown sequence")
            display("sequence {} cannot be found in BAM", chrom)
        }
    }
}


/// Ringbuffer of BAM records. This data structure ensures that no bam record is read twice while
/// extracting observations for given variants.
pub struct RecordBuffer {
    reader: bam::IndexedReader,
    inner: VecDeque<bam::Record>,
    window: u32,
    use_secondary: bool
}


unsafe impl Sync for RecordBuffer {}
unsafe impl Send for RecordBuffer {}


impl RecordBuffer {
    /// Create a new `RecordBuffer`.
    pub fn new(bam: bam::IndexedReader, window: u32, use_secondary: bool) -> Self {
        RecordBuffer {
            reader: bam,
            inner: VecDeque::with_capacity(window as usize * 2),
            window: window as u32,
            use_secondary: use_secondary
        }
    }

    /// Return end position of buffer.
    fn end(&self) -> Option<u32> {
        self.inner.back().map(|rec| rec.pos() as u32)
    }

    fn tid(&self) -> Option<i32> {
        self.inner.back().map(|rec| rec.tid())
    }

    /// Fill buffer around the given interval.
    pub fn fill(&mut self, chrom: &[u8], start: u32, end: u32) -> Result<(), Box<Error>> {
        if let Some(tid) = self.reader.header.tid(chrom) {
            let window_start = cmp::max(start as i32 - self.window as i32 - 1, 0) as u32;
            if self.inner.is_empty() || self.end().unwrap() < window_start || self.tid().unwrap() != tid as i32 {
                let end = self.reader.header.target_len(tid).unwrap();
                try!(self.reader.seek(tid, window_start, end));
                debug!("Clearing ringbuffer");
                self.inner.clear();
            } else {
                // remove records too far left
                let to_remove = self.inner.iter().take_while(|rec| rec.pos() < window_start as i32).count();
                debug!("Removing {} records", to_remove);
                for _ in 0..to_remove {
                    self.inner.pop_front();
                }
            }

            // extend to the right
            loop {
                let mut record = bam::Record::new();
                if let Err(e) = self.reader.read(&mut record) {
                    if e.is_eof() {
                        break;
                    }
                    return Err(Box::new(e));
                }

                let pos = record.pos();
                if record.is_duplicate() || record.is_unmapped() {
                    continue;
                }
                if !self.use_secondary && record.is_secondary() {
                    continue;
                }
                self.inner.push_back(record);
                if pos > end as i32 + self.window as i32 {
                    break;
                }
            }

            debug!("New buffer length: {}", self.inner.len());

            Ok(())
        } else {
            Err(Box::new(RecordBufferError::UnknownSequence(str::from_utf8(chrom).unwrap().to_owned())))
        }
    }

    pub fn iter(&self) -> vec_deque::Iter<bam::Record> {
        self.inner.iter()
    }
}


/// An observation for or against a variant.
#[derive(Clone, Debug, RustcDecodable, RustcEncodable)]
pub struct Observation {
    /// Posterior probability that the read/read-pair has been mapped correctly (1 - MAPQ).
    pub prob_mapping: LogProb,
    /// Probability that the read/read-pair comes from the alternative allele.
    pub prob_alt: LogProb,
    /// Probability that the read/read-pair comes from the reference allele.
    pub prob_ref: LogProb,
    /// Probability of the read/read-pair given that it has been mismapped.
    pub prob_mismapped: LogProb,
    /// Type of evidence.
    pub evidence: Evidence
}


impl Observation {
    pub fn is_alignment_evidence(&self) -> bool {
        if let Evidence::Alignment(_) = self.evidence {
            true
        } else {
            false
        }
    }
}


#[derive(Clone, Debug, RustcDecodable, RustcEncodable)]
pub enum Evidence {
    /// Insert size of fragment
    InsertSize(u32),
    /// Alignment of a single read
    Alignment(String)
}


impl Evidence {
    pub fn dummy_alignment() -> Self {
        Evidence::Alignment("Dummy-Alignment".to_owned())
    }
}


impl<'a, CigarString> From<&'a CigarString> for Evidence
where CigarString: fmt::Display {
    fn from(cigar: &CigarString) -> Self {
        Evidence::Alignment(format!("{}", cigar))
    }
}


/// Expected insert size in terms of mean and standard deviation.
/// This should be estimated from unsorted(!) bam files to avoid positional biases.
#[derive(Copy, Clone, Debug)]
pub struct InsertSize {
    pub mean: f64,
    pub sd: f64
}


/// A sequenced sample, e.g., a tumor or a normal sample.
pub struct Sample {
    record_buffer: RecordBuffer,
    use_fragment_evidence: bool,
    use_mapq: bool,
    adjust_mapq: bool,
    insert_size: InsertSize,
    likelihood_model: model::likelihood::LatentVariableModel,
    prob_spurious_isize: LogProb,
    prob_insertion_artifact: LogProb,
    prob_deletion_artifact: LogProb,
    prob_insertion_extend_artifact: LogProb,
    prob_deletion_extend_artifact: LogProb,
    max_indel_overlap: u32,
    indel_haplotype_window: u32,
    pair_hmm: RefCell<PairHMM>
}


impl Sample {
    /// Create a new `Sample`.
    ///
    /// # Arguments
    ///
    /// * `bam` - BAM file with the aligned and deduplicated sequence reads.
    /// * `pileup_window` - Window around the variant that shall be search for evidence (e.g. 5000).
    /// * `use_fragment_evidence` - Whether to use read pairs that are left and right of variant.
    /// * `use_secondary` - Whether to use secondary alignments.
    /// * `insert_size` - estimated insert size
    /// * `prior_model` - Prior assumptions about allele frequency spectrum of this sample.
    /// * `likelihood_model` - Latent variable model to calculate likelihoods of given observations.
    /// * `prob_spurious_isize` - rate of wrongly reported insert size abberations (mapper dependent, BWA: 0.01332338, LASER: 0.05922201)
    /// * `max_indel_overlap` - maximum number of bases a read may be aligned beyond the start or end of an indel in order to be considered as an observation
    /// * `indel_haplotype_window` - maximum number of considered bases around an indel breakpoint
    pub fn new(
        bam: bam::IndexedReader,
        pileup_window: u32,
        use_fragment_evidence: bool,
        use_secondary: bool,
        use_mapq: bool,
        adjust_mapq: bool,
        insert_size: InsertSize,
        likelihood_model: model::likelihood::LatentVariableModel,
        prob_spurious_isize: Prob,
        prob_insertion_artifact: Prob,
        prob_deletion_artifact: Prob,
        prob_insertion_extend_artifact: Prob,
        prob_deletion_extend_artifact: Prob,
        max_indel_overlap: u32,
        indel_haplotype_window: u32
    ) -> Self {
        Sample {
            record_buffer: RecordBuffer::new(bam, pileup_window, use_secondary),
            use_fragment_evidence: use_fragment_evidence,
            use_mapq: use_mapq,
            adjust_mapq: adjust_mapq,
            insert_size: insert_size,
            likelihood_model: likelihood_model,
            prob_spurious_isize: LogProb::from(prob_spurious_isize),
            prob_insertion_artifact: LogProb::from(prob_insertion_artifact),
            prob_deletion_artifact: LogProb::from(prob_deletion_artifact),
            prob_insertion_extend_artifact: LogProb::from(prob_insertion_extend_artifact),
            prob_deletion_extend_artifact: LogProb::from(prob_deletion_extend_artifact),
            max_indel_overlap: max_indel_overlap,
            indel_haplotype_window: indel_haplotype_window,
            pair_hmm: RefCell::new(PairHMM::new())
        }
    }

    /// Return likelihood model.
    pub fn likelihood_model(&self) -> model::likelihood::LatentVariableModel {
        self.likelihood_model
    }

    /// Extract observations for the given variant.
    pub fn extract_observations(
        &mut self,
        chrom: &[u8],
        start: u32,
        variant: &Variant,
        chrom_seq: &[u8]
    ) -> Result<Vec<Observation>, Box<Error>> {
        let mut observations = Vec::new();
        let (end, centerpoint) = match variant {
            &Variant::Deletion(length)  => (start + length, start + length / 2),
            &Variant::Insertion(_) => (start, start),
            &Variant::SNV(_) => (start, start)
        };
        let mut pairs = HashMap::new();
        let mut n_overlap = 0;

        debug!("variant: {}:{} {:?}", str::from_utf8(chrom).unwrap(), start, variant);

        // move window to the current variant
        debug!("Filling buffer...");
        try!(self.record_buffer.fill(chrom, start, end));
        debug!("Done.");


        match variant {
            &Variant::SNV(_) => {
                // iterate over records
                for record in self.record_buffer.iter() {
                    let pos = record.pos() as u32;
                    let cigar = record.cigar();
                    let end_pos = cigar.end_pos() as u32;

                    if pos <= start && end_pos >= start {
                        let cigar = record.cigar();
                        observations.push(
                            self.read_observation(&record, &cigar, start, variant, chrom_seq)?
                        );
                        n_overlap += 1;
                    }
                }
            },
            &Variant::Insertion(_) | &Variant::Deletion(_) => {

                // iterate over records
                for record in self.record_buffer.iter() {
                    let cigar = record.cigar();
                    let pos = record.pos() as u32;
                    let end_pos = cigar.end_pos() as u32;

                    let overlap = {
                        // consider soft clips for overlap detection
                        let pos = if let Cigar::SoftClip(l) = cigar[0] {
                            pos.saturating_sub(l)
                        } else {
                            pos
                        };
                        let end_pos = if let Cigar::SoftClip(l) = cigar[cigar.len() - 1] {
                            end_pos + l
                        } else {
                            end_pos
                        };

                        if end_pos <= end {
                            cmp::min(end_pos.saturating_sub(start), variant.len())
                        } else {
                            cmp::min(end.saturating_sub(pos), variant.len())
                        }
                    };

                    // read evidence
                    if overlap > 0 {
                        if overlap <= self.max_indel_overlap {
                            observations.push(
                                self.read_observation(&record, &cigar, start, variant, chrom_seq)?
                            );
                            n_overlap += 1;
                        }
                    }

                    // fragment evidence
                    // We have to consider the fragment even if the read has been used for the
                    // overlap evidence above. Otherwise, results would be biased towards
                    // alternative alleles (since reference reads tend to overlap the centerpoint).
                    if self.use_fragment_evidence &&
                       (record.is_first_in_template() || record.is_last_in_template()) {
                        // We ensure fair sampling by checking if the whole fragment overlaps the
                        // centerpoint. Only taking the internal segment would not be fair,
                        // because then the second read of reference fragments tends to cross
                        // the centerpoint and the fragment would be discarded.
                        // The latter would not happen for alt fragments, because the second read
                        // would map right of the variant in that case.
                        if pos <= centerpoint {
                            // need to check mate
                            // since the bam file is sorted by position, we can't see the mate first
                            let insert_size = record.insert_size().abs() as u32;
                            if pos + insert_size >= centerpoint {
                                pairs.insert(record.qname().to_owned(), record.mapq());
                            }
                        } else if let Some(mate_mapq) = pairs.get(record.qname()) {
                            // mate already visited, and this fragment overlaps centerpoint
                            observations.push(self.fragment_observation(&record, *mate_mapq, variant));
                        }
                    }
                }
            }
        }

        debug!("Extracted observations ({} fragments, {} overlapping reads).", pairs.len(), n_overlap);

        if self.adjust_mapq {
            match variant {
                // only adjust for deletion and insertion
                &Variant::Deletion(_) | &Variant::Insertion(_) => adjust_mapq(&mut observations),
                _ => ()
            }

        }
        Ok(observations)
    }

    fn prob_mapping(&self, mapq: u8) -> LogProb {
        if self.use_mapq {
            prob_mapping(mapq)
        } else {
            LogProb::ln_one()
        }
    }

    fn read_observation(
        &self,
        record: &bam::Record,
        cigar: &CigarStringView,
        start: u32,
        variant: &Variant,
        chrom_seq: &[u8]
    ) -> Result<Observation, Box<Error>> {
        let prob_mapping = self.prob_mapping(record.mapq());
        debug!("prob_mapping={}", *prob_mapping);

        let (prob_ref, prob_alt) = match variant {
            &Variant::Deletion(_) | &Variant::Insertion(_) => {
                self.prob_read_indel(record, cigar, start, variant, chrom_seq)?
            },
            &Variant::SNV(_) => {
                prob_read_snv(record, cigar, start, variant, chrom_seq)?
            }
        };

        Ok(Observation {
            prob_mapping: prob_mapping,
            prob_alt: prob_alt,
            prob_ref: prob_ref,
            prob_mismapped: LogProb::ln_one(), // if the read is mismapped, we assume sampling probability 1.0
            evidence: Evidence::from(cigar)
        })
    }

    fn fragment_observation(
        &self,
        record: &bam::Record,
        mate_mapq: u8,
        variant: &Variant
    ) -> Observation {
        let insert_size = record.insert_size().abs();
        let shift = match variant {
            &Variant::Deletion(_)  => variant.len() as f64,
            &Variant::Insertion(_) => -(variant.len() as f64),
            &Variant::SNV(_) => panic!("no fragment observations for SNV")
        };
        let p_alt = (
            // case: correctly called indel
            self.prob_spurious_isize.ln_one_minus_exp() + isize_pmf(
                insert_size as f64,
                self.insert_size.mean + shift,
                self.insert_size.sd
            )
        ).ln_add_exp(
            // case: no indel, false positive call
            self.prob_spurious_isize +
            isize_pmf(
                insert_size as f64,
                self.insert_size.mean,
                self.insert_size.sd
            )
        );

        let obs = Observation {
            prob_mapping: self.prob_mapping(record.mapq()) + self.prob_mapping(mate_mapq),
            prob_alt: p_alt,
            prob_ref: isize_pmf(insert_size as f64, self.insert_size.mean, self.insert_size.sd),
            prob_mismapped: LogProb::ln_one(), // if the fragment is mismapped, we assume sampling probability 1.0
            evidence: Evidence::InsertSize(insert_size as u32)
        };

        obs
    }

    /// For an indel, calculate likelihood of ref (first) or alt (second) for a given read,
    /// i.e. Pr(Z_i | ref) and Pr(Z_i | alt).
    pub fn prob_read_indel(
        &self,
        record: &bam::Record,
        cigar: &CigarStringView,
        start: u32,
        variant: &Variant,
        chrom_seq: &[u8]
    ) -> Result<(LogProb, LogProb), Box<Error>> {
        // TODO we really need the insertion sequence!!
        // Otherwise we loose a lot of discriminative power.

        let read_seq = record.seq();
        let read_qual = record.qual();

        let (read_offset, read_end, breakpoint) = {
            let (varstart, varend) = match variant {
                &Variant::Deletion(_) => (start, start + variant.len()),
                &Variant::Insertion(_) => (start, start + 1),
                &Variant::SNV(_) => panic!("bug: unsupported variant")
            };

            match (
                cigar.read_pos(varstart, true, true)?,
                cigar.read_pos(varend, true, true)?
            ) {
                // read encloses variant
                (Some(qstart), Some(qend)) => {
                    let qstart = qstart as usize;
                    let qend = qend as usize;
                    let read_offset = qstart.saturating_sub(self.indel_haplotype_window as usize);
                    let read_end = cmp::min(
                        qend + self.indel_haplotype_window as usize,
                        read_seq.len()
                    );
                    //println!("window: {} - {} (variant: {} - {})", read_offset, read_end, qstart, qend);
                    (read_offset, read_end, varstart as usize)
                },
                (Some(qstart), None) => {
                    let qstart = qstart as usize;
                    let read_offset = qstart.saturating_sub(self.indel_haplotype_window as usize);
                    let read_end = cmp::min(
                        qstart + self.indel_haplotype_window as usize,
                        read_seq.len()
                    );
                    (read_offset, read_end, varstart as usize)
                },
                (None, Some(qend)) => {
                    let qend = qend as usize;
                    let read_offset = qend.saturating_sub(self.indel_haplotype_window as usize);
                    let read_end = cmp::min(
                        qend + self.indel_haplotype_window as usize,
                        read_seq.len()
                    );
                    (read_offset, read_end, varend as usize)
                },
                (None, None) => {
                    panic!(
                        "bug: read does not overlap breakpoint: pos={}, cigar={}, start={}, len={}",
                        record.pos(),
                        cigar,
                        start,
                        variant.len()
                    );
                }
            }
        };

        let start = start as usize;
        // the window on the reference could be a bit larger to allow some flexibility with close
        // indels. But it should not be so large that the read can align outside of the breakpoint.
        let ref_window = (self.indel_haplotype_window as f64 * 1.5) as usize;


        let prob_emit_xy = |ref_base, j| {
            let j_ = j + read_offset;
            let read_base = read_seq[j_];
            prob_read_base(read_base, ref_base, read_qual[j_])
        };
        let prob_emit_x = |_| LogProb::ln_one();
        let prob_emit_y = |j| prob_read_base_miscall(read_qual[j + read_offset]);

        // ref allele
        let ref_offset = breakpoint.saturating_sub(ref_window);
        let ref_end = cmp::min(breakpoint + ref_window, chrom_seq.len());

        let prob_ref = self.pair_hmm.borrow_mut().prob_related(
            self.prob_insertion_artifact,
            self.prob_deletion_artifact,
            self.prob_insertion_extend_artifact,
            self.prob_deletion_extend_artifact,
            &|i, j| prob_emit_xy(chrom_seq[i + ref_offset], j),
            &prob_emit_x,
            &prob_emit_y,
            ref_end - ref_offset,
            read_end - read_offset,
            true,
            true
        );

        // alt allele
        let prob_alt = match variant {
            &Variant::Deletion(_) => {
                let l = variant.len() as usize;

                let ref_offset = start.saturating_sub(ref_window);
                let ref_end = cmp::min(start + ref_window, chrom_seq.len());
                let project_i = |i| {
                    let i_ = i + ref_offset;
                    if i_ <= start {
                        i_
                    } else {
                        (i_ + l)
                    }
                };

                self.pair_hmm.borrow_mut().prob_related(
                    self.prob_insertion_artifact,
                    self.prob_deletion_artifact,
                    self.prob_insertion_extend_artifact,
                    self.prob_deletion_extend_artifact,
                    &|i, j| {
                        prob_emit_xy(chrom_seq[project_i(i)], j)
                    },
                    &prob_emit_x,
                    &prob_emit_y,
                    ref_end - ref_offset,
                    read_end - read_offset,
                    true,
                    true
                )
            },
            &Variant::Insertion(ref ins_seq) => {
                let l = ins_seq.len() as usize;

                let ref_offset = start.saturating_sub(ref_window) as usize;
                let ref_end = cmp::min(start + l + ref_window, chrom_seq.len());

                self.pair_hmm.borrow_mut().prob_related(
                    self.prob_insertion_artifact,
                    self.prob_deletion_artifact,
                    self.prob_insertion_extend_artifact,
                    self.prob_deletion_extend_artifact,
                    &|i, j| {
                        let i_ = i + ref_offset;
                        let ref_base = if i_ <= start {
                            // before variant
                            chrom_seq[i_]
                        } else if i_ > start + l {
                            // after variant
                            chrom_seq[i_ - l]
                        } else {
                            // within variant
                            ins_seq[i_ - (start + 1)]
                        };
                        prob_emit_xy(ref_base, j)
                    },
                    &prob_emit_x,
                    &prob_emit_y,
                    ref_end - ref_offset,
                    read_end - read_offset,
                    true,
                    true
                )
            },
            _ => {
                panic!("bug: unsupported variant");
            }
        };

        Ok((prob_ref, prob_alt))
    }
}


/// as shown in http://www.milefoot.com/math/stat/pdfc-normaldisc.htm
pub fn isize_pmf(value: f64, mean: f64, sd: f64) -> LogProb {
    // TODO fix density in paper
    LogProb(
        (
            ugaussian_P((value + 0.5 - mean) / sd) -
            ugaussian_P((value - 0.5 - mean) / sd)
        ).ln()// - ugaussian_P(-mean / sd).ln()
    )
}


/// Continuous normal density (obsolete).
pub fn isize_density(value: f64, mean: f64, sd: f64) -> LogProb {
    LogProb(gaussian_pdf(value - mean, sd).ln())
}


/// Manual normal density (obsolete, we can use GSL (see above)).
pub fn isize_density_louis(value: f64, mean: f64, sd: f64) -> LogProb {
    let mut p = 0.5 / (1.0 - 0.5 * erfc((mean + 0.5) / sd * consts::FRAC_1_SQRT_2));
    p *= erfc((-value - 0.5 + mean)/sd * consts::FRAC_1_SQRT_2) - erfc((-value + 0.5 + mean)/sd * consts::FRAC_1_SQRT_2);

    LogProb(p.ln())
}


pub fn isize_mixture_density_louis(value: f64, d: f64, mean: f64, sd: f64, rate: f64) -> LogProb {
    let p = 0.5 / ( rate*(1.0 - 0.5*erfc((mean + 0.5)/sd*consts::FRAC_1_SQRT_2)) + (1.0 - rate)*(1.0 - 0.5*erfc((mean + d + 0.5)/sd* consts::FRAC_1_SQRT_2)) );
    LogProb((p * (
        rate*( erfc((-value - 0.5 + mean)/sd*consts::FRAC_1_SQRT_2) - erfc((-value + 0.5 + mean)/sd*consts::FRAC_1_SQRT_2) ) + (1.0 - rate)*( erfc((-value - 0.5 + mean + d)/sd*consts::FRAC_1_SQRT_2) - erfc((-value + 0.5 + mean + d)/sd*consts::FRAC_1_SQRT_2) )
    )).ln())
}

#[cfg(test)]
mod tests {
    extern crate env_logger;

    use super::*;
    use model;
    use likelihood;
    use constants;

    use csv;
    use std::str;
    use itertools::Itertools;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use bio::stats::{LogProb, PHREDProb, Prob};
    use bio::io::fasta;


    fn read_observations(path: &str) -> Vec<Observation> {
        let mut reader = csv::Reader::from_file(path).expect("error reading example").delimiter(b'\t');
        let obs = reader.decode().collect::<Result<Vec<(String, u32, u32, String, Observation)>, _>>().unwrap();
        let groups = obs.into_iter().group_by(|&(_, _, _, ref sample, _)| {
            sample == "case"
        });
        let case_obs = groups.into_iter().next().unwrap().1.into_iter().map(|(_, _, _, _, obs)| obs).collect_vec();
        case_obs
    }


    #[test]
    fn test_adjust_mapq_with_fragment_evidence() {
        let mut observations = vec![
            Observation {
                prob_mapping: LogProb(0.5f64.ln()),
                prob_alt: LogProb::ln_one(),
                prob_ref: LogProb::ln_zero(),
                prob_mismapped: LogProb::ln_one(),
                evidence: Evidence::dummy_alignment()
            },
            Observation {
                prob_mapping: LogProb::ln_one(),
                prob_alt: LogProb::ln_one(),
                prob_ref: LogProb::ln_zero(),
                prob_mismapped: LogProb::ln_one(),
                evidence: Evidence::InsertSize(300)
            },
            Observation {
                prob_mapping: LogProb::ln_one(),
                prob_alt: LogProb::ln_zero(),
                prob_ref: LogProb::ln_one(),
                prob_mismapped: LogProb::ln_one(),
                evidence: Evidence::InsertSize(300)
            }
        ];

        adjust_mapq(&mut observations);
        println!("{:?}", observations);
        assert_relative_eq!(*observations[0].prob_mapping, *LogProb(0.5f64.ln()));
    }

    #[test]
    fn test_adjust_mapq_without_fragment_evidence() {
        let mut observations = vec![
            Observation {
                prob_mapping: LogProb(0.5f64.ln()),
                prob_alt: LogProb::ln_one(),
                prob_ref: LogProb::ln_zero(),
                prob_mismapped: LogProb::ln_one(),
                evidence: Evidence::dummy_alignment()
            },
            Observation {
                prob_mapping: LogProb::ln_one(),
                prob_alt: LogProb::ln_zero(),
                prob_ref: LogProb::ln_one(),
                prob_mismapped: LogProb::ln_one(),
                evidence: Evidence::InsertSize(300)
            },
            Observation {
                prob_mapping: LogProb::ln_one(),
                prob_alt: LogProb::ln_zero(),
                prob_ref: LogProb::ln_one(),
                prob_mismapped: LogProb::ln_one(),
                evidence: Evidence::InsertSize(300)
            }
        ];

        adjust_mapq(&mut observations);
        println!("{:?}", observations);
        assert_relative_eq!(*observations[0].prob_mapping, *LogProb(0.0f64.ln()));
    }

    #[test]
    fn test_adjust_mapq_weak_fragment_evidence() {
        let mut observations = vec![
            Observation {
                prob_mapping: LogProb(0.5f64.ln()),
                prob_alt: LogProb::ln_one(),
                prob_ref: LogProb::ln_zero(),
                prob_mismapped: LogProb::ln_one(),
                evidence: Evidence::dummy_alignment()
            },
            Observation {
                prob_mapping: LogProb::ln_one(),
                prob_alt: LogProb(0.5f64.ln()),
                prob_ref: LogProb(0.5f64.ln()),
                prob_mismapped: LogProb::ln_one(),
                evidence: Evidence::InsertSize(300)
            },
            Observation {
                prob_mapping: LogProb::ln_one(),
                prob_alt: LogProb::ln_zero(),
                prob_ref: LogProb::ln_one(),
                prob_mismapped: LogProb::ln_one(),
                evidence: Evidence::InsertSize(300)
            }
        ];

        adjust_mapq(&mut observations);
        println!("{:?}", observations);
        assert_eq!(*observations[0].prob_mapping, *LogProb(0.25f64.ln()));
    }

    #[test]
    fn test_adjust_mapq_real() {
        let mut observations = read_observations("tests/example7.obs.txt");

        adjust_mapq(&mut observations);
        for obs in observations {
            println!("{:?}", obs);
        }
    }

    #[test]
    fn test_isize_density() {
        let d1 = isize_density(300.0, 312.0, 15.0);
        let d2 = isize_pmf(300.0, 312.0, 15.0);
        let d3 = isize_density_louis(300.0, 312.0, 15.0);
        println!("{} {} {}", *d1, *d2, *d3);

        let d_mix = isize_mixture_density_louis(212.0, -100.0, 312.0, 15.0, 0.05);
        let rate = LogProb(0.05f64.ln());
        let p_alt = (
            // case: correctly called indel
            rate.ln_one_minus_exp() + isize_pmf(
                212.0,
                312.0 - 100.0,
                15.0
            )
        ).ln_add_exp(
            // case: no indel, false positive call
            rate +
            isize_pmf(
                212.0,
                312.0,
                15.0
            )
        );

        println!("{} {}", d_mix.exp(), p_alt.exp());
    }

    #[test]
    fn test_prob_mapping() {
        assert_relative_eq!(prob_mapping(0).exp(), 0.0);
        assert_relative_eq!(prob_mapping(10).exp(), 0.9);
        assert_relative_eq!(prob_mapping(20).exp(), 0.99);
        assert_relative_eq!(prob_mapping(30).exp(), 0.999);
    }

    fn setup_sample(isize_mean: f64) -> Sample {
        Sample::new(
            bam::IndexedReader::from_path(&"tests/indels.bam").unwrap(),
            2500,
            true,
            true,
            true,
            false,
            InsertSize { mean: isize_mean, sd: 20.0 },
            likelihood::LatentVariableModel::new(1.0),
            Prob(0.0),
            constants::PROB_ILLUMINA_INS,
            constants::PROB_ILLUMINA_DEL,
            Prob(0.0),
            Prob(0.0),
            20,
            10
        )
    }

    #[test]
    fn test_read_observation_indel() {
        let variant = model::Variant::Insertion(b"GCATCCTGCG".to_vec());
        // insertion starts at 546 and has length 10
        let varpos = 546;

        let sample = setup_sample(150.0);
        let bam = bam::Reader::from_path(&"tests/indels.bam").unwrap();
        let records = bam.records().collect_vec();

        let ref_seq = ref_seq();

        let true_alt_probs = [-0.09, -0.02, -73.09, -16.95, -73.09];

        for (record, true_alt_prob) in records.into_iter().zip(true_alt_probs.into_iter()) {
            let record = record.unwrap();
            let cigar = record.cigar();
            let obs = sample.read_observation(&record, &cigar, varpos, &variant, &ref_seq).unwrap();
            println!("{:?}", obs);
            assert_relative_eq!(*obs.prob_alt, *true_alt_prob, epsilon=0.01);
            assert_relative_eq!(*obs.prob_mapping, *(LogProb::from(PHREDProb(60.0)).ln_one_minus_exp()));
            assert_relative_eq!(*obs.prob_mismapped, *LogProb::ln_one());
        }
    }

    #[test]
    fn test_fragment_observation_no_evidence() {
        let sample = setup_sample(150.0);
        let bam = bam::Reader::from_path(&"tests/indels.bam").unwrap();
        let records = bam.records().map(|rec| rec.unwrap()).collect_vec();

        for varlen in &[0, 5, 10, 100] {
            println!("varlen {}", varlen);
            println!("insertion");
            let variant = model::Variant::Insertion(vec![b'A'; *varlen]);
            for record in &records {
                let obs = sample.fragment_observation(record, 60u8, &variant);
                println!("{:?}", obs);
                if *varlen == 0 {
                    assert_relative_eq!(*obs.prob_ref, *obs.prob_alt);
                } else {
                    assert!(obs.prob_ref > obs.prob_alt);
                }
            }
            println!("deletion");
            let variant = model::Variant::Deletion(*varlen as u32);
            for record in &records {
                let obs = sample.fragment_observation(record, 60u8, &variant);
                println!("{:?}", obs);
                if *varlen == 0 {
                    assert_relative_eq!(*obs.prob_ref, *obs.prob_alt);
                } else {
                    assert!(obs.prob_ref > obs.prob_alt);
                }
            }
        }
    }

    #[test]
    fn test_fragment_observation_evidence() {
        let bam = bam::Reader::from_path(&"tests/indels.bam").unwrap();
        let records = bam.records().map(|rec| rec.unwrap()).collect_vec();

        println!("deletion");
        let sample = setup_sample(100.0);
        let variant = model::Variant::Deletion(50);
        for record in &records {
            let obs = sample.fragment_observation(record, 60u8, &variant);
            println!("{:?}", obs);
            assert_relative_eq!(obs.prob_ref.exp(), 0.0, epsilon=0.001);
            assert!(obs.prob_alt > obs.prob_ref);
        }

        println!("insertion");
        let sample = setup_sample(200.0);
        let variant = model::Variant::Insertion(vec![b'A'; 50]);
        for record in &records {
            let obs = sample.fragment_observation(record, 60u8, &variant);
            println!("{:?}", obs);
            assert_relative_eq!(obs.prob_ref.exp(), 0.0, epsilon=0.001);
            assert!(obs.prob_alt > obs.prob_ref);
        }
    }

    #[test]
    fn test_record_buffer() {
        let bam = bam::IndexedReader::from_path(&"tests/indels.bam").unwrap();
        let mut buffer = RecordBuffer::new(bam, 10, true);

        buffer.fill(b"17", 10, 20).unwrap();
        buffer.fill(b"17", 478, 500).unwrap();
        buffer.fill(b"17", 1000, 1700).unwrap();
        // TODO add assertions
    }

    fn ref_seq() -> Vec<u8> {
        let mut fa = fasta::Reader::from_file(&"tests/chr17.prefix.fa").unwrap();
        let mut chr17 = fasta::Record::new();
        fa.read(&mut chr17).unwrap();

        chr17.seq().to_owned()
    }

    #[test]
    fn test_prob_read_indel() {
        let _ = env_logger::init();

        let bam = bam::Reader::from_path(&"tests/indels+clips.bam").unwrap();
        let records = bam.records().map(|rec| rec.unwrap()).collect_vec();
        let ref_seq = ref_seq();
        let window = 100;
        let sample = setup_sample(312.0);

        // truth
        let probs_alt = [-0.09, -16.95, -73.09, -0.022, -0.011, -0.03];
        let probs_ref = [-151.13, -163.03, -0.01, -67.75, -67.74, -67.76];

        // variant (obtained via bcftools)
        let start = 546;
        let variant = model::Variant::Insertion(b"GCATCCTGCG".to_vec());
        for (i, rec) in records.iter().enumerate() {
            println!("{}", str::from_utf8(rec.qname()).unwrap());
            let (prob_ref, prob_alt) = sample.prob_read_indel(rec, &rec.cigar(), start, &variant, &ref_seq).unwrap();
            println!("Pr(ref)={} Pr(alt)={}", *prob_ref, *prob_alt);
            println!("{:?}", rec.cigar());
            assert_relative_eq!(*prob_ref, probs_ref[i], epsilon=0.1);
            assert_relative_eq!(*prob_alt, probs_alt[i], epsilon=0.1);
        }
    }
}

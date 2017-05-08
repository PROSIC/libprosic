use std::str;
use std::collections::{HashMap, VecDeque, vec_deque};
use std::cmp;
use std::error::Error;
use std::f64::consts;
use log::LogLevel::Debug;
use std::ascii::AsciiExt;

use rgsl::randist::gaussian::{gaussian_pdf, ugaussian_P};
use rgsl::error::erfc;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::Cigar;
use bio::stats::{LogProb, PHREDProb, Prob};
use itertools;

use model;
use model::Variant;


lazy_static! {
    static ref PROB_CONFUSION: LogProb = LogProb::from(Prob(0.3333));
}


/// Calculate probability of read_base given ref_base.
pub fn prob_read_base(read_base: u8, ref_base: u8, base_qual: u8) -> LogProb {
    let prob_miscall = LogProb::from(PHREDProb::from(base_qual as f64));

    if read_base.to_ascii_uppercase() == ref_base.to_ascii_uppercase() {
        prob_miscall.ln_one_minus_exp()
    } else {
        // TODO replace the second term with technology specific confusion matrix
        prob_miscall + *PROB_CONFUSION
    }
}

/// For an SNV, calculate likelihood of ref (first) or alt (second) for a given read, i.e. Pr(Z_i | ref) and Pr(Z_i | alt).
/// This follows Samtools and GATK.
/// TODO: Add cigar string parsing tests for SNVs.
/// TODO: Implement more general cigar parsing function analogous to htslib resolve_cigar2 (sam.c) somewhere more central in libprosic or rust-htslib.
pub fn prob_read_snv(record: &bam::Record, cigar: &[Cigar], vpos: u32, variant: Variant, ref_seq: &[u8]) -> Option<(LogProb, LogProb)> {
    let contains_vpos = |pos: u32, cigar_length: u32|  pos <= vpos && pos + cigar_length > vpos;

    if let Variant::SNV(base) = variant {
        let mut rpos = record.pos() as u32; // reference position
        let mut qpos = 0u32; // position within query / read / observation
        let mut c = 0; // index into cigar operation vector
        // find first cigar operation referring to qpos = 0 (and thus bases in record.seq()), because
        // all augmentations of qpos and rpos before that are invalid
        for i in 0..cigar.len() {
            match &cigar[i] {
                &Cigar::Match(_) |
                &Cigar::Diff(_)  |
                &Cigar::Equal(_) |
                &Cigar::Ins(_)   | // this is unexpected, but bwa + GATK indel realignment can produce insertions before matching positions
                &Cigar::SoftClip(_) => {
                    c = i;
                    break;
                },
                &Cigar::Del(_) => panic!("Cigar String: Unexpected Deletion operation found before any reference match."),
                &Cigar::Back(_) => panic!("Back: Unsupported Cigar operation (not in SAMv1 spec)."),
                &Cigar::Pad(_) => panic!("Pad: Unsupported Cigar operation (unclear SAMv1 spec)."),
                &Cigar::RefSkip(_) => panic!("Cigar String: Unexpected Reference Skip operation found before any reference match."),
                &Cigar::HardClip(_) => ()
            }

        }
        while rpos <= vpos {
            match &cigar[c] {
                // potential SNV evidence
                &Cigar::Match(l) | &Cigar::Diff(l) | &Cigar::Equal(l) if contains_vpos(rpos, l) => {
                    qpos += vpos - rpos; // difference between vpos and first position of current cigar operation
                    let read_base = record.seq()[qpos as usize];
                    let base_qual = record.qual()[qpos as usize];
                    let prob_alt = prob_read_base(read_base, base, base_qual);
                    let prob_ref = prob_read_base(read_base, ref_seq[vpos as usize], base_qual);
                    return Some( (prob_ref, prob_alt) );
                },
                // for others, just increase pos and qpos as needed
                &Cigar::Match(l)   |
                &Cigar::Diff(l)    |
                &Cigar::Equal(l)   => {
                    rpos += l;
                    qpos += l;
                    c += 1;
                },
                &Cigar::SoftClip(l) |
                &Cigar::Ins(l)  => {
                    qpos += l;
                    c += 1;
                },
                &Cigar::RefSkip(l) |
                &Cigar::Del(l)  => {
                    rpos += l;
                    c += 1;
                },
                &Cigar::HardClip(_) => return None,
                &Cigar::Back(_) => panic!("Back: Unsupported Cigar operation (not in SAMv1 spec)."),
                &Cigar::Pad(_) => panic!("Pad: Unsupported Cigar operation (unclear SAMv1 spec).")
            }
        }
        // if no SNV evidence is found before passing the vpos, must be
        return None;
    } else {
        panic!("unsupported variant");
    }
}


/// For an indel, calculate likelihood of ref (first) or alt (second) for a given read, i.e. Pr(Z_i | ref) and Pr(Z_i | alt).
///
/// # Idea
///
/// The MAPQ already encodes the uncertainty of the given read position. Hence, all we need to do is
/// to calculate the probability that a read is sampled from the given position under the assumption that
/// (a) there is no variant, and (b) there is the given variant.
/// We do this by multiplying the sampling probability of each read base, see `prob_read_base`.
/// This is of course a simplification: other indels or SNVs in the same region are ignored.
/// However, this is reasonable because it is (a) unlikely and (b) they will be ignored in all cases,
/// such that they simply lead to globally reduced likelihoods, which are normalized away by bayes theorem.
/// Other homo/heterozgous variants on the same haplotype as the investigated are no problem. They will affect
/// both the alt and the ref case equally.
pub fn prob_read_indel(record: &bam::Record, cigar: &[Cigar], start: u32, variant: Variant, ref_seq: &[u8]) -> (LogProb, LogProb) {
    let pos = record.pos() as u32;
    let read_seq = record.seq();
    let m = read_seq.len() as u32;
    let quals = record.qual();

    let prob_read_base = |read_pos, ref_pos| {
        let ref_base = ref_seq[ref_pos as usize];
        let read_base = read_seq[read_pos as usize];
        prob_read_base(read_base, ref_base, quals[read_pos as usize])
    };

    // indels can lead to shifts in the mapping position
    // since we don't know whether any indel comes before our indel of interest, we have to consider
    // all possible shifts due to indels.
    let total_indel_len = cigar.iter().map(|c| match c {
        &Cigar::Ins(l) => l,
        &Cigar::Del(l) => l,
        _ => 0
    }).sum();

    // calculate maximal shift to the left without getting outside of the indel start
    let pos_min = cmp::max(pos.saturating_sub(total_indel_len), start.saturating_sub(m));
    // exclusive upper bound
    let pos_max = pos + total_indel_len + 1;
    debug!("--------------");
    debug!("cigar: {:?}", cigar);
    debug!("calculating indel likelihood for shifts within {} - {}", pos_min, pos_max);
    assert!(pos >= pos_min && pos <= pos_max, "original mapping position should be within the evaluated shifts");

    let capacity = (pos_max - pos_min) as usize;
    let mut prob_alts = Vec::with_capacity(capacity);
    let mut prob_refs = Vec::with_capacity(capacity);

    // debugging
    let mut alt_matches = None;
    let mut ref_matches = None;
    if log_enabled!(Debug) {
        alt_matches = Some(Vec::with_capacity(m as usize));
        ref_matches = Some(Vec::with_capacity(m as usize));
    }
    let debug_match = |read_pos, ref_pos| {
        if ref_seq[ref_pos as usize].to_ascii_uppercase() == read_seq[read_pos as usize].to_ascii_uppercase() {
            '='
        } else {
            read_seq[read_pos as usize] as char
        }
    };

    for p in pos_min..pos_max {
        let mut prob_ref = LogProb::ln_one();
        let mut prob_alt = LogProb::ln_one();

        // TODO consider the case that a deletion is made invisible via an insertion of the same size and vice versa?

        let prefix_end = start.saturating_sub(p);
        // common prefix
        for i in 0..prefix_end {
            let prob = prob_read_base(i, p + i);
            prob_ref = prob_ref + prob;
            prob_alt = prob_alt + prob;

            if log_enabled!(Debug) {
                // debugging
                let x = debug_match(i, p + i);
                alt_matches.as_mut().unwrap().push(x);
                ref_matches.as_mut().unwrap().push(x);
            }
        }

        if log_enabled!(Debug) {
            alt_matches.as_mut().unwrap().push('|');
            ref_matches.as_mut().unwrap().push('|');
            ref_matches.as_mut().unwrap().push('|');
        }

        // ref likelihood
        for i in prefix_end..m {
            let prob = prob_read_base(i, p + i);
            prob_ref = prob_ref + prob;

            if log_enabled!(Debug) {
                let x = debug_match(i, p + i);
                ref_matches.as_mut().unwrap().push(x);
            }
        }

        // alt likelihood
        match variant {
            Variant::Insertion(l) => {
                // reduce length if insertion is left of p
                let l = if start >= p { l as u32 } else { l - (p - start) };
                // TODO this will miss one base in some cases.
                // but it is needed because some callers specify calls as G->GAA and some as *->AA
                // a better place to fix is when parsing the vcf file.
                let suffix_start = (start + l + 1).saturating_sub(p);

                if log_enabled!(Debug) {
                    if suffix_start <= m {
                        for i in prefix_end..suffix_start {
                            alt_matches.as_mut().unwrap().push(read_seq[i as usize] as char);
                        }
                        alt_matches.as_mut().unwrap().push('|');
                    }
                }

                for i in suffix_start..m {
                    let prob = prob_read_base(i, p + i - l);
                    prob_alt = prob_alt + prob;

                    if log_enabled!(Debug) {
                        let x = debug_match(i, p + i - l);
                        alt_matches.as_mut().unwrap().push(x);
                    }
                }
            },
            Variant::Deletion(l) => {
                // reduce length if deletion is left of p
                let l = if start >= p { l as u32 } else { l - (p - start) };
                // TODO this will miss one base in some cases.
                // but it is needed because some callers specify calls as GAA->G and some as AA->*
                // a better place to fix is when parsing the vcf file.
                let suffix_start = (start + 1).saturating_sub(p);

                for i in suffix_start..m {
                    let prob = prob_read_base(i, p + i + l);
                    prob_alt = prob_alt + prob;

                    if log_enabled!(Debug) {
                        let x = debug_match(i, p + i + l);
                        alt_matches.as_mut().unwrap().push(x);
                    }
                }
            },
            _ => panic!("unsupported variant type")
        }
        prob_alts.push(prob_alt);
        prob_refs.push(prob_ref);

        if log_enabled!(Debug) {
            debug!("shift {}:", p as i32 - pos as i32);
            debug!("ref: {:?}", itertools::join(ref_matches.as_ref().unwrap(), ""));
            debug!("alt: {:?}", itertools::join(alt_matches.as_ref().unwrap(), ""));
            ref_matches.as_mut().unwrap().clear();
            alt_matches.as_mut().unwrap().clear();
        }
    }
    let prob_alt = LogProb::ln_sum_exp(&prob_alts);
    let prob_ref = LogProb::ln_sum_exp(&prob_refs);

    // normalize over all events, because we assume that the read comes from here (w_i=1)
    // the uncertainty of the mapping position (w_i) is captured by prob_mapping.
    let total = prob_ref.ln_add_exp(prob_alt);
    debug!("Pr(alt)={}, Pr(ref)={}, Pr(total)={}", *prob_alt, *prob_ref, *total);
    (prob_ref - total, prob_alt - total)
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
        if let Evidence::Alignment = self.evidence {
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
    Alignment
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
    max_indel_dist: u32,
    max_indel_len_diff: u32
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
    /// * `prob_missed_insertion_alignment` - rate of missed insertion alignments if insertion is present (mapper dependent, BWA: 0.2138, LASER: 0.3460)
    /// * `prob_missed_deletion_alignment` - rate of missed deletion alignments if deletion is present (mapper dependent, BWA: 0.0310, LASER: 0.0964)
    /// * `prob_spurious_indel_alignment` - rate of wrongly reported indel alignments if indel is absent (mapper dependent)
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
    ) -> Self {
        Sample {
            record_buffer: RecordBuffer::new(bam, pileup_window, use_secondary),
            use_fragment_evidence: use_fragment_evidence,
            use_mapq: use_mapq,
            adjust_mapq: adjust_mapq,
            insert_size: insert_size,
            likelihood_model: likelihood_model,
            prob_spurious_isize: LogProb::from(prob_spurious_isize),
            max_indel_dist: 50,
            max_indel_len_diff: 20
        }
    }

    pub fn max_indel_dist(mut self, dist: u32) -> Self {
        self.max_indel_dist = dist;
        self
    }

    pub fn max_indel_len_diff(mut self, diff: u32) -> Self {
        self.max_indel_len_diff = diff;
        self
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
        variant: Variant,
        chrom_seq: &[u8]
    ) -> Result<Vec<Observation>, Box<Error>> {
        let mut observations = Vec::new();
        let (end, varpos) = match variant {
            Variant::Deletion(length)  => (start + length, (start + length / 2) as i32), // TODO do we really need two centerpoints?
            Variant::Insertion(length) => (start + length, start as i32),
            Variant::SNV(_) => (start + 1, start as i32)
        };
        let mut pairs = HashMap::new();
        let mut n_overlap = 0;

        debug!("variant: {}:{} {:?}", str::from_utf8(chrom).unwrap(), start, variant);

        // move window to the current variant
        debug!("Filling buffer...");
        try!(self.record_buffer.fill(chrom, start, end));
        debug!("Done.");

        // iterate over records
        for record in self.record_buffer.iter() {
            let pos = record.pos();
            let cigar = record.cigar();
            let end_pos = record.end_pos(&cigar);
            // TODO do we also need to consider reads between start and end?
            if ((pos as u32) < start && (end_pos as u32) > start) || ((pos as u32) < end && (end_pos as u32) > end) {
                match variant {
                    Variant::Deletion(_) | Variant::Insertion(_) => {
                        // overlapping indel alignment
                        let obs = self.read_observation(&record, &cigar, start, variant, chrom_seq);
                        match obs {
                            Some(o) => {
                                observations.push(o);
                                n_overlap += 1;
                            },
                            None => ()
                        }
                    },
                    Variant::SNV(_) => {
                        let obs = self.read_observation(&record, &cigar, start, variant, chrom_seq);
                        match obs {
                            Some(o) => {
                                observations.push(o);
                                n_overlap += 1;
                            },
                            None => {
                                debug!("Did not add record to buffer, SNV position deleted or skipped.");
                            }
                        }
                    }
                }
            } else if variant.has_fragment_evidence() &&
                      self.use_fragment_evidence &&
                      (record.is_first_in_template() || record.is_last_in_template())
            {
                if end_pos <= varpos {
                    // need to check mate
                    // since the bam file is sorted by position, we can't see the mate first
                    if record.mpos() >= varpos {
                        pairs.insert(record.qname().to_owned(), record.mapq());
                    }
                } else if let Some(mate_mapq) = pairs.get(record.qname()) {
                    // mate already visited, and this read maps right of varpos
                    observations.push(self.fragment_observation(&record, *mate_mapq, variant));
                }
            }
        }
        debug!("Extracted observations ({} fragments, {} overlapping reads).", pairs.len(), n_overlap);

        if self.adjust_mapq {
            match variant {
                // only adjust for deletion and insertion
                Variant::Deletion(_) | Variant::Insertion(_) => adjust_mapq(&mut observations),
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
        cigar: &[Cigar],
        start: u32,
        variant: Variant,
        chrom_seq: &[u8]
    ) -> Option<Observation> {
        let prob_mapping = self.prob_mapping(record.mapq());
        let probs: Option<(LogProb, LogProb)>;

        match variant {
            Variant::Deletion(_) | Variant::Insertion(_) => {
                probs = Some(prob_read_indel(record, cigar, start, variant, chrom_seq));
            },
            Variant::SNV(_) => {
                probs = prob_read_snv(record, cigar, start, variant, chrom_seq);
            }
        };

        match probs {
            Some( (prob_ref, prob_alt) ) => {
                Some ( Observation {
                    prob_mapping: prob_mapping,
                    prob_alt: prob_alt,
                    prob_ref: prob_ref,
                    prob_mismapped: LogProb::ln_one(), // if the read is mismapped, we assume sampling probability 1.0
                    evidence: Evidence::Alignment
                } )
            },
            None => None
        }
    }

    fn fragment_observation(
        &self,
        record: &bam::Record,
        mate_mapq: u8,
        variant: Variant
    ) -> Observation {
        let insert_size = record.insert_size().abs();
        let shift = match variant {
            Variant::Deletion(length)  => length as f64,
            Variant::Insertion(length) => -(length as f64),
            Variant::SNV(_) => panic!("no fragment observations for SNV")
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

    use csv;
    use std::str;
    use itertools::Itertools;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use rust_htslib::bam::record::Cigar;
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
                evidence: Evidence::Alignment
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
                evidence: Evidence::Alignment
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
                evidence: Evidence::Alignment
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
            bam::IndexedReader::new(&"tests/indels.bam").unwrap(),
            2500,
            true,
            true,
            true,
            false,
            InsertSize { mean: isize_mean, sd: 20.0 },
            likelihood::LatentVariableModel::new(1.0),
            Prob(0.0)
        )
    }

    #[test]
    fn test_read_observation_indel() {
        let variant = model::Variant::Insertion(10);
        // insertion starts at 546 and has length 10
        let varpos = 546;

        let sample = setup_sample(150.0);
        let bam = bam::Reader::new(&"tests/indels.bam").unwrap();
        let records = bam.records().collect_vec();

        let ref_seq = ref_seq();

        let true_alt_probs = [0.0, 0.0, -483.4099, 0.0, -483.4099];

        for (record, true_alt_prob) in records.into_iter().zip(true_alt_probs.into_iter()) {
            let record = record.unwrap();
            let cigar = record.cigar();
            let obs = sample.read_observation(&record, &cigar, varpos, variant, &ref_seq);
            match obs {
                Some(o) => {
                    assert_relative_eq!(*o.prob_alt, *true_alt_prob, epsilon=0.001);
                    assert_relative_eq!(*o.prob_mapping, *(LogProb::from(PHREDProb(60.0)).ln_one_minus_exp()));
                    assert_relative_eq!(*o.prob_mismapped, *LogProb::ln_one());
                },
                None => panic!("read_observation() in test_read_observation_indel() returned 'None'")
            }
        }
    }

    #[test]
    fn test_fragment_observation_no_evidence() {
        let sample = setup_sample(150.0);
        let bam = bam::Reader::new(&"tests/indels.bam").unwrap();
        let records = bam.records().map(|rec| rec.unwrap()).collect_vec();

        for varlen in &[0, 5, 10, 100] {
            println!("varlen {}", varlen);
            println!("insertion");
            let variant = model::Variant::Insertion(*varlen);
            for record in &records {
                let obs = sample.fragment_observation(record, 60u8, variant);
                println!("{:?}", obs);
                if *varlen == 0 {
                    assert_relative_eq!(*obs.prob_ref, *obs.prob_alt);
                } else {
                    assert!(obs.prob_ref > obs.prob_alt);
                }
            }
            println!("deletion");
            let variant = model::Variant::Deletion(*varlen);
            for record in &records {
                let obs = sample.fragment_observation(record, 60u8, variant);
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
        let bam = bam::Reader::new(&"tests/indels.bam").unwrap();
        let records = bam.records().map(|rec| rec.unwrap()).collect_vec();

        println!("deletion");
        let sample = setup_sample(100.0);
        let variant = model::Variant::Deletion(50);
        for record in &records {
            let obs = sample.fragment_observation(record, 60u8, variant);
            println!("{:?}", obs);
            assert_relative_eq!(obs.prob_ref.exp(), 0.0, epsilon=0.001);
            assert!(obs.prob_alt > obs.prob_ref);
        }

        println!("insertion");
        let sample = setup_sample(200.0);
        let variant = model::Variant::Insertion(50);
        for record in &records {
            let obs = sample.fragment_observation(record, 60u8, variant);
            println!("{:?}", obs);
            assert_relative_eq!(obs.prob_ref.exp(), 0.0, epsilon=0.001);
            assert!(obs.prob_alt > obs.prob_ref);
        }
    }

    #[test]
    fn test_record_buffer() {
        let bam = bam::IndexedReader::new(&"tests/indels.bam").unwrap();
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

        let bam = bam::Reader::new(&"tests/indels+clips.bam").unwrap();
        let records = bam.records().map(|rec| rec.unwrap()).collect_vec();
        let ref_seq = ref_seq();

        // truth
        let probs_alt = [0.0, 0.0, -483.4099, 0.0, 0.0, 0.0];
        let probs_ref = [-250.198, -214.7137, 0.0, -250.198, -621.1132, -250.198];

        // variant (obtained via bcftools)
        let start = 546;
        let variant = model::Variant::Insertion(10);
        for (i, rec) in records.iter().enumerate() {
            println!("{}", str::from_utf8(rec.qname()).unwrap());
            let (prob_ref, prob_alt) = prob_read_indel(rec, &rec.cigar(), start, variant, &ref_seq);
            println!("Pr(ref)={} Pr(alt)={}", *prob_ref, *prob_alt);
            println!("{:?}", rec.cigar());
            assert_relative_eq!(*prob_ref, probs_ref[i], epsilon=0.001);
            assert_relative_eq!(*prob_alt, probs_alt[i], epsilon=0.001);
        }
    }

    #[test]
    fn test_prob_read_snv() {
        let ref_seq: Vec<u8> = b"CCTATACGCGT"[..].to_owned();

        let mut records: Vec<bam::Record> = Vec::new();
        let mut qname: &[u8];
        let mut seq: &[u8];

        // Ignore leading HardClip, skip leading SoftClip, reference nucleotide
        qname = b"HC_SC_M";
        let cigar = [Cigar::HardClip(5), Cigar::SoftClip(2), Cigar::Match(6)];
        seq  = b"AATATACG";
        let qual = [20, 20, 30, 30, 30, 40, 30, 30];
        let mut record1 = bam::Record::new();
        record1.set(qname, &cigar, seq, &qual);
        record1.set_pos(2);
        records.push(record1);

        // Ignore leading HardClip, skip leading Insertion, alternative nucleotide
        qname = b"HC_Ins_M";
        let cigar = [Cigar::HardClip(2), Cigar::Ins(2), Cigar::Match(6)];
        seq  = b"TTTATGCG";
        let qual = [20, 20, 20, 20, 20, 30, 20, 20];
        let mut record2 = bam::Record::new();
        record2.set(qname, &cigar, seq, &qual);
        record2.set_pos(2);
        records.push(record2);

        // Matches and deletion before position, reference nucleotide
        qname = b"Eq_Diff_Del_Eq";
        let cigar = [Cigar::Equal(2), Cigar::Diff(1), Cigar::Del(2), Cigar::Equal(5)];
        seq  = b"CCAACGCG";
        let qual = [30, 30, 30, 50, 30, 30, 30, 30];
        let mut record3 = bam::Record::new();
        record3.set(qname, &cigar, seq, &qual);
        record3.set_pos(0);
        records.push(record3);

        // single nucleotide Deletion covering SNV position
        qname = b"M_Del_M";
        let cigar = [Cigar::Match(4), Cigar::Del(1), Cigar::Match(4)];
        seq  = b"CTATCGCG";
        let qual = [10, 30, 30, 30, 30, 30, 30, 30];
        let mut record4 = bam::Record::new();
        record4.set(qname, &cigar, seq, &qual);
        record4.set_pos(1);
        records.push(record4);

        // three nucleotide RefSkip covering SNV position
        qname = b"M_RefSkip_M";
        let cigar = [Cigar::Equal(1), Cigar::Diff(1), Cigar::Equal(2), Cigar::RefSkip(3), Cigar::Match(4)];
        seq  = b"CTTAGCGT";
        let qual = [10, 30, 30, 30, 30, 30, 30, 30];
        let mut record5 = bam::Record::new();
        record5.set(qname, &cigar, seq, &qual);
        record5.set_pos(0);
        records.push(record5);


        // truth
        let probs_ref = [0.9999,   0.00033, 0.99999  ];
        let probs_alt = [0.000033, 0.999,   0.0000033];
        let eps       = [0.000001, 0.00001, 0.0000001];

        let vpos = 5;
        let variant = model::Variant::SNV(b'G');
        for (i, rec) in records.iter().enumerate() {
            println!("{}", str::from_utf8(rec.qname()).unwrap());
            let res = prob_read_snv(rec, &rec.cigar(), vpos, variant, &ref_seq);
            match res {
                Some( (prob_ref, prob_alt) ) => {
                    println!("{:?}", rec.cigar());
                    println!("Pr(ref)={} Pr(alt)={}", (*prob_ref).exp(), (*prob_alt).exp() );
                    assert_relative_eq!( (*prob_ref).exp(), probs_ref[i], epsilon = eps[i]);
                    assert_relative_eq!( (*prob_alt).exp(), probs_alt[i], epsilon = eps[i]);
                },
                None => {
                    // anything that's tested for the reference position not being covered, should
                    // have 10 as the quality value of the first base in seq
                    assert_eq!(rec.qual()[0], 10);
                }
            }
        }
    }
}

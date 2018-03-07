use std::cmp;
use std::str;
use std::error::Error;
use std::ascii::AsciiExt;
use std::str::FromStr;

use itertools::Itertools;
use regex::Regex;

use bio::stats::{LogProb, PHREDProb, Prob};
use rust_htslib::bam::record::{CigarStringView, Cigar, CigarString};
use rust_htslib::bam;

use model::Variant;
use model::evidence::observation::ProbSampleAlt;
use pairhmm;



pub fn prob_snv(
    record: &bam::Record,
    cigar: &CigarStringView,
    start: u32,
    variant: &Variant,
    ref_seq: &[u8]
) -> Result<Option<(LogProb, LogProb)>, Box<Error>> {
    if let &Variant::SNV(base) = variant {
        if let Some(qpos) = cigar.read_pos(start, false, false)? {
            let read_base = record.seq()[qpos as usize];
            let base_qual = record.qual()[qpos as usize];
            let prob_alt = prob_read_base(read_base, base, base_qual);
            let prob_ref = prob_read_base(read_base, ref_seq[start as usize], base_qual);
            Ok( Some( (prob_ref, prob_alt) ) )
        } else {
            // a read that spans an SNV might have the respective position deleted (Cigar op 'D')
            // or reference skipped (Cigar op 'N'), and the library should not choke on those reads
            // but instead needs to know NOT to add those reads (as observations) further up
            Ok( None )
        }
    } else {
        panic!("bug: unsupported variant");
    }
}


/// Calculate read evindence for an indel.
pub struct IndelEvidence {
    gap_params: IndelGapParams,
    pairhmm: pairhmm::PairHMM,
    window: u32
}


impl IndelEvidence {
    /// Create a new instance.
    pub fn new(
        prob_insertion_artifact: LogProb,
        prob_deletion_artifact: LogProb,
        prob_insertion_extend_artifact: LogProb,
        prob_deletion_extend_artifact: LogProb,
        window: u32
    ) -> Self {
        IndelEvidence {
            gap_params: IndelGapParams {
                prob_insertion_artifact: prob_insertion_artifact,
                prob_deletion_artifact: prob_deletion_artifact,
                prob_insertion_extend_artifact: prob_insertion_extend_artifact,
                prob_deletion_extend_artifact: prob_deletion_extend_artifact
            },
            pairhmm: pairhmm::PairHMM::new(),
            window: window
        }
    }

    /// Calculate probability for reference and alternative allele.
    pub fn prob(&mut self,
        record: &bam::Record,
        cigar: &CigarStringView,
        start: u32,
        variant: &Variant,
        ref_seq: &[u8]
    ) -> Result<(LogProb, LogProb), Box<Error>> {
        let read_seq = record.seq();
        let read_qual = record.qual();

        let (read_offset, read_end, breakpoint, overlap) = {
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
                    let read_offset = qstart.saturating_sub(self.window as usize);
                    let read_end = cmp::min(
                        qend + self.window as usize,
                        read_seq.len()
                    );
                    (read_offset, read_end, varstart as usize, true)
                },
                (Some(qstart), None) => {
                    let qstart = qstart as usize;
                    let read_offset = qstart.saturating_sub(self.window as usize);
                    let read_end = cmp::min(
                        qstart + self.window as usize,
                        read_seq.len()
                    );
                    (read_offset, read_end, varstart as usize, true)
                },
                (None, Some(qend)) => {
                    let qend = qend as usize;
                    let read_offset = qend.saturating_sub(self.window as usize);
                    let read_end = cmp::min(
                        qend + self.window as usize,
                        read_seq.len()
                    );
                    (read_offset, read_end, varend as usize, true)
                },
                (None, None) => {
                    let m = read_seq.len() / 2;
                    let read_offset = m.saturating_sub(self.window as usize);
                    let read_end = cmp::min(m + self.window as usize, read_seq.len());
                    let breakpoint = record.pos() as usize + m;
                    (read_offset, read_end, breakpoint, false)
                }
            }
        };

        let start = start as usize;
        // the window on the reference should be a bit larger to allow some flexibility with close
        // indels. But it should not be so large that the read can align outside of the breakpoint.
        let ref_window = (self.window as f64 * 1.5) as usize;

        // ref allele
        let prob_ref = self.pairhmm.prob_related(
            &self.gap_params,
            &ReferenceEmissionParams {
                ref_seq: ref_seq,
                read_seq: &read_seq,
                read_qual: read_qual,
                read_offset: read_offset,
                read_end: read_end,
                ref_offset: breakpoint.saturating_sub(ref_window),
                ref_end: cmp::min(breakpoint + ref_window, ref_seq.len()),
            }
        );

        // alt allele
        let prob_alt = if overlap {
            match variant {
                &Variant::Deletion(_) => {
                    self.pairhmm.prob_related(
                        &self.gap_params,
                        &DeletionEmissionParams {
                            ref_seq: ref_seq,
                            read_seq: &read_seq,
                            read_qual: read_qual,
                            read_offset: read_offset,
                            read_end: read_end,
                            ref_offset: start.saturating_sub(ref_window),
                            ref_end: cmp::min(start + ref_window, ref_seq.len()),
                            del_start: start,
                            del_len: variant.len() as usize
                        }
                    )
                },
                &Variant::Insertion(ref ins_seq) => {
                    let l = ins_seq.len() as usize;
                    self.pairhmm.prob_related(
                        &self.gap_params,
                        &InsertionEmissionParams {
                            ref_seq: ref_seq,
                            read_seq: &read_seq,
                            read_qual: read_qual,
                            read_offset: read_offset,
                            read_end: read_end,
                            ref_offset: start.saturating_sub(ref_window),
                            ref_end: cmp::min(start + l + ref_window, ref_seq.len()),
                            ins_start: start,
                            ins_len: l,
                            ins_end: start + l,
                            ins_seq: ins_seq
                        }
                    )
                },
                _ => {
                    panic!("bug: unsupported variant");
                }
            }
        } else {
            // if no overlap, we can simply use prob_ref again
            prob_ref
        };

        Ok((prob_ref, prob_alt))
    }

    /// Probability to sample read from alt allele for each possible max softclip up to a given
    /// theoretical maximum.
    /// If variant is small enough to be in CIGAR, max_softclip should be set to None
    /// (i.e., ignored), and the method will only return one value.
    ///
    /// The key idea is calculate the probability as number of valid placements (considering the
    /// max softclip allowed by the mapper) over all possible placements.
    pub fn prob_sample_alt(
        &self,
        read_len: u32,
        enclosing_possible: bool,
        variant: &Variant
    ) -> ProbSampleAlt {
        let delta = match variant {
            &Variant::Deletion(_)  => variant.len() as u32,
            &Variant::Insertion(_) => variant.len() as u32,
            &Variant::SNV(_) => return ProbSampleAlt::One
        };

        let prob = |max_softclip| {
            let n_alt = cmp::min(delta, read_len);
            let n_alt_valid = cmp::min(n_alt, max_softclip);

            LogProb((n_alt_valid as f64).ln() - (n_alt as f64).ln())
        };

        if !enclosing_possible {
            ProbSampleAlt::Dependent((0..read_len + 1).map(&prob).collect_vec())
        } else {
            ProbSampleAlt::Independent(prob(read_len))
        }
    }
}


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


/// Convert MAPQ (from read mapper) to LogProb for the event that the read maps correctly.
pub fn prob_mapping(record: &bam::Record) -> LogProb {
    LogProb::from(PHREDProb(record.mapq() as f64)).ln_one_minus_exp()
}


pub fn prob_mapping_adjusted(
    record: &bam::Record,
    cigar: &bam::record::CigarStringView,
    chrom_name: &[u8],
    chrom_seq: &[u8]
) -> Result<LogProb, Box<Error>> {
    fn likelihood(
        record: &bam::Record,
        cigar: &bam::record::CigarStringView,
        pos: u32,
        chrom_seq: &[u8]
    ) -> LogProb {
        let seq = record.seq();
        let qual = record.qual();
        let mut ref_pos = pos as u32;
        let mut read_pos = 0;
        let mut lh = LogProb::ln_one();
        for c in cigar {
            match c {
                &Cigar::Match(n) |
                &Cigar::Diff(n)  |
                &Cigar::Equal(n) => {
                    for _ in 0..n {
                        lh += prob_read_base(
                            seq[read_pos as usize],
                            chrom_seq[ref_pos as usize],
                            qual[read_pos as usize]
                        );
                        ref_pos += 1;
                        read_pos += 1;
                    }
                },
                &Cigar::Ins(l) => {
                    read_pos += l;
                },
                &Cigar::Del(l) => {
                    ref_pos += l;
                },
                &Cigar::HardClip(_) => {
                    // nothing happens because the read sequence is clipped
                },
                &Cigar::SoftClip(l) | &Cigar::Pad(l) => {
                    read_pos += l;
                },
                &Cigar::RefSkip(l) => {
                    ref_pos += l;
                }
            }
        }
        lh
    };

    let mut adjusted = false;
    if let Some(xa) = record.aux(b"XA") {
        let xa = xa.string();
        lazy_static! {
            // regex for a cigar string operation
            static ref XA_ENTRY: Regex = Regex::new(
                "(?P<chrom>[^,]+),[+-]?(?P<pos>[0-9]+),(?P<cigar>([0-9]+[MIDNSHP=X])+),[0-9]+;"
            ).unwrap();
        }

        let mut summands = Vec::new();
        for entry in XA_ENTRY.captures_iter(str::from_utf8(xa).unwrap()) {
            // sum over all XA entries on same chromosome
            if entry["chrom"].as_bytes() == chrom_name {
                // XA pos is 1-based, we need a 0-based position
                let pos = u32::from_str(&entry["pos"])? - 1;
                let xcigar = CigarString::from_str(&entry["cigar"])?;
                let cigar_view = xcigar.into_view(pos as i32);
                let lh = likelihood(record, &cigar_view, pos, chrom_seq);
                summands.push(lh);
                adjusted = true;
            }
        }
        if adjusted {
            let lh_primary = likelihood(record, cigar, record.pos() as u32, chrom_seq);
            summands.push(lh_primary);
            //println!("MAPQ: {}, {:?} vs {:?} with {}", record.mapq(), lh_primary, summands, str::from_utf8(xa).unwrap());
            let marginal = LogProb::ln_sum_exp(&summands);
            return Ok(lh_primary - marginal);
        }
    }
    // if no XA tag on same chromosome, use MAPQ given by mapper.
    Ok(prob_mapping(record))
}


/// Gap parameters for PairHMM.
pub struct IndelGapParams {
    pub prob_insertion_artifact: LogProb,
    pub prob_deletion_artifact: LogProb,
    pub prob_insertion_extend_artifact: LogProb,
    pub prob_deletion_extend_artifact: LogProb
}


impl pairhmm::GapParameters for IndelGapParams {
    #[inline]
    fn prob_gap_x(&self) -> LogProb {
        self.prob_insertion_artifact
    }

    #[inline]
    fn prob_gap_y(&self) -> LogProb {
        self.prob_deletion_artifact
    }

    #[inline]
    fn prob_gap_x_extend(&self) -> LogProb {
        self.prob_insertion_extend_artifact
    }

    #[inline]
    fn prob_gap_y_extend(&self) -> LogProb {
        self.prob_deletion_extend_artifact
    }
}


impl pairhmm::StartEndGapParameters for IndelGapParams {
    /// Semiglobal alignment: return true.
    #[inline]
    fn free_start_gap_x(&self) -> bool {
        true
    }

    /// Semiglobal alignment: return true.
    #[inline]
    fn free_end_gap_x(&self) -> bool {
        true
    }

    /// Semiglobal alignment: return 1.0.
    #[inline]
    fn prob_start_gap_x(&self, _: usize) -> LogProb {
        LogProb::ln_one()
    }
}


macro_rules! default_emission {
    () => (
        #[inline]
        fn prob_emit_xy(&self, i: usize, j: usize) -> LogProb {
            let r = self.ref_base(i);
            let j_ = self.project_j(j);
            prob_read_base(self.read_seq[j_], r, self.read_qual[j_])
        }

        #[inline]
        fn prob_emit_x(&self, _: usize) -> LogProb {
            LogProb::ln_one()
        }

        #[inline]
        fn prob_emit_y(&self, j: usize) -> LogProb {
            prob_read_base_miscall(self.read_qual[self.project_j(j)])
        }

        #[inline]
        fn len_x(&self) -> usize {
            self.ref_end - self.ref_offset
        }

        #[inline]
        fn len_y(&self) -> usize {
            self.read_end - self.read_offset
        }
    )
}


/// Emission parameters for PairHMM over reference allele.
pub struct ReferenceEmissionParams<'a> {
    ref_seq: &'a [u8],
    read_seq: &'a bam::record::Seq<'a>,
    read_qual: &'a [u8],
    read_offset: usize,
    ref_offset: usize,
    read_end: usize,
    ref_end: usize
}


impl<'a> ReferenceEmissionParams<'a> {
    #[inline]
    fn ref_base(&self, i: usize) -> u8 {
        self.ref_seq[i + self.ref_offset]
    }

    #[inline]
    fn project_j(&self, j: usize) -> usize {
        j + self.read_offset
    }
}


impl<'a> pairhmm::EmissionParameters for ReferenceEmissionParams<'a> {
    default_emission!();
}


/// Emission parameters for PairHMM over deletion allele.
pub struct DeletionEmissionParams<'a> {
    ref_seq: &'a [u8],
    read_seq: &'a bam::record::Seq<'a>,
    read_qual: &'a [u8],
    read_offset: usize,
    ref_offset: usize,
    read_end: usize,
    ref_end: usize,
    del_start: usize,
    del_len: usize
}


impl<'a> DeletionEmissionParams<'a> {
    #[inline]
    fn ref_base(&self, i: usize) -> u8 {
        let i_ = i + self.ref_offset;
        if i_ <= self.del_start {
            self.ref_seq[i_]
        } else {
            self.ref_seq[i_ + self.del_len]
        }
    }

    #[inline]
    fn project_j(&self, j: usize) -> usize {
        j + self.read_offset
    }
}


impl<'a> pairhmm::EmissionParameters for DeletionEmissionParams<'a> {
    default_emission!();
}


/// Emission parameters for PairHMM over insertion allele.
pub struct InsertionEmissionParams<'a> {
    ref_seq: &'a [u8],
    read_seq: &'a bam::record::Seq<'a>,
    read_qual: &'a [u8],
    read_offset: usize,
    ref_offset: usize,
    read_end: usize,
    ref_end: usize,
    ins_start: usize,
    ins_end: usize,
    ins_len: usize,
    ins_seq: &'a [u8]
}


impl<'a> InsertionEmissionParams<'a> {
    #[inline]
    fn ref_base(&self, i: usize) -> u8 {
        let i_ = i + self.ref_offset;
        if i_ <= self.ins_start {
            self.ref_seq[i_]
        } else if i_ > self.ins_end {
            self.ref_seq[i_ - self.ins_len]
        } else {
            self.ins_seq[i_ - (self.ins_start + 1)]
        }
    }

    #[inline]
    fn project_j(&self, j: usize) -> usize {
        j + self.read_offset
    }
}


impl<'a> pairhmm::EmissionParameters for InsertionEmissionParams<'a> {
    default_emission!();
}

#[cfg(test)]
mod tests {

    use super::*;
    use model;

    use std::str;
    use rust_htslib::bam::record::{Cigar, CigarString};

    #[test]
    fn test_prob_snv() {
        let ref_seq: Vec<u8> = b"CCTATACGCGT"[..].to_owned();

        let mut records: Vec<bam::Record> = Vec::new();
        let mut qname: &[u8];
        let mut seq: &[u8];

        // Ignore leading HardClip, skip leading SoftClip, reference nucleotide
        qname = b"HC_SC_M";
        let cigar = CigarString( vec![Cigar::HardClip(5), Cigar::SoftClip(2), Cigar::Match(6)] );
        seq  = b"AATATACG";
        let qual = [20, 20, 30, 30, 30, 40, 30, 30];
        let mut record1 = bam::Record::new();
        record1.set(qname, &cigar, seq, &qual);
        record1.set_pos(2);
        records.push(record1);

        // Ignore leading HardClip, skip leading Insertion, alternative nucleotide
        qname = b"HC_Ins_M";
        let cigar = CigarString( vec![Cigar::HardClip(2), Cigar::Ins(2), Cigar::Match(6)] );
        seq  = b"TTTATGCG";
        let qual = [20, 20, 20, 20, 20, 30, 20, 20];
        let mut record2 = bam::Record::new();
        record2.set(qname, &cigar, seq, &qual);
        record2.set_pos(2);
        records.push(record2);

        // Matches and deletion before position, reference nucleotide
        qname = b"Eq_Diff_Del_Eq";
        let cigar = CigarString( vec![Cigar::Equal(2), Cigar::Diff(1), Cigar::Del(2), Cigar::Equal(5)] );
        seq  = b"CCAACGCG";
        let qual = [30, 30, 30, 50, 30, 30, 30, 30];
        let mut record3 = bam::Record::new();
        record3.set(qname, &cigar, seq, &qual);
        record3.set_pos(0);
        records.push(record3);

        // single nucleotide Deletion covering SNV position
        qname = b"M_Del_M";
        let cigar = CigarString( vec![Cigar::Match(4), Cigar::Del(1), Cigar::Match(4)] );
        seq  = b"CTATCGCG";
        let qual = [10, 30, 30, 30, 30, 30, 30, 30];
        let mut record4 = bam::Record::new();
        record4.set(qname, &cigar, seq, &qual);
        record4.set_pos(1);
        records.push(record4);

        // three nucleotide RefSkip covering SNV position
        qname = b"M_RefSkip_M";
        let cigar = CigarString( vec![Cigar::Equal(1), Cigar::Diff(1), Cigar::Equal(2), Cigar::RefSkip(3), Cigar::Match(4)] );
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
            if let Ok( Some( (prob_ref, prob_alt) ) ) = prob_snv(rec, &rec.cigar(), vpos, &variant, &ref_seq) {
                println!("{:?}", rec.cigar());
                println!("Pr(ref)={} Pr(alt)={}", (*prob_ref).exp(), (*prob_alt).exp() );
                assert_relative_eq!( (*prob_ref).exp(), probs_ref[i], epsilon = eps[i]);
                assert_relative_eq!( (*prob_alt).exp(), probs_alt[i], epsilon = eps[i]);
            } else {
                // anything that's tested for the reference position not being covered, should
                // have 10 as the quality value of the first base in seq
                assert_eq!(rec.qual()[0], 10);
            }
        }
    }
}

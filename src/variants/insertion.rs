use std::cell::RefCell;
use std::cmp;
use std::rc::Rc;
use std::sync::Arc;

use anyhow::Result;
use bio::stats::pairhmm::EmissionParameters;
use bio::stats::LogProb;
use bio_types::genome::{self, AbstractInterval, AbstractLocus};

use crate::estimation::alignment_properties::AlignmentProperties;
use crate::variants::evidence::realignment::pairhmm::{ReadEmission, RefBaseEmission};
use crate::variants::evidence::realignment::{Realignable, Realigner};
use crate::variants::sampling_bias::{ReadSamplingBias, SamplingBias};
use crate::variants::{AlleleProb, MultiLocus, PairedEndEvidence, SingleLocus, Variant};
use crate::{default_emission, default_ref_base_emission};

pub struct Insertion {
    locus: MultiLocus,
    ins_seq: Rc<Vec<u8>>,
    realigner: RefCell<Realigner>,
}

impl Insertion {
    pub fn new(locus: genome::Locus, ins_seq: Vec<u8>, realigner: Realigner) -> Self {
        Insertion {
            locus: MultiLocus::new(vec![SingleLocus::new(genome::Interval::new(
                locus.contig().to_owned(),
                locus.pos()..locus.pos() + 1,
            ))]),
            ins_seq: Rc::new(ins_seq),
            realigner: RefCell::new(realigner),
        }
    }

    pub fn locus(&self) -> &SingleLocus {
        &self.locus[0]
    }
}

impl<'a> Realignable<'a> for Insertion {
    type EmissionParams = InsertionEmissionParams<'a>;

    fn alt_emission_params(
        &self,
        read_emission_params: Rc<ReadEmission<'a>>,
        ref_seq: Arc<Vec<u8>>,
        ref_window: usize,
    ) -> InsertionEmissionParams<'a> {
        let l = self.ins_seq.len() as usize;
        let start = self.locus().range().start as usize;
        let ref_seq_len = ref_seq.len();
        InsertionEmissionParams {
            ref_seq: ref_seq,
            ref_offset: start.saturating_sub(ref_window),
            ref_end: cmp::min(start + l + ref_window, ref_seq_len),
            ins_start: start,
            ins_len: l,
            ins_end: start + l,
            ins_seq: Rc::clone(&self.ins_seq),
            read_emission: read_emission_params,
        }
    }
}

impl<'a> SamplingBias<'a> for Insertion {
    fn len(&self) -> u64 {
        self.ins_seq.len() as u64
    }
}

impl<'a> ReadSamplingBias<'a> for Insertion {}

impl<'a> Variant<'a> for Insertion {
    type Evidence = PairedEndEvidence<'a>;
    type Loci = MultiLocus;

    fn is_valid_evidence(&self, evidence: &Self::Evidence) -> Option<Vec<usize>> {
        if match evidence {
            PairedEndEvidence::SingleEnd(read) => !self.locus().overlap(read, true).is_none(),
            PairedEndEvidence::PairedEnd { left, right } => {
                !self.locus().overlap(left, true).is_none()
                    || !self.locus().overlap(right, true).is_none()
            }
        } {
            Some(vec![0])
        } else {
            None
        }
    }

    /// Return variant loci.
    fn loci(&self) -> &Self::Loci {
        &self.locus
    }

    /// Calculate probability for alt and reference allele.
    fn prob_alleles(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> Result<Option<AlleleProb>> {
        match evidence {
            PairedEndEvidence::SingleEnd(record) => Ok(Some(
                self.realigner
                    .borrow_mut()
                    .prob_alleles(record, self.locus(), self)?,
            )),
            PairedEndEvidence::PairedEnd { left, right } => {
                let prob_left =
                    self.realigner
                        .borrow_mut()
                        .prob_alleles(left, self.locus(), self)?;
                let prob_right =
                    self.realigner
                        .borrow_mut()
                        .prob_alleles(right, self.locus(), self)?;

                Ok(Some(AlleleProb::new(
                    prob_left.ref_allele() + prob_right.ref_allele(),
                    prob_left.alt_allele() + prob_right.alt_allele(),
                )))
            }
        }
    }

    fn prob_sample_alt(
        &self,
        evidence: &Self::Evidence,
        alignment_properties: &AlignmentProperties,
    ) -> LogProb {
        match evidence {
            PairedEndEvidence::PairedEnd { left, right } => {
                // METHOD: we do not require the fragment to enclose the variant.
                // Hence, we treat both reads independently.
                (self
                    .prob_sample_alt_read(left.seq().len() as u64, alignment_properties)
                    .ln_one_minus_exp()
                    + self
                        .prob_sample_alt_read(right.seq().len() as u64, alignment_properties)
                        .ln_one_minus_exp())
                .ln_one_minus_exp()
            }
            PairedEndEvidence::SingleEnd(read) => {
                self.prob_sample_alt_read(read.seq().len() as u64, alignment_properties)
            }
        }
    }
}

/// Emission parameters for PairHMM over insertion allele.
#[derive(Debug)]
pub struct InsertionEmissionParams<'a> {
    ref_seq: Arc<Vec<u8>>,
    ref_offset: usize,
    ref_end: usize,
    ins_start: usize,
    ins_end: usize,
    ins_len: usize,
    ins_seq: Rc<Vec<u8>>,
    read_emission: Rc<ReadEmission<'a>>,
}

impl<'a> RefBaseEmission for InsertionEmissionParams<'a> {
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

    default_ref_base_emission!();
}

impl<'a> EmissionParameters for InsertionEmissionParams<'a> {
    default_emission!();

    #[inline]
    fn len_x(&self) -> usize {
        self.ref_end - self.ref_offset + self.ins_len
    }
}

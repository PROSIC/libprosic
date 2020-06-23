// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp::Ordering;
use std::fmt::Debug;
use std::ops::{Deref, Range};
use std::str;

use bio::stats::LogProb;
use itertools::Itertools;
use ordered_float::NotNan;
use strum::IntoEnumIterator;
use strum_macros::{EnumIter, EnumString, IntoStaticStr};

use crate::grammar;

pub(crate) mod likelihood;
pub(crate) mod modes;

#[derive(Debug, Clone)]
pub(crate) struct Contamination {
    pub(crate) by: usize,
    pub(crate) fraction: f64,
}

#[derive(Ord, Eq, PartialOrd, PartialEq, Clone, Debug)]
pub(crate) struct Event {
    pub(crate) name: String,
    pub(crate) vafs: grammar::VAFTree,
    pub(crate) strand_bias: StrandBias,
}

impl Event {
    pub(crate) fn is_artifact(&self) -> bool {
        self.strand_bias != StrandBias::None
    }
}

pub(crate) type AlleleFreq = NotNan<f64>;

#[derive(Copy, Clone, PartialOrd, PartialEq, Eq, Debug, Ord)]
pub(crate) enum StrandBias {
    None,
    Forward,
    Reverse,
}

impl Default for StrandBias {
    fn default() -> Self {
        StrandBias::None
    }
}

impl StrandBias {
    pub(crate) fn is_some(&self) -> bool {
        if let StrandBias::None = self {
            false
        } else {
            true
        }
    }

    pub(crate) fn forward_rate(&self) -> LogProb {
        match self {
            StrandBias::None => LogProb(0.5_f64.ln()),
            StrandBias::Forward => LogProb::ln_one(),
            StrandBias::Reverse => LogProb::ln_zero(),
        }
    }

    pub(crate) fn reverse_rate(&self) -> LogProb {
        match self {
            StrandBias::None => LogProb(0.5_f64.ln()),
            StrandBias::Forward => LogProb::ln_zero(),
            StrandBias::Reverse => LogProb::ln_one(),
        }
    }
}

#[allow(non_snake_case)]
pub(crate) fn AlleleFreq(af: f64) -> AlleleFreq {
    NotNan::new(af).unwrap()
}

pub(crate) trait AlleleFreqs: Debug {
    fn is_absent(&self) -> bool;
}
impl AlleleFreqs for DiscreteAlleleFreqs {
    fn is_absent(&self) -> bool {
        self.inner.len() == 1 && self.inner[0] == AlleleFreq(0.0)
    }
}
impl AlleleFreqs for ContinuousAlleleFreqs {
    fn is_absent(&self) -> bool {
        self.is_singleton() && self.start == AlleleFreq(0.0)
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct DiscreteAlleleFreqs {
    inner: Vec<AlleleFreq>,
}

impl DiscreteAlleleFreqs {
    /// Create spectrum of discrete allele frequencies with given values.
    pub(crate) fn new(spectrum: Vec<AlleleFreq>) -> Self {
        DiscreteAlleleFreqs { inner: spectrum }
    }

    /// Return spectrum of all feasible allele frequencies for a given ploidy and maximum
    /// number of amplification.
    ///
    /// In principle, a ploidy of, e.g., 2 allows allele frequencies 0, 0.5, 1.0. However,
    /// if a locus is amplified e.g. 1 time, allele frequencies can be effectively
    /// 0, 0.25, 0.5, 0.75, 1.0, because all reads are projected to the original location in the
    /// reference genome, and it is unclear whether the first, the second, or both loci contain
    /// the variant.
    ///
    /// # Arguments
    /// * ploidy - the assumed overall ploidy
    /// * max_amplification - the maximum amplification factor (1 means no amplification, 2 means
    ///   at most one duplicate, ...).
    fn _feasible(ploidy: u32, max_amplification: u32) -> Self {
        let n = ploidy * max_amplification;
        DiscreteAlleleFreqs {
            inner: (0..n + 1)
                .map(|m| AlleleFreq(m as f64 / n as f64))
                .collect_vec(),
        }
    }

    /// Return spectrum of possible allele frequencies given a ploidy.
    pub(crate) fn feasible(ploidy: u32) -> Self {
        Self::_feasible(ploidy, 1)
    }

    /// Return all frequencies except 0.0.
    pub(crate) fn not_absent(&self) -> Self {
        DiscreteAlleleFreqs {
            inner: self.inner[1..].to_owned(),
        }
    }

    pub(crate) fn absent() -> Self {
        DiscreteAlleleFreqs {
            inner: vec![AlleleFreq(0.0)],
        }
    }
}

impl Deref for DiscreteAlleleFreqs {
    type Target = Vec<AlleleFreq>;

    fn deref(&self) -> &Vec<AlleleFreq> {
        &self.inner
    }
}

/// An allele frequency range
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ContinuousAlleleFreqs {
    inner: Range<AlleleFreq>,
    pub(crate) left_exclusive: bool,
    pub(crate) right_exclusive: bool,
    /// offset to add when calculating the smallest observable value for a left-exclusive 0.0 bound
    zero_offset: NotNan<f64>,
}

impl ContinuousAlleleFreqs {
    pub(crate) fn absent() -> Self {
        Self::singleton(0.0)
    }

    pub(crate) fn singleton(value: f64) -> Self {
        ContinuousAlleleFreqs {
            inner: AlleleFreq(value)..AlleleFreq(value),
            left_exclusive: false,
            right_exclusive: false,
            zero_offset: NotNan::from(1.0),
        }
    }

    pub(crate) fn is_singleton(&self) -> bool {
        self.start == self.end
    }
}

impl Default for ContinuousAlleleFreqs {
    fn default() -> Self {
        Self::absent()
    }
}

impl Ord for ContinuousAlleleFreqs {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.inner.start.cmp(&other.start) {
            Ordering::Equal => self.inner.end.cmp(&other.end),
            ord @ _ => ord,
        }
    }
}

impl PartialOrd for ContinuousAlleleFreqs {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Deref for ContinuousAlleleFreqs {
    type Target = Range<AlleleFreq>;

    fn deref(&self) -> &Range<AlleleFreq> {
        &self.inner
    }
}

#[derive(
    Display,
    Debug,
    Clone,
    Serialize,
    Deserialize,
    EnumString,
    EnumIter,
    IntoStaticStr,
    EnumVariantNames,
)]
pub enum VariantType {
    #[strum(serialize = "INS")]
    Insertion(Option<Range<u64>>),
    #[strum(serialize = "DEL")]
    Deletion(Option<Range<u64>>),
    #[strum(serialize = "SNV")]
    SNV,
    #[strum(serialize = "MNV")]
    MNV,
    #[strum(serialize = "REF")]
    None, // site with no suggested alternative allele
}

impl From<&str> for VariantType {
    fn from(string: &str) -> VariantType {
        match string {
            "INS" => VariantType::Insertion(None),
            "DEL" => VariantType::Deletion(None),
            "SNV" => VariantType::SNV,
            "REF" => VariantType::None,
            _ => panic!("bug: given string does not describe a valid variant type"),
        }
    }
}

#[derive(Clone, Debug)]
pub(crate) enum Variant {
    Deletion(u64),
    Insertion(Vec<u8>),
    SNV(u8),
    MNV(Vec<u8>),
    None,
}

impl Variant {
    pub(crate) fn is_type(&self, vartype: &VariantType) -> bool {
        match (self, vartype) {
            (&Variant::Deletion(l), &VariantType::Deletion(Some(ref range))) => {
                l >= range.start && l < range.end
            }
            (&Variant::Insertion(_), &VariantType::Insertion(Some(ref range))) => {
                self.len() >= range.start && self.len() < range.end
            }
            (&Variant::Deletion(_), &VariantType::Deletion(None)) => true,
            (&Variant::Insertion(_), &VariantType::Insertion(None)) => true,
            (&Variant::SNV(_), &VariantType::SNV) => true,
            (&Variant::MNV(_), &VariantType::MNV) => true,
            (&Variant::None, &VariantType::None) => true,
            _ => false,
        }
    }

    pub(crate) fn len(&self) -> u64 {
        match self {
            &Variant::Deletion(l) => l,
            &Variant::Insertion(ref s) => s.len() as u64,
            &Variant::SNV(_) => 1,
            &Variant::MNV(ref alt) => alt.len() as u64,
            &Variant::None => 1,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::variants::evidence::observation::{Observation, ObservationBuilder};

    use bio::stats::LogProb;

    pub(crate) fn observation(prob_mapping: LogProb, prob_alt: LogProb, prob_ref: LogProb) -> Observation {
        ObservationBuilder::default()
            .prob_mapping_mismapping(prob_mapping)
            .prob_alt(prob_alt)
            .prob_ref(prob_ref)
            .prob_missed_allele(prob_ref.ln_add_exp(prob_alt) - LogProb(2.0_f64.ln()))
            .prob_sample_alt(LogProb::ln_one())
            .prob_overlap(LogProb::ln_one())
            .prob_any_strand(LogProb::ln_one())
            .forward_strand(true)
            .reverse_strand(true)
            .build()
            .unwrap()
    }
}

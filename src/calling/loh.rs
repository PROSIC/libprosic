// Copyright 2019 Johannes Köster, Jan Forster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::ops::RangeInclusive;
use std::iter;
use std::mem;
use std::path::Path;

use anyhow::Result;
use bio::stats::{bayesian::bayes_factors::BayesFactor, hmm, hmm::Model, LogProb, PHREDProb, Prob};
use derive_builder::Builder;
use itertools::join;
use itertools::Itertools;
use itertools_num::linspace;
use rayon::prelude::*;
use rust_htslib::bcf;
use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::{Read, HeaderView};

use crate::model::modes::tumor::TumorNormalPairView;
use crate::model::AlleleFreq;
use crate::utils;
use lp_modeler::dsl::*;
use lp_modeler::solvers::{CbcSolver, SolverTrait};
use std::cmp::Ordering;

impl Ord for RangeInclusive<u64> {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.start().cmp(&other.start() ) {
            Ordering::Less => Ordering::Less,
            Ordering::Greater => Ordering::Greater,
            Ordering::Equal => {
                self.end().cmp(&other.end() )
            },
        }
    }
}

#[derive(Builder)]
#[builder(pattern = "owned")]
pub(crate) struct Caller {
    #[builder(private)]
    bcf_reader: bcf::IndexedReader,
    #[builder(private)]
    bcf_writer: bcf::Writer,
    #[builder(private)]
    contig_lens: HashMap<u32, u32>,
    #[builder(private)]
    alpha: f64
}

impl CallerBuilder {
    pub(crate) fn bcfs<P: AsRef<Path>>(mut self, in_path: Option<P>, out_path: Option<P>, alpha: f64) -> Result<Self> {
        self.alpha = alpha;

        self = self.bcf_reader(if let Some(path) = in_path {
            bcf::Reader::from_path(path)?
        } else {
            bcf::Reader::from_stdin()?
        });

        let bcf_reader = self.bcf_reader.as_ref().unwrap();

        let mut header = bcf::Header::new();
        for sample in bcf_reader.header().samples() {
            header.push_sample(sample);
        }

        header.push_record(
            "##INFO=<ID=LOHEND,Number=1,Type=Integer,Description=\"Last variant position supporting loss-of-heterozygosity region.\">"
                .as_bytes(),
        );

        let mut contig_lens = HashMap::new();
        // register sequences
        for rec in bcf_reader.header().header_records() {
            if let bcf::header::HeaderRecord::Contig { values, .. } = rec {
                let name = values.get("ID").unwrap();
                let len = values.get("length").unwrap();
                contig_lens.insert(HeaderView::name2id(name), len.parse()?);
                header.push_record(format!("##contig=<ID={},length={}>", name, len).as_bytes());
            }
        }

        self = self.contig_lens(contig_lens);

        Ok(self.bcf_writer(if let Some(path) = out_path {
            bcf::Writer::from_path(path, &header, false, bcf::Format::BCF)?
        } else {
            bcf::Writer::from_stdout(&header, false, bcf::Format::BCF)?
        }))
    }
}


impl Caller {
    pub(crate) fn call(&mut self) -> Result<()> {
        for (contig_id, contig_length) in self.contig_lens {

            // Problem Data
            let intervals = ContigLOHProbs::new(&mut self, &contig_id).create_all_intervals();

            // Define problem and objective sense
            let mut problem = LpProblem::new("LOH segmentation", LpObjective::Maximize);

            // Define Variables
            let interval_loh_indicator: BtreeMap<RangeInclusive<u64>, LpInteger> =
                intervals.iter()
                    .map(| (&range, &val)| {
                        let key = range;
                        let loh_indicator = LpInteger::new(&format!("{}", key))
                            .lower_bound(Some(0))
                            .upper_bound(Some(1));
                        (key, loh_indicator)
                    } )
                    .collect();

            // Define problem variables
            let ref n_intervals = LpInteger::new("number of intervals");

            // Define Objective Function
            let obj_vec: Vec<LpExpression> = {
                interval_loh_indicator.iter()
                    .map( |(&interval, loh_indicator)| {
                        loh_indicator * (interval.end() - interval.start() + 1)
                } )
            }.collect();
            problem += obj_vec.sum();

            // Constraint: no overlapping intervals
            for ( current_interval, _) in intervals {
                problem += sum(
                    &Vec::from(interval_loh_indicator.keys()),
                    | &interval | {
                        if current_interval.contains(interval.start()) | current_interval.contains(interval.end()) {
                            interval_loh_indicator.get(interval)?
                        } else {
                            0
                        }
                    }
                ).le(
                    if interval_loh_indicator.get(current_interval)? == 1 {
                        0
                    } else {
                        contig_length - ( current_interval.end() - current_interval.start() + 1)
                    }
                );
            }

            // Constraint: control false discovery rate at alpha
            problem += sum(
                    &Vec::from(interval_loh_indicator.keys()),
                    | &interval | {
                        interval_loh_indicator.get(interval)? * Prob::from( intervals.get(interval)?.ln_one_minus_exp() )
                    }
                ).le(sum(
                    &Vec::from(interval_loh_indicator.keys()),
                     | &interval | {
                         interval_loh_indicator.get(interval)?
                     }
                ) * self.alpha
            );

            // Specify solver
            let solver = CbcSolver::new();

            // Run optimisation and process output hashmap
            // Write result to bcf
            match solver.run(&problem) {
                Ok((status, var_values)) => {
                    println!("Status {:?}", status);
                    for (name, value) in var_values.iter() {
                        println!("value of {} = {}", name, value);
                    }
                },
                Err(msg) => println!("{}", msg),
            }
        }
        Ok(())
    }
}

pub(crate) fn info_phred_to_log_prob(
    record: &mut bcf::Record,
    info_field_name: &[u8]
) -> LogProb {
    let prob = record.info(info_field_name).float()?;
    if let Some(_prob) = prob {
        if !_prob[0].is_missing() && !_prob[0].is_nan() {
            let log_prob = LogProb::from(PHREDProb(_prob[0] as f64));
            assert!(
                * log_prob.is_valid(),
                "invalid PHRED probability '{}': {}, at pos: {}",
                info_field_name,
                _prob[0],
                record.pos()
            );
        }
        log_prob
    }
}

pub(crate) struct Interval {
    range: RangeInclusive<u64>,
    is_loh: Bool,
    posterior_prob_loh: LogProb,
}


#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq, PartialOrd, Ord)]
pub(crate) struct ContigLOHProbs {
    contig_id: u32,
    length: u64,
    cum_log_prob_loh: Vec<LogProb>,
}

impl ContigLOHProbs {
    pub(crate) fn new(
        caller: &mut Caller,
        contig_id: &u32
    ) -> ContigLOHProbs {
        let mut record = caller.bcf_reader.empty_record();
        let contig_length = caller.contig_lens.get(contig_id)?;
        let mut cum_log_prob_loh = Vec::with_capacity(contig_length as usize + 1 );
        caller.bcf_reader.fetch(*contig_id, 0, contig_length - 1);
        cum_log_prob_loh.push(LogProb::ln_one() );
        if caller.bcf_reader.read(&mut record) {
        }
        while caller.bcf_reader.read(&mut record) {
            cum_log_prob_loh.push(cum_log_prob_loh.last() + log_prob_loh_given_germ_het );
        }
        ContigLOHProbs {
            contig_id: *contig_id,
            length: contig_length as u64,
            cum_log_prob_loh: cum_log_prob_loh,
        }
    }

    pub(crate) fn create_all_intervals(
        &self
    ) -> BTreeMap<RangeInclusive<u64>, LogProb> {
        let mut intervals= BTreeMap::new();
        for start in 1..=self.length {
            for end in start..=self.length {
                intervals.insert(
                    start..=end,
                    self.cum_log_prob_loh[end] - self.cum_log_prob_loh[start - 1]
                )
            }
        }
        intervals
    }

    fn log_probs_loh_given_germ_het(
        record: &mut bcf::Record,
    ) -> LogProb {
        let log_prob_loh = info_phred_to_log_prob(record, b"PROB_LOH");
        let log_prob_no_loh = info_phred_to_log_prob(record, b"PROB_NO_LOH");
        let log_prob_germline_het = log_prob_loh.ln_add_exp(log_prob_no_loh);
        let log_prob_loh_given_germ_het = log_prob_not_germline_het.ln_add_exp(log_prob_loh + log_prob_germline_het);
        log_prob_loh_given_germ_het
    }
}


impl<'a> LOHRegion<'a> {
    pub(crate) fn write(
        &self,
        record: &mut bcf::Record,
    ) -> Result<()> {
        record.set_rid(Some(self.contig_id));
        record.set_pos(self.pos as i64);
        record.set_alleles(&[b"N", b"<LOH>"])?;
        record.push_info_integer(b"LOH_END", &[self.end as i32])?;
        record.push_info_integer(b"LOH_LEN", &[self.len() as i32])?;
        record.push_info_integer(b"N_LOCI", &[self.calls as i32])?;

        Ok(())
    }

    pub(crate) fn len(&self) -> u32 {
        (self.end - self.start + 1) as u32
    }
}

// //#[cfg(test)]
// mod tests {
//     use super::*;
//
//     #[test]
//     fn test_allele_freq_pdf() {
//         assert_eq!(
//             allele_freq_pdf(AlleleFreq(0.64), AlleleFreq(1.0), 10),
//             LogProb::ln_zero()
//         );
//         assert_eq!(
//             allele_freq_pdf(AlleleFreq(0.1), AlleleFreq(0.0), 10),
//             LogProb::ln_zero()
//         );
//     }
//
// }

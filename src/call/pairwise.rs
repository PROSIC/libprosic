use std::path::Path;
use std::error::Error;
use std::f32;
use std::str;

use itertools::Itertools;
use ndarray::prelude::*;
use csv;
use rust_htslib::bcf;
use rust_htslib::bcf::record::Numeric;
use bio::stats::{PHREDProb, LogProb};
use bio::io::fasta;

use model::AlleleFreqs;
use model::priors;
use model::PairCaller;
use model;
use ComplementEvent;
use Event;
use utils;


fn phred_scale<'a, I: IntoIterator<Item=&'a Option<LogProb>>>(probs: I) -> Vec<f32> {
    probs.into_iter().map(|&p| {
        match p {
            Some(p) => PHREDProb::from(p).abs() as f32,
            None    => f32::missing()
        }
    }).collect_vec()
}


pub struct PairEvent<A: AlleleFreqs, B: AlleleFreqs> {
    /// event name
    pub name: String,
    /// allele frequencies for case sample
    pub af_case: A,
    /// allele frequencies for control sample
    pub af_control: B
}


impl<A: AlleleFreqs, B: AlleleFreqs> Event for PairEvent<A, B> {
    fn name(&self) -> &str {
        &self.name
    }
}


fn pileups<'a, A, B, P>(
    inbcf: &bcf::Reader,
    record: &mut bcf::Record,
    joint_model: &'a mut PairCaller<A, B, P>,
    reference_buffer: &mut utils::ReferenceBuffer,
    omit_snvs: bool,
    omit_indels: bool,
    max_indel_len: Option<u32>,
    exclusive_end: bool
) -> Result<Vec<Option<model::PairPileup<'a, A, B, P>>>, Box<Error>> where
    A: AlleleFreqs,
    B: AlleleFreqs,
    P: priors::PairModel<A, B>
{
    let chrom = chrom(&inbcf, &record);
    let variants = try!(utils::collect_variants(record, omit_snvs, omit_indels, max_indel_len.map(|l| 0..l), exclusive_end));

    let chrom_seq = try!(reference_buffer.seq(&chrom));

    let mut pileups = Vec::with_capacity(variants.len());
    for variant in variants {
        pileups.push(if let Some(variant) = variant {
            Some(try!(joint_model.pileup(chrom, record.pos(), variant, chrom_seq)))
        } else {
            None
        });
    }

    Ok(pileups)
}


/// Call variants with the given model.
///
/// # Arguments
///
/// * `inbcf` - path to BCF/VCF with preprocessed variant calls (`"-"` for STDIN).
/// * `outbcf` - path to BCF/VCF with results (`"-"` for STDOUT).
/// * `events` - events to call
/// * `complement_event` - optional complementary event to call (e.g. absent)
/// * `joint_model` - calling model to use
/// * `omit_snvs` - omit single nucleotide variants
/// * `omit_indels` - omit indels
/// * `outobs` - optional path where to store observations as JSON
///
/// # Returns
///
/// `Result` object with eventual error message.
pub fn call<A, B, P, M, R, W, X, F>(
    inbcf: &R,
    outbcf: &W,
    fasta: &F,
    events: &[PairEvent<A, B>],
    complement_event: Option<&ComplementEvent>,
    pair_model: &mut PairCaller<A, B, P>,
    omit_snvs: bool,
    omit_indels: bool,
    max_indel_len: Option<u32>,
    outobs: Option<&X>,
    exclusive_end: bool
) -> Result<(), Box<Error>> where
    A: AlleleFreqs,
    B: AlleleFreqs,
    P: priors::PairModel<A, B>,
    R: AsRef<Path>,
    W: AsRef<Path>,
    X: AsRef<Path>,
    F: AsRef<Path>
{
    let fasta = try!(fasta::IndexedReader::from_file(fasta));
    let mut reference_buffer = utils::ReferenceBuffer::new(fasta);

    let inbcf = try!(bcf::Reader::new(inbcf));
    let mut header = bcf::Header::with_template(&inbcf.header);
    for event in events {
        header.push_record(
            event.header_entry("PROB", "PHRED-scaled probability for").as_bytes()
        );
    }
    if let Some(complement_event) = complement_event {
        header.push_record(complement_event.header_entry("PROB", "PHRED-scaled probability for").as_bytes());
    }
    // add tags for expected allele frequency
    header.push_record(
        b"##INFO=<ID=CASE_AF,Number=A,Type=Float,\
        Description=\"Maximum a posteriori probability estimate of allele frequency in case sample.\">"
    );
    header.push_record(
        b"##INFO=<ID=CONTROL_AF,Number=A,Type=Float,\
        Description=\"Maximum a posteriori probability estimate of allele frequency in control sample.\">"
    );

    let mut outbcf = try!(bcf::Writer::new(outbcf, &header, false, false));
    let mut outobs = if let Some(f) = outobs {
        let mut writer = try!(csv::Writer::from_file(f)).delimiter(b'\t');
        // write header for observations
        try!(writer.write(["chrom", "pos", "allele", "sample", "prob_mapping", "prob_alt", "prob_ref", "prob_mismapped", "evidence"].iter()));
        Some(writer)
    } else { None };
    let mut record = bcf::Record::new();
    let mut i = 0;
    loop {
        if let Err(e) = inbcf.read(&mut record) {
            if e.is_eof() {
                return Ok(())
            } else {
                return Err(Box::new(e));
            }
        }
        i += 1;
        // translate to header of the writer
        outbcf.translate(&mut record);
        let pileups = try!(pileups(&inbcf, &mut record, pair_model, &mut reference_buffer, omit_snvs, omit_indels, max_indel_len, exclusive_end));

        if !pileups.is_empty() && {
            let mut non_empty = true;
            for pileup in pileups.iter() {
                if let &Some(ref pileup) = pileup {
                    non_empty = non_empty && pileup.case_observations().len() > 0 && pileup.control_observations().len() > 0;
                }
            }
            non_empty
        }
/*        pileups.iter().fold(true, |non_empty, let &Some(ref x) = x| non_empty && (x > 0) ) &&
            pileups.iter().fold(true, |non_empty, &x| non_empty && (x.unwrap().control_observations().len() > 0) )*/{
            if let Some(ref mut outobs) = outobs {
                let chrom = str::from_utf8(chrom(&inbcf, &record)).unwrap();
                for (i, pileup) in pileups.iter().enumerate() {
                    if let &Some(ref pileup) = pileup {
                        for obs in pileup.case_observations() {
                            try!(outobs.encode((chrom, record.pos(), i, "case", obs)));
                        }
                        for obs in pileup.control_observations() {
                            try!(outobs.encode((chrom, record.pos(), i, "control", obs)));
                        }
                    }
                }
                try!(outobs.flush());
            }

            let mut posterior_probs = Array::default((events.len(), pileups.len()));
            for (i, event) in events.iter().enumerate() {
                for (j, pileup) in pileups.iter().enumerate() {
                    let p = if let &Some(ref pileup) = pileup {
                        Some(pileup.posterior_prob(&event.af_case, &event.af_control))
                    } else {
                        // indicate missing value
                        None
                    };

                    posterior_probs[(i, j)] = p;
                }
                try!(record.push_info_float(
                    event.tag_name("PROB").as_bytes(),
                    &phred_scale(posterior_probs.row(i).iter())
                ));
            }
            if let Some(complement_event) = complement_event {
                let mut complement_probs = Vec::with_capacity(pileups.len());
                for (j, pileup) in pileups.iter().enumerate() {
                    let p = if pileup.is_some() {
                        let event_probs = posterior_probs.column(j).iter().cloned().collect_vec();
                        let total = LogProb::ln_sum_exp(&event_probs.iter().map(|v| v.unwrap()).collect_vec());
                        debug!("Total probability over all defined Events: {}.", total.exp());
                        // total can slightly exceed 1 due to the numerical integration
                        Some(
                            if *total >= 0.0 {
                                LogProb::ln_zero()
                            } else {
                                total.ln_one_minus_exp()
                            }
                        )
                    } else {
                        // indicate missing value
                        None
                    };
                    complement_probs.push(p);
                }
                try!(record.push_info_float(
                    complement_event.tag_name("PROB").as_bytes(),
                    &phred_scale(complement_probs.iter())
                ));
            }
            let mut case_afs = Vec::with_capacity(pileups.len());
            let mut control_afs = Vec::with_capacity(pileups.len());
            for pileup in &pileups {
                if let &Some(ref pileup) = pileup {
                    let (case_af, control_af) = pileup.map_allele_freqs();
                    case_afs.push(*case_af as f32);
                    control_afs.push(*control_af as f32);
                } else {
                    case_afs.push(f32::missing());
                    control_afs.push(f32::missing());
                }
            }
            try!(record.push_info_float(b"CASE_AF", &case_afs));
            try!(record.push_info_float(b"CONTROL_AF", &control_afs));
        }
        try!(outbcf.write(&record));
        if i % 1000 == 0 {
            info!("{} records processed.", i);
        }
    }
}

fn chrom<'a>(inbcf: &'a bcf::Reader, record: &bcf::Record) -> &'a [u8] {
    inbcf.header.rid2name(record.rid().unwrap())
}

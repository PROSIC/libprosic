// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::HashMap;
use std::convert::{From, TryFrom};
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::error::Error;

use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use bio::stats::{LogProb, Prob};
use itertools::Itertools;
use rayon;
use rust_htslib::bam;
use serde_yaml;
use structopt;
use structopt::StructOpt;

use crate::calling;
use crate::conversion;
use crate::errors;
use crate::estimation;
use crate::estimation::alignment_properties::AlignmentProperties;
use crate::filtration;
use crate::grammar;
use crate::model::modes::generic::{FlatPrior, GenericModelBuilder};
use crate::model::sample::{estimate_alignment_properties, SampleBuilder};
use crate::model::{Contamination, VariantType};
use crate::testcase::TestcaseBuilder;
use crate::SimpleEvent;

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
#[structopt(
    name = "varlociraptor",
    about = "A caller for SNVs and indels in tumor-normal pairs.",
    setting = structopt::clap::AppSettings::ColoredHelp,
)]
pub enum Varlociraptor {
    #[structopt(
        name = "call",
        about = "Call variants.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Call {
        #[structopt(subcommand)]
        kind: CallKind,
    },
    #[structopt(
        name = "filter-calls",
        about = "Filter calls by either controlling the false discovery rate (FDR) at given level, or by posterior odds against the given events.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    FilterCalls {
        #[structopt(subcommand)]
        method: FilterMethod,
    },
    #[structopt(
        name = "decode-phred",
        about = "Decode PHRED-scaled values to human readable probabilities.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    DecodePHRED,
    #[structopt(
        name = "estimate",
        about = "Perform estimations.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    #[structopt(
        name = "estimate",
        about = "Perform estimations.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Estimate{
        #[structopt(subcommand)]
        kind: EstimateKind
    },
}

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
pub enum EstimateKind {
    #[structopt(
        name = "tmb",
        about = "Estimate tumor mutational burden. Takes Varlociraptor calls (must be annotated with e.g. snpEFF) from STDIN, prints TMB estimate in Vega-lite JSON format to STDOUT.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    TMB {
        #[structopt(
            long = "somatic-tumor-events",
            help = "Events to consider (e.g. SOMATIC_TUMOR).",
            required = true
        )]
        somatic_tumor_events: Vec<String>,
        #[structopt(
            long = "tumor-sample",
            help = "Name of the tumor sample in the given VCF/BCF."
        )]
        tumor_sample: String,
        #[structopt(
            long = "coding-genome-size",
            help = "Size of the covered coding genome."
        )]
        coding_genome_size: f64,
    }
}

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
pub enum CallKind {
    #[structopt(
        name = "variants",
        about = "Call variants.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Variants {
        #[structopt(subcommand)]
        mode: VariantCallMode,
        #[structopt(
            parse(from_os_str),
            help = "FASTA file with reference genome. Has to be indexed with samtools faidx."
        )]
        reference: PathBuf,
        #[structopt(
            parse(from_os_str),
            long,
            help = "VCF/BCF file to process (if omitted, read from STDIN)."
        )]
        candidates: Option<PathBuf>,
        #[structopt(
            parse(from_os_str),
            long,
            help = "BCF file that shall contain the results (if omitted, write to STDOUT)."
        )]
        output: Option<PathBuf>,
        #[structopt(
            long = "spurious-ins-rate",
            default_value = "2.8e-6",
            help = "Rate of spuriously inserted bases by the sequencer (Illumina: 2.8e-6, see Schirmer et al. BMC Bioinformatics 2016)."
        )]
        spurious_ins_rate: f64,
        #[structopt(
            long = "spurious-del-rate",
            default_value = "5.1e-6",
            help = "Rate of spuriosly deleted bases by the sequencer (Illumina: 5.1e-6, see Schirmer et al. BMC Bioinformatics 2016)."
        )]
        spurious_del_rate: f64,
        #[structopt(
            long = "spurious-insext-rate",
            default_value = "0.0",
            help = "Extension rate of spurious insertions by the sequencer (Illumina: 0.0, see Schirmer et al. BMC Bioinformatics 2016)"
        )]
        spurious_insext_rate: f64,
        #[structopt(
            long = "spurious-delext-rate",
            default_value = "0.0",
            help = "Extension rate of spurious deletions by the sequencer (Illumina: 0.0, see Schirmer et al. BMC Bioinformatics 2016)"
        )]
        spurious_delext_rate: f64,
        #[structopt(long = "omit-snvs", help = "Don't call SNVs.")]
        omit_snvs: bool,
        #[structopt(long = "omit-indels", help = "Don't call Indels.")]
        omit_indels: bool,
        #[structopt(
            long = "max-indel-len",
            default_value = "1000",
            help = "Omit longer indels when calling."
        )]
        max_indel_len: u32,
        #[structopt(
            long = "indel-window",
            default_value = "64",
            help = "Number of bases to consider left and right of indel breakpoint when \
                    calculating read support. This number should not be too large in order to \
                    avoid biases caused by other close variants. Currently implemented maximum \
                    value is 64."
        )]
        indel_window: u32,
        #[structopt(
            long = "max-depth",
            default_value = "200",
            help = "Maximum number of observations to use for calling. If locus is exceeding this \
                    number, downsampling is performed."
        )]
        max_depth: usize,
        #[structopt(
            long = "testcase-locus",
            help = "Create a test case for the given locus. Locus must be given in the form \
                    CHROM:POS[:IDX]. IDX is thereby an optional value to select a particular \
                    variant at the locus, counting from 1. If IDX is not specified, the first \
                    variant will be chosen. Alternatively, for single variant VCFs, you can \
                    specify 'all'."
        )]
        testcase_locus: Option<String>,
        #[structopt(
            long = "testcase-prefix",
            help = "Create test case files in the given directory."
        )]
        testcase_prefix: Option<String>,
    },
    #[structopt(
        name = "cnvs",
        about = "Call CNVs in tumor-normal sample pairs. This is experimental.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    CNVs {
        #[structopt(
            parse(from_os_str),
            long,
            help = "VCF/BCF file (generated by varlociraptor call-tumor-normal) to process \
                    (if omitted, read from STDIN)."
        )]
        calls: Option<PathBuf>,
        #[structopt(
            parse(from_os_str),
            long,
            help = "BCF file that shall contain the results (if omitted, write to STDOUT)."
        )]
        output: Option<PathBuf>,
        #[structopt(long, short = "p", help = "Tumor purity.")]
        purity: f64,
        #[structopt(
            long = "min-bayes-factor",
            default_value = "1.01",
            help = "Minimum bayes factor (> 1.0) between likelihoods of CNV and no CNV to consider. \
                    The higher this value, the fewer candidate CNVs will be investigated. \
                    Note that this can be usually left unchanged, because every CNV is provided \
                    with a posterior probability that can be used for filtering, e.g., via \
                    'varlociraptor control-fdr'."
        )]
        min_bayes_factor: f64,
        #[structopt(
            long,
            default_value = "1000",
            help = "Maximum distance between supporting loci in a CNV."
        )]
        max_dist: u32,
        #[structopt(long, short = "t", help = "Number of threads to use.")]
        threads: usize,
    },
}

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
pub enum VariantCallMode {
    #[structopt(
        name = "tumor-normal",
        about = "Call somatic and germline variants from a tumor-normal sample pair and a VCF/BCF with candidate variants.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    TumorNormal {
        #[structopt(parse(from_os_str), help = "BAM file with reads from tumor sample.")]
        tumor: PathBuf,
        #[structopt(parse(from_os_str), help = "BAM file with reads from normal sample.")]
        normal: PathBuf,
        #[structopt(short, long, default_value = "1.0", help = "Purity of tumor sample.")]
        purity: f64,
        #[structopt(
            parse(from_os_str),
            long = "tumor-alignment-properties",
            help = "Alignment properties JSON file for tumor sample. If not provided, properties \
                    will be estimated from the given BAM file."
        )]
        tumor_alignment_properties: Option<PathBuf>,
        #[structopt(
            parse(from_os_str),
            long = "normal-alignment-properties",
            help = "Alignment properties JSON file for normal sample. If not provided, properties \
                    will be estimated from the given BAM file."
        )]
        normal_alignment_properties: Option<PathBuf>,
    },
    #[structopt(
        name = "generic",
        about = "Call variants for a given scenario specified with the varlociraptor calling \
                 grammar and a VCF/BCF with candidate variants.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Generic {
        #[structopt(
            parse(from_os_str),
            long,
            help = "Scenario defined in the varlociraptor calling grammar."
        )]
        scenario: PathBuf,
        #[structopt(long, help = "BAM files with aligned reads for each sample.")]
        bams: Vec<String>,
        #[structopt(
            long = "alignment-properties",
            help = "Alignment properties JSON file for normal sample. If not provided, properties \
                    will be estimated from the given BAM file."
        )]
        alignment_properties: Vec<String>,
    },
}

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
pub enum FilterMethod {
    #[structopt(
        name = "control-fdr",
        about = "Filter variant calls by controlling FDR. Filtered calls a printed to STDOUT.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    ControlFDR {
        #[structopt(parse(from_os_str), help = "BCF file with varlociraptor calls.")]
        calls: PathBuf,
        #[structopt(
            long = "var",
            possible_values = { use strum::IntoEnumIterator; &VariantType::iter().map(|v| v.into()).collect_vec() },
            help = "Variant type to consider."
        )]
        vartype: VariantType,
        #[structopt(long, help = "FDR to control for.")]
        fdr: f64,
        #[structopt(long, help = "Events to consider.")]
        events: Vec<String>,
        #[structopt(long, help = "Minimum indel length to consider.")]
        minlen: Option<u32>,
        #[structopt(long, help = "Maximum indel length to consider (exclusive).")]
        maxlen: Option<u32>,
    },
    #[structopt(
        name = "posterior-odds",
        about = "Filter variant calls by posterior odds of given events against the rest of events.",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    PosteriorOdds {
        #[structopt(
            possible_values = { use strum::IntoEnumIterator; &KassRaftery::iter().map(|v| v.into()).collect_vec() },
            help = "Kass-Raftery score to filter against."
        )]
        odds: KassRaftery,
        #[structopt(long, help = "Events to consider.")]
        events: Vec<String>,
    },
}

fn parse_key_values(values: &[String]) -> Option<HashMap<String, PathBuf>> {
    let mut map = HashMap::new();
    for value in values {
        let kv = value.split_terminator('=').collect_vec();
        if kv.len() == 2 {
            map.insert(kv[0].to_owned(), PathBuf::from(kv[1]));
        } else {
            return None;
        }
    }
    Some(map)
}

impl Default for Varlociraptor {
    fn default() -> Self {
        Varlociraptor::from_iter(vec!["--help"])
    }
}

pub fn run(opt: Varlociraptor) -> Result<(), Box<dyn Error>> {
    let opt_clone = opt.clone();
    match opt {
        Varlociraptor::Call { kind } => {
            match kind {
                CallKind::Variants {
                    mode,
                    spurious_ins_rate,
                    spurious_del_rate,
                    spurious_insext_rate,
                    spurious_delext_rate,
                    indel_window,
                    omit_snvs,
                    omit_indels,
                    max_indel_len,
                    max_depth,
                    reference,
                    candidates,
                    output,
                    testcase_locus,
                    testcase_prefix,
                } => {
                    let spurious_ins_rate = Prob::checked(spurious_ins_rate)?;
                    let spurious_del_rate = Prob::checked(spurious_del_rate)?;
                    let spurious_insext_rate = Prob::checked(spurious_insext_rate)?;
                    let spurious_delext_rate = Prob::checked(spurious_delext_rate)?;
                    if indel_window > (128 / 2) {
                        Err(structopt::clap::Error::with_description( "Command-line option --indel-window requires a value <= 64 with the current implementation.", structopt::clap::ErrorKind::ValueValidation))?;
                    };
                    dbg!(indel_window);

                    let sample_builder = || {
                        SampleBuilder::default()
                            .error_probs(
                                spurious_ins_rate,
                                spurious_del_rate,
                                spurious_insext_rate,
                                spurious_delext_rate,
                                indel_window as u32,
                            )
                            .max_depth(max_depth)
                    };

                    let testcase_builder = if let Some(testcase_locus) = testcase_locus {
                        if let Some(testcase_prefix) = testcase_prefix {
                            if let Some(candidates) = candidates.as_ref() {
                                // just write a testcase and quit
                                Some(
                                    TestcaseBuilder::default()
                                        .prefix(PathBuf::from(testcase_prefix))
                                        .options(opt_clone)
                                        .locus(&testcase_locus)?
                                        .reference(&reference)?
                                        .candidates(candidates)?,
                                )
                            } else {
                                Err(errors::Error::MissingCandidates)?;
                                None
                            }
                        } else {
                            Err(errors::Error::MissingPrefix)?;
                            None
                        }
                    } else {
                        None
                    };

                    match mode {
                        VariantCallMode::Generic {
                            ref scenario,
                            ref bams,
                            ref alignment_properties,
                        } => {
                            if let Some(bams) = parse_key_values(bams) {
                                if let Some(alignment_properties) =
                                    parse_key_values(alignment_properties)
                                {
                                    if let Some(mut testcase_builder) = testcase_builder {
                                        for (name, bam) in &bams {
                                            testcase_builder =
                                                testcase_builder.register_bam(name, bam);
                                        }

                                        let mut testcase = testcase_builder
                                            .scenario(Some(scenario.to_owned()))
                                            .build()?;
                                        testcase.write()?;
                                        return Ok(());
                                    }

                                    let mut scenario_content = String::new();
                                    File::open(scenario)?.read_to_string(&mut scenario_content)?;

                                    let scenario: grammar::Scenario =
                                        serde_yaml::from_str(&scenario_content)?;
                                    let mut contaminations = scenario.sample_info();
                                    let mut resolutions = scenario.sample_info();
                                    let mut samples = scenario.sample_info();

                                    // parse samples
                                    for (sample_name, sample) in scenario.samples().iter() {
                                        let contamination =
                                            if let Some(contamination) = sample.contamination() {
                                                let contaminant = scenario
                                                .idx(contamination.by())
                                                .ok_or(
                                                errors::Error::InvalidContaminationSampleName {
                                                    name: sample_name.to_owned(),
                                                },
                                            )?;
                                                Some(Contamination {
                                                    by: contaminant,
                                                    fraction: *contamination.fraction(),
                                                })
                                            } else {
                                                None
                                            };
                                        contaminations =
                                            contaminations.push(sample_name, contamination);
                                        resolutions =
                                            resolutions.push(sample_name, *sample.resolution());

                                        let bam = bams.get(sample_name).ok_or(
                                            errors::Error::InvalidBAMSampleName {
                                                name: sample_name.to_owned(),
                                            },
                                        )?;
                                        let alignment_properties =
                                            est_or_load_alignment_properites(
                                                &alignment_properties.get(sample_name).as_ref(),
                                                bam,
                                            )?;
                                        let bam_reader = bam::IndexedReader::from_path(bam)?;
                                        let sample = sample_builder()
                                            .name(sample_name.to_owned())
                                            .alignments(bam_reader, alignment_properties)
                                            .build()?;
                                        samples = samples.push(sample_name, sample);
                                    }

                                    // register groups
                                    // for (sample_name, sample) in scenario.samples().iter() {
                                    //     if let Some(group) = sample.group() {
                                    //         sample_idx.insert(
                                    //             group,
                                    //             *sample_idx.get(sample_name).unwrap(),
                                    //         );
                                    //     }
                                    // }

                                    let model = GenericModelBuilder::default()
                                        // TODO allow to define prior in the grammar
                                        .prior(FlatPrior::new())
                                        .contaminations(contaminations.build())
                                        .resolutions(resolutions.build())
                                        .build()?;

                                    // setup caller
                                    let mut caller_builder = calling::variants::CallerBuilder::default()
                                        .samples(samples.build())
                                        .reference(reference)?
                                        .inbcf(candidates.as_ref())?
                                        .model(model)
                                        .omit_snvs(omit_snvs)
                                        .omit_indels(omit_indels)
                                        .max_indel_len(max_indel_len);
                                    for (event_name, vaftree) in scenario.vaftrees()? {
                                        caller_builder = caller_builder.event(&event_name, vaftree);
                                    }
                                    caller_builder = caller_builder.outbcf(output.as_ref())?;

                                    let mut caller = caller_builder.build()?;

                                    caller.call()?;
                                } else {
                                    Err(errors::Error::InvalidAlignmentPropertiesSpec)?
                                }
                            } else {
                                Err(errors::Error::InvalidBAMSpec)?
                            }
                        }
                        VariantCallMode::TumorNormal {
                            ref tumor,
                            ref normal,
                            purity,
                            ref tumor_alignment_properties,
                            ref normal_alignment_properties,
                        } => {
                            if let Some(testcase_builder) = testcase_builder {
                                let mut testcase = testcase_builder
                                    .register_bam("tumor", tumor)
                                    .register_bam("normal", normal)
                                    .scenario(None)
                                    .build()?;
                                testcase.write()?;
                                return Ok(());
                            }

                            let scenario = grammar::Scenario::try_from(
                                format!(
                                    r#"
                            samples:
                              tumor:
                                resolution: 100
                                contamination:
                                  by: normal
                                  fraction: {impurity}
                                universe: "[0.0,1.0]"
                              normal:
                                resolution: 5
                                universe: "[0.0,0.5[ | 0.5 | 1.0"
                            events:
                              somatic_tumor:  "tumor:]0.0,1.0] & normal:0.0"
                              somatic_normal: "tumor:]0.0,1.0] & normal:]0.0,0.5["
                              germline_het:   "tumor:]0.0,1.0] & normal:0.5"
                              germline_hom:   "tumor:]0.0,1.0] & normal:1.0"
                            "#,
                                    impurity = 1.0 - purity
                                )
                                .as_str(),
                            )?;

                            let tumor_alignment_properties = est_or_load_alignment_properites(
                                tumor_alignment_properties,
                                tumor,
                            )?;
                            let normal_alignment_properties = est_or_load_alignment_properites(
                                normal_alignment_properties,
                                normal,
                            )?;
                            info!("Estimated alignment properties:");
                            info!("{:?}", tumor_alignment_properties);
                            info!("{:?}", normal_alignment_properties);

                            let tumor_bam = bam::IndexedReader::from_path(tumor)?;
                            let normal_bam = bam::IndexedReader::from_path(normal)?;

                            let tumor_sample = sample_builder()
                                .name("tumor".to_owned())
                                .alignments(tumor_bam, tumor_alignment_properties)
                                .build()?;
                            let normal_sample = sample_builder()
                                .name("normal".to_owned())
                                .alignments(normal_bam, normal_alignment_properties)
                                .build()?;

                            let samples = scenario
                                .sample_info()
                                .push("tumor", tumor_sample)
                                .push("normal", normal_sample)
                                .build();
                            let contaminations = scenario
                                .sample_info()
                                .push(
                                    "tumor",
                                    Some(Contamination {
                                        by: scenario.idx("normal").unwrap(),
                                        fraction: 1.0 - purity,
                                    }),
                                )
                                .push("normal", None)
                                .build();
                            let resolutions = scenario
                                .sample_info()
                                .push("tumor", 100)
                                .push("normal", 5)
                                .build();

                            let model = GenericModelBuilder::default()
                                .prior(FlatPrior::new())
                                .contaminations(contaminations)
                                .resolutions(resolutions)
                                .build()?;

                            let mut caller_builder = calling::variants::CallerBuilder::default()
                                .samples(samples)
                                .reference(reference)?
                                .inbcf(candidates.as_ref())?
                                .model(model)
                                .omit_snvs(omit_snvs)
                                .omit_indels(omit_indels)
                                .max_indel_len(max_indel_len);

                            for (event_name, vaftree) in scenario.vaftrees()? {
                                caller_builder = caller_builder.event(&event_name, vaftree);
                            }

                            let mut caller = caller_builder.outbcf(output.as_ref())?.build()?;

                            caller.call()?;
                        }
                    }
                }
                CallKind::CNVs {
                    calls,
                    output,
                    min_bayes_factor,
                    threads,
                    purity,
                    max_dist,
                } => {
                    rayon::ThreadPoolBuilder::new()
                        .num_threads(threads)
                        .build_global()?;

                    if min_bayes_factor <= 1.0 {
                        Err(errors::Error::InvalidMinBayesFactor)?
                    }

                    let mut caller = calling::cnvs::CallerBuilder::default()
                        .bcfs(calls.as_ref(), output.as_ref())?
                        .min_bayes_factor(min_bayes_factor)
                        .purity(purity)
                        .max_dist(max_dist)
                        .build()?;
                    caller.call()?;
                }
            }
        }
        Varlociraptor::FilterCalls { method } => match method {
            FilterMethod::ControlFDR {
                calls,
                events,
                fdr,
                vartype,
                minlen,
                maxlen,
            } => {
                let events = events
                    .into_iter()
                    .map(|event| SimpleEvent {
                        name: event.to_owned(),
                    })
                    .collect_vec();
                let vartype = match (vartype, minlen, maxlen) {
                    (VariantType::Insertion(None), Some(minlen), Some(maxlen)) => {
                        VariantType::Insertion(Some(minlen..maxlen))
                    }
                    (VariantType::Deletion(None), Some(minlen), Some(maxlen)) => {
                        VariantType::Deletion(Some(minlen..maxlen))
                    }
                    (vartype @ _, _, _) => vartype.clone(),
                };

                filtration::fdr::control_fdr::<_, &PathBuf, &str>(
                    &calls,
                    None,
                    &events,
                    &vartype,
                    LogProb::from(Prob::checked(fdr)?),
                )?;
            }
            FilterMethod::PosteriorOdds { ref events, odds } => {
                let events = events
                    .into_iter()
                    .map(|event| SimpleEvent {
                        name: event.to_owned(),
                    })
                    .collect_vec();

                filtration::posterior_odds::filter_by_odds::<_, &PathBuf, &PathBuf>(
                    None, None, &events, odds,
                )?;
            }
        },
        Varlociraptor::DecodePHRED => {
            conversion::decode_phred::decode_phred()?;
        },
        Varlociraptor::Estimate { kind } => {
            match kind {
                EstimateKind::TMB {
                    somatic_tumor_events,
                    tumor_sample,
                    coding_genome_size,
                } => {
                    estimation::tumor_mutational_burden::estimate(&somatic_tumor_events, &tumor_sample, coding_genome_size as u64)?
                },
            }
        }
    }
    Ok(())
}

pub fn est_or_load_alignment_properites(
    alignment_properties_file: &Option<impl AsRef<Path>>,
    bam_file: impl AsRef<Path>,
) -> Result<AlignmentProperties, Box<dyn Error>> {
    if let Some(alignment_properties_file) = alignment_properties_file {
        Ok(serde_json::from_reader(File::open(
            alignment_properties_file,
        )?)?)
    } else {
        estimate_alignment_properties(bam_file)
    }
}

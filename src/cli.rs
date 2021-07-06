// Copyright 2016-2019 Johannes Köster, David Lähnemann.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::collections::HashMap;
use std::convert::{From, TryFrom};
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use anyhow::{Context, Result};
use bio::io::fasta;
use bio::stats::bayesian::bayes_factors::evidence::KassRaftery;
use bio::stats::{LogProb, Prob};
use itertools::Itertools;
use rust_htslib::bcf;
use structopt::StructOpt;
use strum::IntoEnumIterator;

use crate::calling;
use crate::conversion;
use crate::errors;
use crate::estimation;
use crate::estimation::alignment_properties::AlignmentProperties;
//use crate::estimation::sample_variants;
//use crate::estimation::tumor_mutational_burden;
use crate::filtration;
use crate::grammar;
use crate::reference;
use crate::testcase;
use crate::variants::evidence::realignment;
use crate::variants::evidence::realignment::pairhmm::GapParams;

use crate::variants::model::prior::CheckablePrior;
use crate::variants::model::prior::{Inheritance, Prior};
use crate::variants::model::{Contamination, VariantType};
use crate::variants::sample::{estimate_alignment_properties, ProtocolStrandedness};
use crate::variants::types::breakends::BreakendIndex;
use crate::SimpleEvent;

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
#[structopt(
    name = "varlociraptor",
    about = "A caller for SNVs and indels in tumor-normal pairs.",
    setting = structopt::clap::AppSettings::ColoredHelp,
)]
pub enum Varlociraptor {
    #[structopt(
        name = "preprocess",
        about = "Preprocess variants",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Preprocess {
        #[structopt(subcommand)]
        kind: PreprocessKind,
    },
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
    Estimate {
        #[structopt(subcommand)]
        kind: EstimateKind,
    },
    #[structopt(
        name = "plot",
        about = "Create plots",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Plot {
        #[structopt(subcommand)]
        kind: PlotKind,
    },
}

pub struct PreprocessInput {
    reference: PathBuf,
    bam: PathBuf,
}

impl Varlociraptor {
    fn preprocess_input(&self) -> PreprocessInput {
        if let Varlociraptor::Preprocess {
            kind:
                PreprocessKind::Variants {
                    ref reference,
                    ref bam,
                    ..
                },
        } = &self
        {
            PreprocessInput {
                reference: reference.to_owned(),
                bam: bam.to_owned(),
            }
        } else {
            panic!("bug: these are not preprocess options.");
        }
    }
}

fn default_reference_buffer_size() -> usize {
    10
}

fn default_pairhmm_mode() -> String {
    "exact".to_owned()
}

fn default_min_bam_refetch_distance() -> u64 {
    1
}

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
pub enum PreprocessKind {
    #[structopt(
        name = "variants",
        about = "Preprocess given variants by obtaining internal observations (allele likelihoods, strand information, ...)\
                 for each fragment. \
                 The obtained observations are printed to STDOUT in BCF format. Note that the resulting BCFs \
                 will be very large and are only intended for internal use (e.g. for piping into 'varlociraptor \
                 call variants generic').",
        usage = "varlociraptor preprocess variants reference.fasta --candidates candidates.bcf --bam sample.bam > sample.observations.bcf",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Variants {
        #[structopt(
            parse(from_os_str),
            help = "FASTA file with reference genome. Has to be indexed with samtools faidx."
        )]
        reference: PathBuf,
        #[structopt(
            parse(from_os_str),
            long,
            required = true,
            help = "VCF/BCF file to process (if omitted, read from STDIN)."
        )]
        candidates: PathBuf,
        #[structopt(
            long,
            required = true,
            help = "BAM file with aligned reads from a single sample."
        )]
        bam: PathBuf,
        #[structopt(
            long = "reference-buffer-size",
            short = "b",
            default_value = "10",
            help = "Number of reference sequences to keep in buffer. Use a smaller value \
                    to save memory at the expense of sometimes reduced parallelization."
        )]
        #[serde(default = "default_reference_buffer_size")]
        reference_buffer_size: usize,
        #[structopt(
            long = "min-bam-refetch-distance",
            default_value = "1",
            help = "Base pair distance to last fetched BAM interval such that a refetching is performed \
                  instead of reading through until the next interval is reached. Making this too small \
                  can cause unnecessary random access. Making this too large can lead to unneccessary \
                  iteration over irrelevant records. Benchmarking has shown that at least for short reads, \
                  a value of 1 (e.g. always refetch) does not incur additional costs and is a reasonable \
                  default."
        )]
        #[serde(default = "default_min_bam_refetch_distance")]
        min_bam_refetch_distance: u64,
        #[structopt(
            long = "alignment-properties",
            help = "Alignment properties JSON file for sample. If not provided, properties \
                    will be estimated from the given BAM file."
        )]
        alignment_properties: Option<PathBuf>,
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
        #[structopt(
            long = "strandedness",
            default_value = "opposite",
            possible_values = &ProtocolStrandedness::iter().map(|v| v.into()).collect_vec(),
            help = "Strandedness of sequencing protocol in case of paired-end (opposite strand as usual or same strand as with mate-pair sequencing.)"
        )]
        protocol_strandedness: ProtocolStrandedness,
        #[structopt(
            long = "indel-window",
            default_value = "64",
            help = "Number of bases to consider left and right of breakpoint when \
                    calculating read support. Currently implemented maximum \
                    value is 64."
        )]
        realignment_window: u64,
        #[structopt(
            long = "max-depth",
            default_value = "200",
            help = "Maximum number of observations to use for calling. If locus is exceeding this \
                    number, downsampling is performed."
        )]
        max_depth: usize,
        #[structopt(
            long = "omit-insert-size",
            help = "Do not consider insert size when calculating support for a variant. Use this flag when \
                    processing amplicon data, where indels do not impact the observed insert size"
        )]
        #[serde(default)]
        omit_insert_size: bool,
        #[structopt(
            long = "pairhmm-mode",
            possible_values = &["fast", "exact"],
            default_value = "exact",
            help = "PairHMM computation mode (either fast or exact). Fast mode means that only the best \
                    alignment path is considered for probability calculation. In rare cases, this can lead \
                    to wrong results for single reads. Hence, we advice to not use it when \
                    discrete allele frequences are of interest (0.5, 1.0). For continuous \
                    allele frequencies, fast mode should cause almost no deviations from the \
                    exact results. Also, if per sample allele frequencies are irrelevant (e.g. \
                    in large cohorts), fast mode can be safely used."
        )]
        #[serde(default = "default_pairhmm_mode")]
        pairhmm_mode: String,
    },
}

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
pub enum PlotKind {
    #[structopt(
        name = "variant-calling-prior",
        about = "Plot variant calling prior given a scenario. Plot is printed to STDOUT in Vega-lite format.",
        usage = "varlociraptor plot variant-calling-prior --scenario scenario.yaml > plot.vl.json",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    VariantCallingPrior {
        #[structopt(
            parse(from_os_str),
            long = "scenario",
            required = true,
            help = "Variant calling scenario that configures the prior."
        )]
        scenario: PathBuf,
        #[structopt(
            long = "contig",
            required = true,
            help = "Contig to consider for ploidy information."
        )]
        contig: String,
        #[structopt(long = "sample", required = true, help = "Sample to plot.")]
        sample: String,
    },
    #[structopt(
        name = "scatter",
        about = "Plot variant allelic fraction scatter plot overlayed with a contour plot between two sample groups",
        usage = "varlociraptor plot scatter --somatic-tumor-events SOMATIC_TUMOR \
        --sample-y sample1 --sample-x sample2 sample3 < calls.bcf | vg2svg > scatter.svg",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Scatter {
        #[structopt(
            long = "somatic-tumor-events",
            default_value = "SOMATIC_TUMOR",
            help = "Events to consider (e.g. SOMATIC_TUMOR)."
        )]
        somatic_tumor_events: Vec<String>,
        #[structopt(
            long = "sample-x",
            help = "Name of the first sample in the given VCF/BCF."
        )]
        sample_x: String,
        #[structopt(
            long = "sample-y",
            help = "Name(s) of the alternative sample(s) in the given VCF/BCF. Multiple samples can be given."
        )]
        sample_y: Vec<String>,
    },
}

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
pub enum EstimateKind {
    #[structopt(
        name = "mutational-burden",
        about = "Estimate mutational burden. Takes Varlociraptor calls (must be annotated \
                 with e.g. VEP but using ANN instead of CSQ) from STDIN, prints mutational burden estimate in Vega-lite JSON format to STDOUT. \
                 It can be converted to an image via vega-lite-cli (see conda package).",
        usage = "varlociraptor estimate mutational-burden --coding-genome-size 3e7 --events SOMATIC_TUMOR \
                 --sample tumor < calls.bcf | vg2svg > tmb.svg",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    MutationalBurden {
        #[structopt(long = "events", help = "Events to consider (e.g. SOMATIC_TUMOR).")]
        events: Vec<String>,
        #[structopt(
            long = "sample",
            help = "Name(s) of the sample(s) in the given VCF/BCF. Multiple samples can be given when using the multibar plot mode."
        )]
        sample: Vec<String>,
        #[structopt(
            long = "coding-genome-size",
            default_value = "3e7",
            help = "Size (in bases) of the covered coding genome."
        )]
        coding_genome_size: f64,
        #[structopt(
            long = "plot-mode",
            possible_values = &estimation::mutational_burden::PlotMode::iter().map(|v| v.into()).collect_vec(),
            help = "How to plot (as stratified curve, histogram or multi-sample barplot)."
        )]
        mode: estimation::mutational_burden::PlotMode,
        #[structopt(
            long = "vaf-cutoff",
            default_value = "0.2",
            help = "Minimal variant allelic fraction to consider for mutli-sample barplot"
        )]
        cutoff: f64,
    },
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
            long = "omit-strand-bias",
            help = "Do not consider strand bias when calculating the probability of an artifact.\
                    Use this flag when processing (panel) sequencing data, where the wet-lab \
                    methodology leads to strand bias in the coverage of genuine variants."
        )]
        #[serde(default)]
        omit_strand_bias: bool,
        #[structopt(
            long = "omit-read-orientation-bias",
            help = "Do not consider read orientation bias when calculating the probability of an \
                    artifact."
        )]
        #[serde(default)]
        omit_read_orientation_bias: bool,
        #[structopt(
            long = "omit-read-position-bias",
            help = "Do not consider read position bias when calculating the probability of an \
                    artifact. Use this flag when processing (panel) sequencing data, where the \
                    wet-lab methodology leads to stacks of reads starting at the same position."
        )]
        #[serde(default)]
        omit_read_position_bias: bool,
        #[structopt(
            long = "omit-softclip-bias",
            help = "Do not consider softclip bias when calculating the probability of an \
                    artifact. Use this flag when processing (panel) sequencing data, where the \
                    wet-lab methodology leads to stacks of reads starting at the same position."
        )]
        #[serde(default)]
        omit_softclip_bias: bool,
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
        #[structopt(
            long = "testcase-anonymous",
            help = "Anonymize any identifiers (via uuid4) and the sequences (by randomly permuting the alphabet) in the test case."
        )]
        testcase_anonymous: bool,
        #[structopt(
            long,
            short,
            help = "Output variant calls to given path (in BCF format). If omitted, prints calls to STDOUT."
        )]
        output: Option<PathBuf>,
    },

    #[structopt(
        name = "haplotype-abundances",
        about = "Call haplotype abundances (e.g. for HLA typing or viral strain quantification).",
        usage = "varlociraptor call haplotype-abundances --haplotype-counts counts.hdf5 \
                 --haplotype-variants calls.vcf --output results.tsv",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    HaplotypeCalls {
        #[structopt(
            parse(from_os_str),
            long = "haplotype-counts",
            required = true,
            help = "HDF5 haplotype counts calculated by Kallisto."
        )]
        haplotype_counts: PathBuf,
        #[structopt(
                parse(from_os_str),
                long = "haplotype-variants",
                required = true,
                help = "Haplotype variants compared to a common reference.", // TODO later, we will add a subcommand to generate this file with Varlociraptor as well
            )]
        haplotype_variants: PathBuf,
        #[structopt(
            default_value = "0.1",
            help = "Minimum value for normalized Kallisto counts."
        )]
        min_norm_counts: f64,
        #[structopt(
            long,
            help = "Folder to store quality control plots for the inference of a CDF from Kallisto bootstraps for each haplotype of interest."
        )]
        qc_plot: PathBuf,
    },
    // #[structopt(
    //     name = "cnvs",
    //     about = "Call CNVs in tumor-normal sample pairs. This is experimental (do not use it yet).",
    //     setting = structopt::clap::AppSettings::ColoredHelp,
    // )]
    // CNVs {
    //     #[structopt(
    //         parse(from_os_str),
    //         long,
    //         help = "VCF/BCF file (generated by varlociraptor call-tumor-normal) to process \
    //                 (if omitted, read from STDIN)."
    //     )]
    //     calls: Option<PathBuf>,
    //     #[structopt(
    //         parse(from_os_str),
    //         long,
    //         help = "BCF file that shall contain the results (if omitted, write to STDOUT)."
    //     )]
    //     output: Option<PathBuf>,
    //     #[structopt(long, short = "p", help = "Tumor purity.")]
    //     purity: f64,
    //     #[structopt(
    //         long = "min-bayes-factor",
    //         default_value = "1.01",
    //         help = "Minimum bayes factor (> 1.0) between likelihoods of CNV and no CNV to consider. \
    //                 The higher this value, the fewer candidate CNVs will be investigated. \
    //                 Note that this can be usually left unchanged, because every CNV is provided \
    //                 with a posterior probability that can be used for filtering, e.g., via \
    //                 'varlociraptor filter-calls control-fdr'."
    //     )]
    //     min_bayes_factor: f64,
    //     #[structopt(
    //         long,
    //         default_value = "1000",
    //         help = "Maximum distance between supporting loci in a CNV."
    //     )]
    //     max_dist: u32,
    //     #[structopt(long, short = "t", help = "Number of threads to use.")]
    //     threads: usize,
    // }
}

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
pub enum VariantCallMode {
    #[structopt(
        name = "tumor-normal",
        about = "Call somatic and germline variants from a tumor-normal sample pair and a VCF/BCF with candidate variants.",
        usage = "varlociraptor call variants tumor-normal --purity 0.75 --tumor tumor.bcf --normal normal.bcf > calls.bcf",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    TumorNormal {
        #[structopt(
            parse(from_os_str),
            long = "tumor",
            required = true,
            help = "BCF file with varlociraptor preprocess results for the tumor sample."
        )]
        tumor_observations: PathBuf,
        #[structopt(
            parse(from_os_str),
            long = "normal",
            required = true,
            help = "BCF file with varlociraptor preprocess results for the normal sample."
        )]
        normal_observations: PathBuf,
        #[structopt(short, long, default_value = "1.0", help = "Purity of tumor sample.")]
        purity: f64,
    },
    #[structopt(
        name = "generic",
        about = "Call variants for a given scenario specified with the varlociraptor calling \
                 grammar and a VCF/BCF with candidate variants.",
        usage = "varlociraptor call variants --output calls.bcf generic --scenario scenario.yaml \
                 --obs relapse=relapse.bcf tumor=tumor.bcf normal=normal.bcf",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    Generic {
        #[structopt(
            parse(from_os_str),
            long,
            required = true,
            help = "Scenario defined in the varlociraptor calling grammar."
        )]
        scenario: PathBuf,
        #[structopt(
            long = "obs",
            required = true,
            help = "BCF file with varlociraptor preprocess results for samples defined in the given scenario (given as samplename=path/to/calls.bcf). It is possible to omit a sample here (e.g. model tumor/normal in the scenario, but only call on the tumor sample when there is no normal sequenced). In that case, the resulting probabilities will be accordingly uncertain, because observations of the omitted sample are missing (which is equivalent to having no coverage in the sample)."
        )]
        sample_observations: Vec<String>,
    },
}

#[derive(Debug, StructOpt, Serialize, Deserialize, Clone)]
pub enum FilterMethod {
    #[structopt(
        name = "control-fdr",
        about = "Filter variant calls by controlling FDR. Filtered calls are printed to STDOUT.",
        usage = "varlociraptor filter-calls control-fdr calls.bcf --events SOMATIC_TUMOR --fdr 0.05 \
                 --var SNV > calls.filtered.bcf",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    ControlFDR {
        #[structopt(parse(from_os_str), help = "BCF file with varlociraptor calls.")]
        calls: PathBuf,
        #[structopt(
            long = "var",
            possible_values = &VariantType::iter().map(|v| v.into()).collect_vec(),
            help = "Variant type to consider. When controlling global FDR (not using --local) this should \
            be used to control FDR for each type separately. Otherwise, less certain variant types will be \
            underrepresented."
        )]
        vartype: Option<VariantType>,
        #[structopt(long, help = "FDR to control for.")]
        fdr: f64,
        #[structopt(
            long = "local",
            help = "Control local FDR instead of global FDR. This means that for each record, the posterior \
            of the selected events has to be at least 1-fdr."
        )]
        local: bool,
        #[structopt(long, help = "Events to consider.")]
        events: Vec<String>,
        #[structopt(long, help = "Minimum indel length to consider.")]
        minlen: Option<u64>,
        #[structopt(long, help = "Maximum indel length to consider (exclusive).")]
        maxlen: Option<u64>,
    },
    #[structopt(
        name = "posterior-odds",
        about = "Filter variant calls by posterior odds of given events against the rest of events. \
                 Calls are taken from STDIN, filtered calls are printed to STDOUT.",
        usage = "varlociraptor filter-calls posterior-odds --events SOMATIC_TUMOR --odds strong < calls.bcf",
        setting = structopt::clap::AppSettings::ColoredHelp,
    )]
    PosteriorOdds {
        #[structopt(
            long,
            possible_values = &KassRaftery::iter().map(|v| v.into()).collect_vec(),
            help = "Kass-Raftery score to filter against."
        )]
        odds: KassRaftery,
        #[structopt(long, help = "Events to consider.")]
        events: Vec<String>,
    },
}

pub(crate) type PathMap = HashMap<String, PathBuf>;

fn parse_key_values(values: &[String]) -> Option<PathMap> {
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

pub fn run(opt: Varlociraptor) -> Result<()> {
    let opt_clone = opt.clone();
    match opt {
        Varlociraptor::Preprocess { kind } => {
            match kind {
                PreprocessKind::Variants {
                    reference,
                    candidates,
                    bam,
                    alignment_properties,
                    output,
                    spurious_ins_rate,
                    spurious_del_rate,
                    spurious_insext_rate,
                    spurious_delext_rate,
                    protocol_strandedness,
                    realignment_window,
                    max_depth,
                    omit_insert_size,
                    reference_buffer_size,
                    min_bam_refetch_distance,
                    pairhmm_mode,
                } => {
                    // TODO: handle testcases

                    let spurious_ins_rate = Prob::checked(spurious_ins_rate)?;
                    let spurious_del_rate = Prob::checked(spurious_del_rate)?;
                    let spurious_insext_rate = Prob::checked(spurious_insext_rate)?;
                    let spurious_delext_rate = Prob::checked(spurious_delext_rate)?;
                    if realignment_window > (128 / 2) {
                        return Err(
                            structopt::clap::Error::with_description(
                                "Command-line option --indel-window requires a value <= 64 with the current implementation.", 
                                structopt::clap::ErrorKind::ValueValidation
                            ).into()
                        );
                    };

                    // If we omit the insert size information for calculating the evidence, we can savely allow hardclips here.
                    let allow_hardclips = omit_insert_size;
                    let alignment_properties = est_or_load_alignment_properties(
                        &alignment_properties,
                        &bam,
                        omit_insert_size,
                        allow_hardclips,
                    )?;

                    let gap_params = GapParams {
                        prob_insertion_artifact: LogProb::from(spurious_ins_rate),
                        prob_deletion_artifact: LogProb::from(spurious_del_rate),
                        prob_insertion_extend_artifact: LogProb::from(spurious_insext_rate),
                        prob_deletion_extend_artifact: LogProb::from(spurious_delext_rate),
                    };

                    let reference_buffer = Arc::new(reference::Buffer::new(
                        fasta::IndexedReader::from_file(&reference)
                            .context("Unable to read genome reference.")?,
                        reference_buffer_size,
                    ));

                    if pairhmm_mode == "fast" {
                        let mut processor =
                            calling::variants::preprocessing::ObservationProcessor::builder()
                                .alignment_properties(alignment_properties)
                                .protocol_strandedness(protocol_strandedness)
                                .max_depth(max_depth)
                                .inbam(bam)
                                .min_bam_refetch_distance(min_bam_refetch_distance)
                                .reference_buffer(Arc::clone(&reference_buffer))
                                .breakend_index(BreakendIndex::new(&candidates)?)
                                .inbcf(candidates)
                                .options(opt_clone)
                                .outbcf(output)
                                .realigner(realignment::PathHMMRealigner::new(
                                    gap_params,
                                    realignment_window,
                                    reference_buffer,
                                ))
                                .build();

                        processor.process()?;
                    } else {
                        let mut processor =
                            calling::variants::preprocessing::ObservationProcessor::builder()
                                .alignment_properties(alignment_properties)
                                .protocol_strandedness(protocol_strandedness)
                                .max_depth(max_depth)
                                .inbam(bam)
                                .min_bam_refetch_distance(min_bam_refetch_distance)
                                .reference_buffer(Arc::clone(&reference_buffer))
                                .breakend_index(BreakendIndex::new(&candidates)?)
                                .inbcf(candidates)
                                .options(opt_clone)
                                .outbcf(output)
                                .realigner(realignment::PairHMMRealigner::new(
                                    reference_buffer,
                                    gap_params,
                                    realignment_window,
                                ))
                                .build();

                        processor.process()?;
                    }
                }
            }
        }
        Varlociraptor::Call { kind } => {
            match kind {
                CallKind::Variants {
                    mode,
                    omit_strand_bias,
                    omit_read_orientation_bias,
                    omit_read_position_bias,
                    omit_softclip_bias,
                    testcase_locus,
                    testcase_prefix,
                    testcase_anonymous,
                    output,
                } => {
                    let testcase_builder = if let Some(testcase_locus) = testcase_locus {
                        if let Some(testcase_prefix) = testcase_prefix {
                            // TODO obtain sample information from input bcfs?
                            Some(
                                testcase::TestcaseBuilder::default()
                                    .prefix(PathBuf::from(testcase_prefix))
                                    .anonymize(testcase_anonymous)
                                    .locus(&testcase_locus)?,
                            )
                        } else {
                            return Err(errors::Error::MissingPrefix.into());
                        }
                    } else {
                        None
                    };

                    let call_generic = |scenario: grammar::Scenario,
                                        observations: PathMap|
                     -> Result<()> {
                        let sample_infos = SampleInfos::try_from(&scenario)?;

                        // record observation paths
                        let mut sample_observations = scenario.sample_info();
                        for (sample_name, _) in scenario.samples().iter() {
                            if let Some(obs) = observations.get(sample_name) {
                                sample_observations =
                                    sample_observations.push(sample_name, Some(obs.to_owned()));
                            } else {
                                sample_observations = sample_observations.push(sample_name, None);
                            }
                        }
                        let sample_observations = sample_observations.build();

                        for obs_sample_name in observations.keys() {
                            if !sample_infos.names.as_slice().contains(obs_sample_name) {
                                return Err(errors::Error::InvalidObservationSampleName {
                                    name: obs_sample_name.to_owned(),
                                }
                                .into());
                            }
                        }

                        let breakend_index =
                            BreakendIndex::new(sample_observations.first_not_none()?)?;

                        let prior = Prior::builder()
                            .ploidies(None)
                            .universe(None)
                            .uniform(sample_infos.uniform_prior)
                            .germline_mutation_rate(sample_infos.germline_mutation_rates)
                            .somatic_effective_mutation_rate(
                                sample_infos.somatic_effective_mutation_rates,
                            )
                            .inheritance(sample_infos.inheritance)
                            .genome_size(
                                scenario
                                    .species()
                                    .as_ref()
                                    .map_or(None, |species| *species.genome_size()),
                            )
                            .heterozygosity(scenario.species().as_ref().map_or(None, |species| {
                                species.heterozygosity().map(|het| LogProb::from(Prob(het)))
                            }))
                            .variant_type_fractions(scenario.variant_type_fractions())
                            .build();

                        // setup caller
                        let caller = calling::variants::CallerBuilder::default()
                            .samplenames(sample_infos.names)
                            .observations(sample_observations)
                            .omit_strand_bias(omit_strand_bias)
                            .omit_read_orientation_bias(omit_read_orientation_bias)
                            .omit_read_position_bias(omit_read_position_bias)
                            .omit_softclip_bias(omit_softclip_bias)
                            .scenario(scenario)
                            .prior(prior)
                            .contaminations(sample_infos.contaminations)
                            .resolutions(sample_infos.resolutions)
                            .breakend_index(breakend_index)
                            .outbcf(output)
                            .build()
                            .unwrap();

                        // call
                        caller.call()?;

                        Ok(())
                    };

                    match mode {
                        VariantCallMode::Generic {
                            scenario,
                            sample_observations,
                        } => {
                            if let Some(sample_observations) =
                                parse_key_values(&sample_observations)
                            {
                                if let Some(mut testcase_builder) = testcase_builder {
                                    for (i, (sample_name, obspath)) in
                                        sample_observations.iter().enumerate()
                                    {
                                        let options = calling::variants::preprocessing::read_preprocess_options(obspath)?;
                                        let preprocess_input = options.preprocess_input();
                                        testcase_builder = testcase_builder.register_sample(
                                            &sample_name,
                                            preprocess_input.bam,
                                            options,
                                        )?;
                                        if i == 0 {
                                            testcase_builder =
                                                testcase_builder.candidates(obspath.to_owned());
                                            testcase_builder = testcase_builder
                                                .reference(preprocess_input.reference)?;
                                        }
                                    }

                                    let mut testcase = testcase_builder
                                        .scenario(Some(scenario))
                                        .mode(testcase::Mode::Generic)
                                        .build()
                                        .unwrap();
                                    info!("Writing testcase.");
                                    testcase.write()?;
                                    return Ok(());
                                }

                                let scenario = grammar::Scenario::from_path(scenario)?;

                                call_generic(scenario, sample_observations)?;
                            } else {
                                return Err(errors::Error::InvalidObservationsSpec.into());
                            }
                        }
                        VariantCallMode::TumorNormal {
                            tumor_observations,
                            normal_observations,
                            purity,
                        } => {
                            if let Some(testcase_builder) = testcase_builder {
                                let tumor_options =
                                    calling::variants::preprocessing::read_preprocess_options(
                                        &tumor_observations,
                                    )?;
                                let normal_options =
                                    calling::variants::preprocessing::read_preprocess_options(
                                        &normal_observations,
                                    )?;
                                let mut testcase = testcase_builder
                                    .candidates(tumor_observations)
                                    .reference(tumor_options.preprocess_input().reference)?
                                    .register_sample(
                                        "tumor",
                                        tumor_options.preprocess_input().bam,
                                        tumor_options,
                                    )?
                                    .register_sample(
                                        "normal",
                                        normal_options.preprocess_input().bam,
                                        normal_options,
                                    )?
                                    .scenario(None)
                                    .mode(testcase::Mode::TumorNormal)
                                    .purity(Some(purity))
                                    .build()
                                    .unwrap();

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

                            let mut observations = PathMap::default();
                            observations.insert("tumor".to_owned(), tumor_observations);
                            observations.insert("normal".to_owned(), normal_observations);

                            call_generic(scenario, observations)?;
                        }
                    }
                } // CallKind::CNVs {
                //     calls,
                //     output,
                //     min_bayes_factor,
                //     threads,
                //     purity,
                //     max_dist,
                // } => {
                //     rayon::ThreadPoolBuilder::new()
                //         .num_threads(threads)
                //         .build_global()?;

                //     if min_bayes_factor <= 1.0 {
                //         Err(errors::Error::InvalidMinBayesFactor)?
                //     }

                //     let mut caller = calling::cnvs::CallerBuilder::default()
                //         .bcfs(calls.as_ref(), output.as_ref())?
                //         .min_bayes_factor(min_bayes_factor)
                //         .purity(purity)
                //         .max_dist(max_dist)
                //         .build()
                //         .unwrap();
                //     caller.call()?;
                // }
                CallKind::HaplotypeCalls {
                    haplotype_counts,
                    haplotype_variants,
                    qc_plot,
                    min_norm_counts,
                } => {
                    let caller = calling::haplotype_abundances::CallerBuilder::default()
                        .hdf5_reader(hdf5::File::open(&haplotype_counts)?)
                        .vcf_reader(bcf::Reader::from_path(&haplotype_variants)?)
                        .min_norm_counts(min_norm_counts)
                        .build()
                        .unwrap();
                    let seqnames = caller.filter_seqnames().unwrap();
                    for seqname in seqnames {
                        let ecdf = caller.cdf(seqname).unwrap();
                        ecdf.plot_qc(qc_plot.clone()).unwrap();
                    }
                }
            }
        }

        Varlociraptor::FilterCalls { method } => match method {
            FilterMethod::ControlFDR {
                calls,
                events,
                fdr,
                local,
                vartype,
                minlen,
                maxlen,
            } => {
                let events = events
                    .into_iter()
                    .map(|event| SimpleEvent { name: event })
                    .collect_vec();
                let vartype = match (vartype, minlen, maxlen) {
                    (Some(VariantType::Insertion(None)), Some(minlen), Some(maxlen)) => {
                        Some(VariantType::Insertion(Some(minlen..maxlen)))
                    }
                    (Some(VariantType::Deletion(None)), Some(minlen), Some(maxlen)) => {
                        Some(VariantType::Deletion(Some(minlen..maxlen)))
                    }
                    (vartype, _, _) => vartype,
                };

                filtration::fdr::control_fdr::<_, &PathBuf, &str>(
                    &calls,
                    None,
                    &events,
                    vartype.as_ref(),
                    LogProb::from(Prob::checked(fdr)?),
                    local,
                )?;
            }
            FilterMethod::PosteriorOdds { ref events, odds } => {
                let events = events
                    .iter()
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
        }
        Varlociraptor::Estimate { kind } => match kind {
            EstimateKind::MutationalBurden {
                events,
                sample,
                coding_genome_size,
                mode,
                cutoff,
            } => estimation::mutational_burden::collect_estimates(
                &events,
                &sample,
                coding_genome_size as u64,
                mode,
                cutoff as f64,
            )?,
        },
        Varlociraptor::Plot { kind } => match kind {
            PlotKind::VariantCallingPrior {
                scenario,
                contig,
                sample,
            } => {
                let scenario = grammar::Scenario::from_path(scenario)?;
                let sample_infos = SampleInfos::try_from(&scenario)?;

                let mut universes = scenario.sample_info();
                let mut ploidies = scenario.sample_info();
                for (sample_name, sample) in scenario.samples().iter() {
                    universes = universes.push(
                        sample_name,
                        sample.contig_universe(&contig, scenario.species())?,
                    );
                    ploidies = ploidies.push(
                        sample_name,
                        sample.contig_ploidy(&contig, scenario.species())?,
                    );
                }
                let universes = universes.build();
                let ploidies = ploidies.build();

                let prior = Prior::builder()
                    .variant_type_fractions(scenario.variant_type_fractions())
                    .ploidies(Some(ploidies))
                    .universe(Some(universes))
                    .uniform(sample_infos.uniform_prior)
                    .germline_mutation_rate(sample_infos.germline_mutation_rates)
                    .somatic_effective_mutation_rate(sample_infos.somatic_effective_mutation_rates)
                    .inheritance(sample_infos.inheritance)
                    .genome_size(
                        scenario
                            .species()
                            .as_ref()
                            .map_or(None, |species| *species.genome_size()),
                    )
                    .heterozygosity(scenario.species().as_ref().map_or(None, |species| {
                        species.heterozygosity().map(|het| LogProb::from(Prob(het)))
                    }))
                    .build();
                prior.check()?;

                prior.plot(&sample, &sample_infos.names)?;
            }
            PlotKind::Scatter {
                somatic_tumor_events,
                sample_x,
                sample_y,
            } => estimation::sample_variants::vaf_scatter(
                &somatic_tumor_events,
                &sample_x,
                &sample_y,
            )?,
        },
    }
    Ok(())
}

pub(crate) fn est_or_load_alignment_properties(
    alignment_properties_file: &Option<impl AsRef<Path>>,
    bam_file: impl AsRef<Path>,
    omit_insert_size: bool,
    allow_hardclips: bool,
) -> Result<AlignmentProperties> {
    if let Some(alignment_properties_file) = alignment_properties_file {
        Ok(serde_json::from_reader(File::open(
            alignment_properties_file,
        )?)?)
    } else {
        estimate_alignment_properties(bam_file, omit_insert_size, allow_hardclips)
    }
}

struct SampleInfos {
    uniform_prior: grammar::SampleInfo<bool>,
    contaminations: grammar::SampleInfo<Option<Contamination>>,
    resolutions: grammar::SampleInfo<usize>,
    germline_mutation_rates: grammar::SampleInfo<Option<f64>>,
    somatic_effective_mutation_rates: grammar::SampleInfo<Option<f64>>,
    inheritance: grammar::SampleInfo<Option<Inheritance>>,
    names: grammar::SampleInfo<String>,
}

impl<'a> TryFrom<&'a grammar::Scenario> for SampleInfos {
    type Error = anyhow::Error;

    fn try_from(scenario: &grammar::Scenario) -> Result<Self> {
        let mut contaminations = scenario.sample_info();
        let mut resolutions = scenario.sample_info();
        let mut sample_names = scenario.sample_info();
        let mut germline_mutation_rates = scenario.sample_info();
        let mut somatic_effective_mutation_rates = scenario.sample_info();
        let mut inheritance = scenario.sample_info();
        let mut uniform_prior = scenario.sample_info();

        for (sample_name, sample) in scenario.samples().iter() {
            let contamination = if let Some(contamination) = sample.contamination() {
                let contaminant = scenario.idx(contamination.by()).ok_or(
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
            uniform_prior = uniform_prior.push(sample_name, sample.has_uniform_prior());
            contaminations = contaminations.push(sample_name, contamination);
            resolutions = resolutions.push(sample_name, *sample.resolution());
            sample_names = sample_names.push(sample_name, sample_name.to_owned());
            germline_mutation_rates = germline_mutation_rates.push(
                sample_name,
                sample.germline_mutation_rate(scenario.species()),
            );
            somatic_effective_mutation_rates = somatic_effective_mutation_rates.push(
                sample_name,
                sample.somatic_effective_mutation_rate(scenario.species()),
            );
            inheritance = inheritance.push(
                sample_name,
                if let Some(inheritance) = sample.inheritance() {
                    let parent_idx = |parent| {
                        scenario
                            .idx(parent)
                            .ok_or(errors::Error::InvalidInheritanceSampleName {
                                name: sample_name.to_owned(),
                            })
                    };
                    Some(match inheritance {
                        grammar::Inheritance::Mendelian { from: parents } => {
                            Inheritance::Mendelian {
                                from: (parent_idx(&parents.0)?, parent_idx(&parents.1)?),
                            }
                        }
                        grammar::Inheritance::Clonal {
                            from: parent,
                            somatic,
                        } => Inheritance::Clonal {
                            from: parent_idx(parent)?,
                            somatic: *somatic,
                        },
                        grammar::Inheritance::Subclonal {
                            from: parent,
                            origin,
                        } => Inheritance::Subclonal {
                            from: parent_idx(parent)?,
                            origin: *origin,
                        },
                    })
                } else {
                    None
                },
            );
        }

        Ok(SampleInfos {
            uniform_prior: uniform_prior.build(),
            contaminations: contaminations.build(),
            resolutions: resolutions.build(),
            germline_mutation_rates: germline_mutation_rates.build(),
            somatic_effective_mutation_rates: somatic_effective_mutation_rates.build(),
            inheritance: inheritance.build(),
            names: sample_names.build(),
        })
    }
}

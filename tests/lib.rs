use std::error::Error;

use std::fs;
use std::fs::File;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::str;

use bio::io::fasta;
use bio::stats::{LogProb, Prob};
use eval::Expr;
use itertools::Itertools;
use rust_htslib::bcf::Read as BCFRead;
use rust_htslib::{bam, bcf};
use serde_json;
use tempfile::{self, NamedTempFile};
use yaml_rust::{Yaml, YamlLoader};

use varlociraptor::cli::{run, CallKind, VariantCallMode, Varlociraptor};

struct Testcase {
    inner: Vec<Yaml>,
    path: PathBuf,
}

impl Testcase {
    fn new(path: impl AsRef<Path>) -> Result<Self, Box<Error>> {
        let mut reader = File::open(path.as_ref().join("testcase.yaml"))?;
        let mut content2 = String::new();
        reader.read_to_string(&mut content2)?;
        Ok(Testcase {
            inner: YamlLoader::load_from_str(&content2)?,
            path: path.as_ref().to_owned(),
        })
    }

    fn yaml(&self) -> &Yaml {
        &self.inner[0]
    }

    fn run(&self) -> Result<(), Box<Error>> {
        let mut options = serde_json::from_str(self.yaml()["options"].as_str().unwrap())?;
        let temp_ref = Self::reference(
            self.yaml()["reference"]["name"].as_str().unwrap(),
            self.yaml()["reference"]["seq"].as_str().unwrap(),
        )?;

        match &mut options {
            Varlociraptor::Call { ref mut kind } => match kind {
                CallKind::Variants {
                    ref mut mode,
                    ref mut reference,
                    ref mut candidates,
                    ref mut output,
                    ref mut testcase_locus,
                    ref mut testcase_prefix,
                    ..
                } => {
                    *reference = temp_ref.path().to_owned();
                    *candidates = Some(self.path.join(self.yaml()["candidate"].as_str().unwrap()));
                    *output = Some(self.output());
                    *testcase_prefix = None;
                    *testcase_locus = None;

                    match mode {
                        VariantCallMode::Generic {
                            ref mut scenario,
                            ref mut bams,
                            ref mut alignment_properties,
                        } => {
                            *scenario = self.path.join(self.yaml()["scenario"].as_str().unwrap());
                            bams.clear();
                            alignment_properties.clear();
                            let mut temp_props = Vec::new();
                            for (sample_name, sample) in
                                self.yaml()["samples"].as_hash().unwrap().iter()
                            {
                                let sample_name = sample_name.as_str().unwrap();
                                let bam = self.path.join(sample["path"].as_str().unwrap());
                                bam::index::build(&bam, None, bam::index::Type::BAI, 1).unwrap();
                                bams.push(format!("{}={}", sample_name, bam.to_str().unwrap()));
                                let props = Self::alignment_properties(
                                    sample["properties"].as_str().unwrap(),
                                )?;
                                alignment_properties.push(format!(
                                    "{}={}",
                                    sample_name,
                                    props.path().to_str().unwrap()
                                ));
                                temp_props.push(props);
                            }
                            run(options)
                        }
                        VariantCallMode::TumorNormal {
                            ref mut tumor,
                            ref mut normal,
                            ref mut tumor_alignment_properties,
                            ref mut normal_alignment_properties,
                            ..
                        } => {
                            *tumor = self
                                .path
                                .join(self.yaml()["samples"]["tumor"]["path"].as_str().unwrap());
                            *normal = self
                                .path
                                .join(self.yaml()["samples"]["normal"]["path"].as_str().unwrap());

                            let temp_tumor_props = Self::alignment_properties(
                                self.yaml()["samples"]["tumor"]["properties"]
                                    .as_str()
                                    .unwrap(),
                            )?;
                            let temp_normal_props = Self::alignment_properties(
                                self.yaml()["samples"]["normal"]["properties"]
                                    .as_str()
                                    .unwrap(),
                            )?;
                            *tumor_alignment_properties = Some(temp_tumor_props.path().to_owned());
                            *normal_alignment_properties =
                                Some(temp_normal_props.path().to_owned());

                            bam::index::build(tumor, None, bam::index::Type::BAI, 1).unwrap();
                            bam::index::build(normal, None, bam::index::Type::BAI, 1).unwrap();

                            run(options)
                        }
                    }
                }
                _ => panic!("unsupported subcommand"),
            },
            _ => panic!("unsupported subcommand"),
        }
    }

    fn output(&self) -> PathBuf {
        self.path.join("calls.bcf")
    }

    fn check(&self) {
        let mut reader = bcf::Reader::from_path(self.output()).unwrap();
        let mut calls = reader.records().map(|r| r.unwrap()).collect_vec();
        assert_eq!(calls.len(), 1, "unexpected number of calls");
        let mut call = calls.pop().unwrap();

        let afs = call.format(b"AF").float().unwrap();
        if let Some(exprs) = self.yaml()["expected"]["allelefreqs"].as_vec() {
            for expr in exprs.iter() {
                let mut expr = Expr::new(expr.as_str().unwrap());

                for (sample, af) in reader.header().samples().into_iter().zip(afs.iter()) {
                    expr = expr.value(str::from_utf8(sample).unwrap(), af[0]);
                }
                assert!(
                    expr.exec()
                        .map(|v| v.as_bool().unwrap_or(false))
                        .unwrap_or(false),
                    "{:?} did not return true",
                    expr
                );
            }
        }

        if let Some(exprs) = self.yaml()["expected"]["posteriors"].as_vec() {
            for expr in exprs.iter() {
                let mut expr = Expr::new(expr.as_str().unwrap());

                for rec in reader.header().header_records() {
                    match rec {
                        bcf::HeaderRecord::Info { values, .. } => {
                            let id = values.get("ID").unwrap().clone();
                            if id.starts_with("PROB_") {
                                let values = call.info(id.as_bytes()).float().unwrap().unwrap();
                                expr = expr.value(id.clone(), values[0])
                            }
                        }
                        _ => (), // ignore other tags
                    }
                }
                assert!(
                    expr.exec()
                        .map(|v| v.as_bool().unwrap_or(false))
                        .unwrap_or(false),
                    "{:?} did not return true",
                    expr
                );
            }
        }
    }

    fn reference(ref_name: &str, ref_seq: &str) -> Result<NamedTempFile, Box<Error>> {
        let mut tmp_ref = tempfile::Builder::new().suffix(".fasta").tempfile()?;
        {
            let mut writer = fasta::Writer::new(&mut tmp_ref);
            writer.write(ref_name, None, ref_seq.as_bytes())?;
        }
        Command::new("samtools")
            .args(&["faidx", tmp_ref.path().to_str().unwrap()])
            .status()
            .expect("failed to create fasta index");

        Ok(tmp_ref)
    }

    fn alignment_properties(properties: &str) -> Result<NamedTempFile, Box<Error>> {
        let mut tmp_props = tempfile::Builder::new().suffix(".json").tempfile()?;
        tmp_props.as_file_mut().write_all(properties.as_bytes())?;

        Ok(tmp_props)
    }
}

macro_rules! testcase {
    ($name:ident) => {
        #[test]
        fn $name() {
            let name = stringify!($name);
            let testcase = Testcase::new(
                &Path::new(file!())
                    .parent()
                    .unwrap()
                    .join("resources/testcases")
                    .join(name),
            )
            .unwrap();
            testcase.run().unwrap();
            testcase.check();
        }
    };
}

testcase!(test01);
testcase!(test02);
testcase!(test03);
testcase!(test04);
testcase!(test05);
testcase!(test06);
testcase!(test07);
testcase!(test08);
testcase!(test09);
testcase!(test10);
testcase!(test11);
testcase!(test12);
testcase!(test13);
testcase!(test14);
testcase!(test15);
testcase!(test16);
testcase!(test17);
testcase!(test18);
testcase!(test19);
testcase!(test20);
// skip the next test because this insertion cannot currently be resolved properly
// TODO find a way to fix this.
// testcase!(test21);
testcase!(test22);
testcase!(test23);
testcase!(test24);
testcase!(test25);
testcase!(test26);
testcase!(test27);
testcase!(test28);
testcase!(test29);
testcase!(test30);
testcase!(test31);
testcase!(test32);
testcase!(test33);
testcase!(test34);
testcase!(test35);
testcase!(test36);
testcase!(pattern_too_long);
testcase!(test_wgbs01);
testcase!(test_long_pattern);

fn basedir(test: &str) -> String {
    format!("tests/resources/{}", test)
}

fn cleanup_file(f: &str) {
    if Path::new(f).exists() {
        fs::remove_file(f).unwrap();
    }
}

fn control_fdr(test: &str, event_str: &str, alpha: f64) {
    let basedir = basedir(test);
    let output = format!("{}/calls.filtered.bcf", basedir);
    cleanup_file(&output);
    varlociraptor::filtration::fdr::control_fdr(
        &format!("{}/calls.matched.bcf", basedir),
        Some(&output),
        &[varlociraptor::SimpleEvent {
            name: event_str.to_owned(),
        }],
        &varlociraptor::model::VariantType::Deletion(Some(1..30)),
        LogProb::from(Prob(alpha)),
    )
    .unwrap();
}

fn assert_call_number(test: &str, expected_calls: usize) {
    let basedir = basedir(test);

    let mut reader = bcf::Reader::from_path(format!("{}/calls.filtered.bcf", basedir)).unwrap();

    let calls = reader.records().map(|r| r.unwrap()).collect_vec();
    // allow one more or less, in order to be robust to numeric fluctuations
    assert!(
        (calls.len() as i32 - expected_calls as i32).abs() <= 1,
        "unexpected number of calls ({} vs {})",
        calls.len(),
        expected_calls
    );
}

#[test]
fn test_fdr_control1() {
    control_fdr("test_fdr_ev_1", "SOMATIC", 0.05);
    //assert_call_number("test_fdr_ev_1", 974);
}

#[test]
fn test_fdr_control2() {
    control_fdr("test_fdr_ev_2", "SOMATIC", 0.05);
    assert_call_number("test_fdr_ev_2", 985);
}

/// same test, but low alpha
#[test]
fn test_fdr_control3() {
    control_fdr("test_fdr_ev_3", "ABSENT", 0.001);
    assert_call_number("test_fdr_ev_3", 0);
}

#[test]
fn test_fdr_control4() {
    control_fdr("test_fdr_ev_4", "SOMATIC_TUMOR", 0.05);
    assert_call_number("test_fdr_ev_4", 0);
}

use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::fs;
use std::cmp;
use std::str;

use askama::Template;
use derive_builder::Builder;
use regex::Regex;
use rust_htslib::bam::Read as BamRead;
use rust_htslib::{bam, bcf, bcf::Read};
use bio::io::fasta;
use structopt::StructOpt;
use serde_json;
use serde::Serialize;

use crate::errors;
use crate::model::sample;
use crate::model::Variant;
use crate::utils;

lazy_static! {
    static ref TESTCASE_RE: Regex =
        Regex::new(r"^(?P<chrom>[^:]+):(?P<pos>\d+)(:(?P<idx>\d+))?$").unwrap();
}

#[derive(Template)]
#[template(path = "testcase.yml", escape = "none")]
struct TestcaseTemplate {
    samples: HashMap<String, Sample>,
    candidate: String,
    ref_name: String,
    ref_seq: String,
    options: String,
}

#[derive(Debug)]
struct Sample {
    path: String,
    properties: String,
}

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct Testcase<T>
where
    T: StructOpt
{
    #[builder(setter(into))]
    prefix: PathBuf,
    #[builder(private)]
    chrom_name: Vec<u8>,
    #[builder(private)]
    pos: u32,
    #[builder(private)]
    idx: usize,
    #[builder(private)]
    reference_reader: fasta::IndexedReader<File>,
    #[builder(private)]
    candidate_reader: bcf::Reader,
    #[builder(private)]
    bams: HashMap<String, PathBuf>,
    options: T
}

impl<T> TestcaseBuilder<T>
where
    T: StructOpt
{
    pub fn reference(self, path: impl AsRef<Path>) -> Result<Self, Box<Error>> {
        Ok(self.reference_reader(fasta::IndexedReader::from_file(&path)?))
    }

    pub fn candidates(self, path: impl AsRef<Path>) -> Result<Self, Box<Error>> {
        Ok(self.candidate_reader(bcf::Reader::from_path(path)?))
    }

    pub fn locus(self, locus: &str) -> Result<Self, Box<Error>> {
        if let Some(captures) = TESTCASE_RE.captures(locus) {
            let chrom_name = captures
                .name("chrom")
                .unwrap()
                .as_str()
                .as_bytes()
                .to_owned();
            let mut pos: u32 = captures.name("pos").unwrap().as_str().parse()?;
            pos -= 1;
            let idx = if let Some(m) = captures.name("idx") {
                let idx: usize = m.as_str().parse()?;
                idx - 1
            } else {
                0
            };
            Ok(self.chrom_name(chrom_name).pos(pos).idx(idx))
        } else {
            Err(errors::TestcaseError::InvalidLocus)?
        }
    }

    pub fn register_bam(mut self, name: &str, path: impl AsRef<Path>) -> Self {
        if self.bams.is_none() {
            self = self.bams(HashMap::new());
        }
        self.bams.as_mut().unwrap().insert(name.to_owned(), path.as_ref().to_owned());

        self
    }
}

impl<T> Testcase<T>
where
    T: StructOpt + Serialize
{
    fn variants(&mut self) -> Result<Vec<bcf::Record>, Box<Error>> {
        // get variant
        let rid = self.candidate_reader.header().name2rid(&self.chrom_name)?;
        let mut found = vec![];
        for res in self.candidate_reader.records() {
            let rec = res?;
            if let Some(rec_rid) = rec.rid() {
                if rec_rid == rid && rec.pos() == self.pos {
                    found.push(rec);
                }
            }
        }
        if found.len() == 0 {
            Err(errors::TestcaseError::NoCandidateFound)?
        } else {
            Ok(found)
        }
    }

    pub fn write(&mut self) -> Result<(), Box<Error>> {
        fs::create_dir_all(&self.prefix)?;

        let candidate_filename = Path::new("candidates.bcf");

        // get and write candidate
        let mut i = 0;
        let mut candidate = None;
        for mut record in self.variants()? {
            let variants = utils::collect_variants(&mut record, false, false, None)?;
            for variant in variants {
                if let Some(variant) = variant {
                    if i == self.idx {
                        candidate = Some((variant, record));
                        break;
                    }
                }
                i += 1;
            }
        }
        if candidate.is_none() {
            return Err(errors::TestcaseError::InvalidIndex)?;
        }
        let candidate = candidate.unwrap();

        let (start, end) = match candidate {
            (Variant::Deletion(l), _) => (self.pos.saturating_sub(1000), self.pos + l + 1000),
            (Variant::Insertion(ref seq), _) => (
                self.pos.saturating_sub(1000),
                self.pos + seq.len() as u32 + 1000,
            ),
            (Variant::SNV(_), _) => (self.pos.saturating_sub(100), self.pos + 1 + 100),
            (Variant::None, _) => (self.pos.saturating_sub(100), self.pos + 1 + 100),
        };

        let mut ref_start = start;
        let mut ref_end = end;
        // first pass, extend reference interval
        for path in self.bams.values() {
            let mut bam_reader = bam::IndexedReader::from_path(path)?;
            let tid = bam_reader.header().tid(&self.chrom_name).unwrap();
            bam_reader.fetch(tid, start, end)?;
            for res in bam_reader.records() {
                let rec = res?;
                ref_start = cmp::min(rec.pos() as u32, ref_start);
                ref_end = cmp::max(rec.cigar().end_pos()? as u32, ref_end);
            }
        }

        // second pass, write samples
        let mut samples = HashMap::new();
        for (name, path) in &self.bams {
            let properties = sample::estimate_alignment_properties(path)?;
            let mut bam_reader = bam::IndexedReader::from_path(path)?;
            let filename = Path::new(name).with_extension("bam");
            let mut bam_writer = bam::Writer::from_path(
                self.prefix.join(&filename),
                &bam::Header::from_template(bam_reader.header()),
            )?;
            let tid = bam_reader.header().tid(&self.chrom_name).unwrap();

            bam_reader.fetch(tid, start, end)?;
            for res in bam_reader.records() {
                let mut rec = res?;
                // update mapping position to interval
                rec.set_pos(rec.pos() - ref_start as i32);
                bam_writer.write(&rec)?;
            }
            samples.insert(
                name.to_owned(),
                Sample {
                    path: filename.to_str().unwrap().to_owned(),
                    properties: serde_json::to_string(&properties)?,
                },
            );
        }

        // write candidate
        let mut candidate_writer = bcf::Writer::from_path(
            self.prefix.join(candidate_filename),
            &bcf::Header::from_template(self.candidate_reader.header()),
            false,
            false,
        )?;
        let (_, mut candidate_record) = candidate;
        candidate_record.set_pos((candidate_record.pos() - ref_start) as i32);
        candidate_writer.write(&candidate_record)?;

        // fetch reference
        let ref_name = str::from_utf8(&self.chrom_name)?;
        self.reference_reader.fetch(ref_name, ref_start as u64, ref_end as u64)?;
        let mut ref_seq = Vec::new();
        self.reference_reader.read(&mut ref_seq)?;

        let mut desc = File::create(self.prefix.join("testcase.yaml"))?;
        desc.write_all(
            TestcaseTemplate {
                samples,
                options: serde_json::to_string(&self.options)?,
                candidate: candidate_filename.to_str().unwrap().to_owned(),
                ref_seq: String::from_utf8(ref_seq)?.to_owned(),
                ref_name: ref_name.to_owned(),
            }
            .render()?
            .as_bytes(),
        )?;

        Ok(())
    }
}

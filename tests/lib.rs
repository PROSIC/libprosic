use std::error::Error;
use std::path::{Path, PathBuf};
use std::fs::File;
use std::io::Read;
use std::process::Command;

use rust_htslib::bam;
use yaml_rust::{YamlLoader, Yaml};
use serde_json;
use tempfile::{self, NamedTempFile};
use bio::io::fasta;

use varlociraptor::cli::{Varlociraptor, run};

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

    fn run(&self) -> Result<(), Box<Error>>{
        dbg!(self.path.join("calls.bcf"));
        let mut options = serde_json::from_str(self.yaml()["options"].as_str().unwrap())?;
        // TODO alignment properties!
        match &mut options {
            Varlociraptor::CallTumorNormal { ref mut reference, ref mut tumor, ref mut normal, ref mut candidates, ref mut output, ref mut testcase_locus, ref mut testcase_prefix, .. } => {
                let temp_ref = Self::reference(self.yaml()["reference"]["name"].as_str().unwrap(), self.yaml()["reference"]["seq"].as_str().unwrap())?;
                *reference = temp_ref.path().to_owned();
                dbg!(reference);
                *tumor = self.path.join(self.yaml()["samples"]["tumor"]["path"].as_str().unwrap());
                *normal = self.path.join(self.yaml()["samples"]["normal"]["path"].as_str().unwrap());
                *candidates = Some(self.path.join(self.yaml()["candidate"].as_str().unwrap()));
                *output = Some(self.path.join("calls.bcf"));
                *testcase_prefix = None;
                *testcase_locus = None;

                bam::index::build(tumor, None, bam::index::Type::BAI, 1).unwrap();
                bam::index::build(normal, None, bam::index::Type::BAI, 1).unwrap();

                run(options)
            },
            _ => panic!("unsupported subcommand")
        }

        // check results

    }

    fn reference(ref_name: &str, ref_seq: &str) -> Result<NamedTempFile, Box<Error>> {
        let mut tmp_ref = tempfile::Builder::new().suffix(".fasta").tempfile()?;
        {
            dbg!(&tmp_ref);
            let mut writer = fasta::Writer::new(&mut tmp_ref);
            writer.write(ref_name, None, ref_seq.as_bytes())?;
        }
        Command::new("samtools")
            .args(&["faidx", tmp_ref.path().to_str().unwrap()])
            .status()
            .expect("failed to create fasta index");

        Ok(tmp_ref)
    }
}

macro_rules! testcase {
    ($name:ident) => {
        #[test]
        fn $name() {
            let name = stringify!($name);
            let testcase = Testcase::new(&Path::new(file!()).parent().unwrap().join("resources/testcases").join(name)).unwrap();
            testcase.run().unwrap();
        }
    }
}

testcase!(test01);

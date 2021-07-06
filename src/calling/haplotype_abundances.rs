use anyhow::Result;
use derive_builder::Builder;
use hdf5;
use kernel_density;
use rust_htslib::bcf;
use serde_json::json;
use std::fs::File;
use std::path::PathBuf;

#[derive(Builder)]
#[builder(pattern = "owned")]
pub(crate) struct Caller {
    hdf5_reader: hdf5::File,
    vcf_reader: bcf::Reader,
    min_norm_counts: f64,
}

#[derive(Serialize)]
pub(crate) struct ECDF {
    seqname: String,
    bootstraps: Vec<u64>,
    ecdf: Vec<f64>,
}

impl Caller {
    //Below function outputs a vector of filtered seqnames according to --min-norm-counts.
    pub(crate) fn filter_seqnames(&self) -> Result<Vec<String>> {
        let hdf5 = &self.hdf5_reader;
        let est_counts = hdf5.dataset("est_counts")?.read_1d::<f64>()?;
        let seq_length = hdf5.dataset("aux/lengths")?.read_1d::<f64>()?; //these two variables arrays have the same length.

        let norm_counts = est_counts / seq_length;
        let mut indices = Vec::new();
        for (i, num) in norm_counts.iter().enumerate() {
            if num > &self.min_norm_counts {
                indices.push(i);
            }
        }
        //let len = hdf5.dataset("aux/ids")?.dtype().unwrap().size();
        let ids = hdf5.dataset("aux/ids")?.read_1d::<hdf5::types::FixedAscii<255>>()?;
    
        let mut filtered: Vec<String> = Vec::new();
        for i in indices {
            filtered.push(ids[i].to_string());
        }
        Ok(filtered)
    }

    //Below function outputs an ECDF struct for one haplotype.
    pub(crate) fn cdf(&self, seqname: String) -> Result<ECDF> {
        let hdf5 = &self.hdf5_reader;
        let ids = hdf5
            .dataset("aux/ids")?
            .read_1d::<hdf5::types::FixedAscii<255>>()?;
        let index = ids.iter().position(|x| x.as_str() == seqname).unwrap();
        let num_bootstraps = hdf5.dataset("aux/num_bootstrap")?.read_1d::<i32>()?;
        let seq_length = hdf5.dataset("aux/lengths")?.read_1d::<u64>()?;
        let mut bootstraps = Vec::new();
        for i in 0..num_bootstraps[0] {
            let dataset = hdf5.dataset(&format!("bootstrap/bs{i}", i = i))?;
            let est_counts = dataset.read_1d::<u64>()?;
            let norm_counts = est_counts / &seq_length;
            let norm_counts = norm_counts[index];
            bootstraps.push(norm_counts);
        }
        let ecdf = kernel_density::ecdf::Ecdf::new(&bootstraps);
        let ecdf: Vec<f64> = bootstraps.iter().map(|&x| ecdf.value(x)).collect();
        Ok(ECDF {
            seqname,
            bootstraps,
            ecdf,
        })
    }

    pub(crate) fn call(&self) {
        todo!()
    }
}

impl ECDF {
    pub(crate) fn plot_qc(&self, mut qc_plot: PathBuf) -> std::io::Result<()> {
        let json = include_str!("../../templates/plots/qc_cdf.json");
        let mut blueprint: serde_json::Value = serde_json::from_str(json)?;
        let bootstraps = &self.bootstraps;
        let ecdf = &self.ecdf;

        let data = json!([bootstraps, ecdf]); //need to play with the values to have discrete x and y points
    
        blueprint["data"]["values"] = data;
        println!{"{:?}", blueprint};

        qc_plot.push(self.seqname.clone() + &".json".to_string());
        let file = File::create(qc_plot)?;
        serde_json::to_writer(file, &blueprint);
        Ok(())
    }
}

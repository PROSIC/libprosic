use std::str;
use std::collections::{HashMap, VecDeque, vec_deque};
use std::marker::PhantomData;
use std::cmp;

use rgsl::randist::gaussian::gaussian_pdf;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::record::Cigar;
use bio::stats::{logprobs, LogProb};

use model::priors::AlleleFreq;
use model;
use model::Variant;


pub fn prob_mapping(mapq: u8) -> LogProb {
    logprobs::ln_1m_exp(logprobs::phred_to_log(mapq as f64))
}


/// Ringbuffer of BAM records. This data structure ensures that no bam record is read twice while
/// extracting observations for given variants.
pub struct RecordBuffer {
    reader: bam::IndexedReader,
    inner: VecDeque<bam::Record>,
    window: u32
}


unsafe impl Sync for RecordBuffer {}
unsafe impl Send for RecordBuffer {}


impl RecordBuffer {
    /// Create a new `RecordBuffer`.
    pub fn new(bam: bam::IndexedReader, window: u32) -> Self {
        RecordBuffer {
            reader: bam,
            inner: VecDeque::with_capacity(window as usize * 2),
            window: window as u32
        }
    }

    /// Return end position of buffer.
    fn end(&self) -> Option<u32> {
        self.inner.back().map(|rec| rec.pos() as u32)
    }

    /// Fill buffer around the given interval.
    pub fn fill(&mut self, chrom: &[u8], start: u32, end: u32) -> Result<(), String> {
        if let Some(tid) = self.reader.header.tid(chrom) {
            let window_start = cmp::max(start as i32 - self.window as i32 - 1, 0) as u32;
            if self.inner.is_empty() || self.end().unwrap() < window_start {
                let end = self.reader.header.target_len(tid).unwrap();
                if self.reader.seek(tid, window_start, end).is_err() {
                    return Err(format!("Unable to seek to variant at {}:{}", str::from_utf8(chrom).unwrap(), start));
                }
                self.inner.clear();
            } else {
                // remove records too far left
                let to_remove = self.inner.iter().take_while(|rec| rec.pos() < window_start as i32).count();
                self.inner.drain(..to_remove);
            }

            // extend to the right
            for record in self.reader.records() {
                let record = match record {
                    Err(_) => return Err("Error reading BAM record.".to_owned()),
                    Ok(rec) => rec
                };
                let pos = record.pos();
                if record.is_duplicate() || record.is_unmapped() {
                    continue;
                }
                self.inner.push_back(record);
                if pos > end as i32 + self.window as i32 {
                    break;
                }
            }

            Ok(())
        } else {
            Err(format!("Sequence {} cannot be found in BAM", str::from_utf8(chrom).unwrap()))
        }
    }

    pub fn iter(&self) -> vec_deque::Iter<bam::Record> {
        self.inner.iter()
    }
}


/// An observation for or against a variant.
#[derive(Clone, Debug)]
pub struct Observation {
    /// Posterior probability that the read/read-pair has been mapped correctly (1 - MAPQ).
    pub prob_mapping: LogProb,
    /// Probability that the read/read-pair comes from the alternative allele.
    pub prob_alt: LogProb,
    /// Probability that the read/read-pair comes from the reference allele.
    pub prob_ref: LogProb,
    /// Probability of the read/read-pair given that it has been mismapped.
    pub prob_mismapped: LogProb
}


/// Expected insert size in terms of mean and standard deviation.
/// This should be estimated from unsorted(!) bam files to avoid positional biases.
#[derive(Copy, Clone, Debug)]
pub struct InsertSize {
    pub mean: f64,
    pub sd: f64
}


/// A sequenced sample, e.g., a tumor or a normal sample.
pub struct Sample<A: AlleleFreq, P: model::priors::Model<A>> {
    record_buffer: RecordBuffer,
    insert_size: InsertSize,
    prior_model: P,
    likelihood_model: model::likelihood::LatentVariableModel,
    a: PhantomData<A>
}


impl<A: AlleleFreq, P: model::priors::Model<A>> Sample<A, P> {
    /// Create a new `Sample`.
    ///
    /// # Arguments
    ///
    /// * `bam` - BAM file with the aligned and deduplicated sequence reads.
    /// * `pileup_window` - Window around the variant that shall be search for evidence (e.g. 5000).
    /// * `insert_size` - estimated insert size
    /// * `prior_model` - Prior assumptions about allele frequency spectrum of this sample.
    /// * `likelihood_model` - Latent variable model to calculate likelihoods of given observations.
    pub fn new(bam: bam::IndexedReader, pileup_window: u32, insert_size: InsertSize, prior_model: P, likelihood_model: model::likelihood::LatentVariableModel) -> Self {
        Sample {
            record_buffer: RecordBuffer::new(bam, pileup_window),
            insert_size: insert_size,
            prior_model: prior_model,
            likelihood_model: likelihood_model,
            a: PhantomData
        }
    }

    /// Return prior model.
    pub fn prior_model(&self) -> &P {
        &self.prior_model
    }

    /// Return likelihood model.
    pub fn likelihood_model(&self) -> &model::likelihood::LatentVariableModel {
        &self.likelihood_model
    }

    /// Extract observations for the given variant.
    pub fn extract_observations(&mut self, chrom: &[u8], start: u32, variant: Variant) -> Result<Vec<Observation>, String> {
        let mut observations = Vec::new();
        let (end, varpos) = match variant {
            Variant::Deletion(length)  => (start + length, (start + length / 2) as i32), // TODO do we really need two centerpoints?
            Variant::Insertion(length) => (start + length, start as i32)
        };
        let mut pairs = HashMap::new();


        // move window to the current variant
        debug!("Filling buffer...");
        try!(self.record_buffer.fill(chrom, start, end));
        debug!("Done.");

        // iterate over records
        for record in self.record_buffer.iter() {
            let pos = record.pos();
            let cigar = record.cigar();
            let end_pos = record.end_pos(&cigar);
            if pos <= varpos && end_pos >= varpos {
                // overlapping alignment
                observations.push(self.read_observation(&record, &cigar, varpos, variant));
            } else if end_pos <= varpos {
                // need to check mate
                // since the bam file is sorted by position, we can't see the mate first
                if record.mpos() >= varpos {
                    pairs.insert(record.qname().to_owned(), record.mapq());
                }
            } else if let Some(mate_mapq) = pairs.get(record.qname()) {
                // mate already visited, and this read maps right of varpos
                observations.push(self.fragment_observation(&record, *mate_mapq, variant));
            }
        }
        debug!("Extracted observations.");
        debug!("{} pairs", pairs.len());
        Ok(observations)
    }

    fn read_observation(&self, record: &bam::Record, cigar: &[Cigar], varpos: i32, variant: Variant) -> Observation {
        let (length, is_del) = match variant {
            Variant::Deletion(length) => (length, true),
            Variant::Insertion(length) => (length, false)
        };
        let is_close = |qpos: i32, qlength: u32| (varpos - (qpos + qlength as i32 / 2)).abs() <= 50;
        let is_similar_length = |qlength: u32| (length as i32 - qlength as i32).abs() <= 20;

        let mut qpos = record.pos();
        let prob_mapping = prob_mapping(record.mapq());
        let mut obs = Observation {
            prob_mapping: prob_mapping,
            prob_alt: 1.0f64.ln(),
            prob_ref: 0.0f64.ln(),
            prob_mismapped: 0.0 // if the read is mismapped, we assume sampling probability 1.0
        };
        for c in cigar {
            match c {
                &Cigar::Match(l) | &Cigar::RefSkip(l) | &Cigar::Equal(l) | &Cigar::Diff(l) => qpos += l as i32,
                &Cigar::Back(l) => qpos -= l as i32,
                &Cigar::Del(l) if is_del && is_similar_length(l) && is_close(qpos, l) => {
                    // supports alt allele
                    return obs;
                },
                &Cigar::Ins(l) if !is_del && is_similar_length(l) && is_close(qpos, l) => {
                    // supports alt allele
                    return obs;
                }
                _ => ()
            }
        }

        // support ref allele
        obs.prob_alt = 0.0f64.ln();
        obs.prob_ref = 1.0f64.ln();
        obs
    }

    fn fragment_observation(&self, record: &bam::Record, mate_mapq: u8, variant: Variant) -> Observation {
        let insert_size = record.insert_size();
        let shift = match variant {
            Variant::Deletion(length)  => length as f64,
            Variant::Insertion(length) => -(length as f64)
        };
        let obs = Observation {
            prob_mapping: prob_mapping(record.mapq()) + prob_mapping(mate_mapq),
            prob_alt: gaussian_pdf(insert_size as f64 - self.insert_size.mean, self.insert_size.sd).ln(),
            prob_ref: gaussian_pdf(insert_size as f64 - self.insert_size.mean + shift, self.insert_size.sd).ln(),
            prob_mismapped: 0.0 // if the fragment is mismapped, we assume sampling probability 1.0
        };
        debug!("{:?}", obs);

        obs
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prob_mapping() {
        assert_relative_eq!(prob_mapping(0).exp(), 0.0);
        assert_relative_eq!(prob_mapping(10).exp(), 0.9);
        assert_relative_eq!(prob_mapping(20).exp(), 0.99);
        assert_relative_eq!(prob_mapping(30).exp(), 0.999);
    }
}

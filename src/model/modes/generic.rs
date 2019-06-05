use std::cmp;

use bio::stats::bayesian::model::{Likelihood, Model, Posterior, Prior};
use bio::stats::LogProb;
use itertools::Itertools;
use vec_map::VecMap;

use crate::grammar;
use crate::model;
use crate::model::likelihood;
use crate::model::sample::Pileup;
use crate::model::{AlleleFreq, Contamination, StrandBias};

#[derive(Debug)]
pub enum CacheEntry {
    ContaminatedSample(likelihood::ContaminatedSampleCache),
    SingleSample(likelihood::SingleSampleCache),
}

impl CacheEntry {
    fn new(contaminated: bool) -> Self {
        if contaminated {
            CacheEntry::ContaminatedSample(likelihood::ContaminatedSampleCache::default())
        } else {
            CacheEntry::SingleSample(likelihood::SingleSampleCache::default())
        }
    }
}

pub type Cache = VecMap<CacheEntry>;

#[derive(Default, Debug, Clone)]
pub struct GenericModelBuilder<P> {
    resolutions: Vec<usize>,
    contaminations: Vec<Option<Contamination>>,
    prior: P,
}

impl<P: Prior> GenericModelBuilder<P>
where
    P: Prior<Event = Vec<likelihood::Event>>,
{
    pub fn push_sample(mut self, resolution: usize, contamination: Option<Contamination>) -> Self {
        self.contaminations.push(contamination);
        self.resolutions.push(resolution);

        self
    }

    pub fn prior(mut self, prior: P) -> Self {
        self.prior = prior;

        self
    }

    pub fn build(self) -> Result<Model<GenericLikelihood, P, GenericPosterior, Cache>, String> {
        let posterior = GenericPosterior::new(self.resolutions);
        let likelihood = GenericLikelihood::new(self.contaminations);
        Ok(Model::new(likelihood, self.prior, posterior))
    }
}

#[derive(new, Default, Clone, Debug)]
pub struct GenericPosterior {
    resolutions: Vec<usize>,
}

impl GenericPosterior {
    fn grid_points(&self, pileups: &[Pileup]) -> Vec<usize> {
        pileups
            .iter()
            .zip(self.resolutions.iter())
            .map(|(pileup, res)| {
                let n_obs = pileup.len();
                let mut n = cmp::min(cmp::max(n_obs + 1, 5), *res);
                if n % 2 == 0 {
                    n += 1;
                }
                n
            })
            .collect()
    }

    fn density<F: FnMut(&<Self as Posterior>::BaseEvent, &<Self as Posterior>::Data) -> LogProb>(
        &self,
        vaf_tree_node: &grammar::vaftree::Node,
        base_events: &mut VecMap<likelihood::Event>,
        sample_grid_points: &[usize],
        pileups: &<Self as Posterior>::Data,
        strand_bias: StrandBias,
        joint_prob: &mut F,
    ) -> LogProb {
        let sample = *vaf_tree_node.sample();
        let mut subdensity = |base_events: &mut VecMap<likelihood::Event>| {
            if vaf_tree_node.is_leaf() {
                joint_prob(&base_events.values().cloned().collect(), pileups)
            } else {
                if vaf_tree_node.is_branching() {
                    LogProb::ln_sum_exp(
                        &vaf_tree_node
                            .children()
                            .iter()
                            .map(|child| {
                                self.density(
                                    child,
                                    &mut base_events.clone(),
                                    sample_grid_points,
                                    pileups,
                                    strand_bias,
                                    joint_prob,
                                )
                            })
                            .collect_vec(),
                    )
                } else {
                    self.density(
                        &vaf_tree_node.children()[0],
                        base_events,
                        sample_grid_points,
                        pileups,
                        strand_bias,
                        joint_prob,
                    )
                }
            }
        };

        let push_base_event = |allele_freq, base_events: &mut VecMap<likelihood::Event>| {
            base_events.insert(
                sample,
                likelihood::Event {
                    allele_freq: allele_freq,
                    strand_bias: strand_bias,
                },
            );
        };

        match vaf_tree_node.vafs() {
            grammar::VAFSpectrum::Set(vafs) => {
                if vafs.len() == 1 {
                    push_base_event(vafs.iter().next().unwrap().clone(), base_events);
                    subdensity(base_events)
                } else {
                    LogProb::ln_sum_exp(
                        &vafs
                            .iter()
                            .map(|vaf| {
                                let mut base_events = base_events.clone();
                                push_base_event(*vaf, &mut base_events);
                                subdensity(&mut base_events)
                            })
                            .collect_vec(),
                    )
                }
            }
            grammar::VAFSpectrum::Range(vafs) => {
                let n_obs = pileups[sample].len();
                LogProb::ln_simpsons_integrate_exp(
                    |_, vaf| {
                        let mut base_events = base_events.clone();
                        push_base_event(AlleleFreq(vaf), &mut base_events);
                        subdensity(&mut base_events)
                    },
                    *vafs.observable_min(n_obs),
                    *vafs.observable_max(n_obs),
                    sample_grid_points[sample],
                )
            }
        }
    }
}

impl Posterior for GenericPosterior {
    type BaseEvent = Vec<likelihood::Event>;
    type Event = model::Event;
    type Data = Vec<Pileup>;

    fn compute<F: FnMut(&Self::BaseEvent, &Self::Data) -> LogProb>(
        &self,
        event: &Self::Event,
        pileups: &Self::Data,
        joint_prob: &mut F,
    ) -> LogProb {
        let grid_points = self.grid_points(pileups);
        let vaf_tree = &event.vafs;
        LogProb::ln_sum_exp(
            &vaf_tree
                .iter()
                .map(|node| {
                    let mut base_events = VecMap::with_capacity(pileups.len());
                    self.density(
                        node,
                        &mut base_events,
                        &grid_points,
                        pileups,
                        event.strand_bias,
                        joint_prob,
                    )
                })
                .collect_vec(),
        )
    }
}

#[derive(Clone, Debug)]
enum SampleModel {
    Contaminated {
        likelihood_model: likelihood::ContaminatedSampleLikelihoodModel,
        by: usize,
    },
    Normal(likelihood::SampleLikelihoodModel),
}

#[derive(Default, Clone, Debug)]
pub struct GenericLikelihood {
    inner: Vec<SampleModel>,
}

impl GenericLikelihood {
    pub fn new(contaminations: Vec<Option<Contamination>>) -> Self {
        let mut inner = Vec::new();
        for contamination in contaminations.iter() {
            if let Some(contamination) = contamination {
                inner.push(SampleModel::Contaminated {
                    likelihood_model: likelihood::ContaminatedSampleLikelihoodModel::new(
                        1.0 - contamination.fraction,
                    ),
                    by: contamination.by,
                });
            } else {
                inner.push(SampleModel::Normal(likelihood::SampleLikelihoodModel::new()));
            }
        }
        GenericLikelihood { inner }
    }
}

impl Likelihood<Cache> for GenericLikelihood {
    type Event = Vec<likelihood::Event>;
    type Data = Vec<Pileup>;

    fn compute(&self, events: &Self::Event, pileups: &Self::Data, cache: &mut Cache) -> LogProb {
        let mut p = LogProb::ln_one();

        for (((sample, event), pileup), inner) in events
            .iter()
            .enumerate()
            .zip(pileups.iter())
            .zip(self.inner.iter())
        {
            p += match inner {
                &SampleModel::Contaminated {
                    ref likelihood_model,
                    by,
                } => {
                    if let CacheEntry::ContaminatedSample(ref mut cache) =
                        cache.entry(sample).or_insert_with(|| CacheEntry::new(true))
                    {
                        likelihood_model.compute(
                            &likelihood::ContaminatedSampleEvent {
                                primary: event.clone(),
                                secondary: events[by].clone(),
                            },
                            pileup,
                            cache,
                        )
                    } else {
                        unreachable!();
                    }
                }
                &SampleModel::Normal(ref likelihood_model) => {
                    if let CacheEntry::SingleSample(ref mut cache) = cache
                        .entry(sample)
                        .or_insert_with(|| CacheEntry::new(false))
                    {
                        likelihood_model.compute(event, pileup, cache)
                    } else {
                        unreachable!();
                    }
                }
            }
        }

        p
    }
}

#[derive(Default, Clone, Debug)]
pub struct FlatPrior {}

impl FlatPrior {
    pub fn new() -> Self {
        FlatPrior {}
    }
}

impl Prior for FlatPrior {
    type Event = Vec<likelihood::Event>;

    fn compute(&self, _event: &Self::Event) -> LogProb {
        LogProb::ln_one()
    }
}

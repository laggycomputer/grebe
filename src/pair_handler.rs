use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::io;
use std::io::BufWriter;

use bio::bio_types::sequence::SequenceRead;
use bio::io::fastq;
use clap::builder::PossibleValue;
use clap::ValueEnum;
use strum::VariantArray;
use strum_macros::VariantArray;

use crate::types::{FastqPair, UMIVec};
use crate::writer::WriterMaybeGzip;

#[derive(Clone, Copy, PartialEq, VariantArray)]
pub(crate) enum UMICollisionResolutionMethod {
    None,
    KeepFirst,
    KeepLast,
    KeepLongestLeft,
    KeepLongestRight,
    KeepLongestExtend,
    QualityVote,
}

impl UMICollisionResolutionMethod {
    fn _compare_for_extension(&self, old: &fastq::Record, new: &fastq::Record) -> fastq::Record {
        match self {
            UMICollisionResolutionMethod::KeepLongestLeft | UMICollisionResolutionMethod::KeepLongestRight |
            UMICollisionResolutionMethod::KeepLongestExtend => {
                match old.seq().len().cmp(&new.seq().len()) {
                    Ordering::Less => match self {
                        // if keeping longest with extension, need to check content
                        UMICollisionResolutionMethod::KeepLongestExtend => {
                            if new.seq().starts_with(old.seq()) { (*new).clone() } else { (*old).clone() }
                        }
                        // if simply keeping longest, content does not matter
                        _ => (*new).clone()
                    },
                    Ordering::Equal => match &self {
                        // maybe allow user to choose here for KeepLongestExtend
                        UMICollisionResolutionMethod::KeepLongestLeft |
                        UMICollisionResolutionMethod::KeepLongestExtend => (*old).clone(),
                        UMICollisionResolutionMethod::KeepLongestRight => (*new).clone(),
                        _ => unreachable!()
                    },
                    Ordering::Greater => (*old).clone()
                }
            }
            _ => unimplemented!()
        }
    }
}

impl ValueEnum for UMICollisionResolutionMethod {
    fn value_variants<'a>() -> &'a [Self] {
        &Self::VARIANTS
    }

    fn to_possible_value(&self) -> Option<PossibleValue> {
        Some(match self {
            Self::None => PossibleValue::new("none")
                .help("keep duplicates, prepend their assigned UMI to their sequence names"),
            UMICollisionResolutionMethod::KeepFirst => PossibleValue::new("keep-first")
                .help("keep the first sequence matched to this UMI, ignore any sequences that follow"),
            UMICollisionResolutionMethod::KeepLast => PossibleValue::new("keep-last")
                .help("keep the last sequence matched, ignore any sequences that come before"),
            UMICollisionResolutionMethod::KeepLongestLeft => PossibleValue::new("keep-longest-left")
                .help("keep the longest sequence matched, favor the earlier sequence when tied"),
            UMICollisionResolutionMethod::KeepLongestRight => PossibleValue::new("keep-longest-right")
                .help("keep the longest sequence matched, favor the later sequence when tied"),
            UMICollisionResolutionMethod::KeepLongestExtend => PossibleValue::new("keep-longest-extend")
                .help("keep the longest sequence, overwrite it if a longer, later read agrees completely"),
            UMICollisionResolutionMethod::QualityVote => PossibleValue::new("quality-vote")
                .help("create one final sequence by combining base calls and qualities from all matched reads"),
        })
    }
}


pub(crate) struct PairHandler {
    pub(crate) record_writers: (fastq::Writer<WriterMaybeGzip>, fastq::Writer<WriterMaybeGzip>),
    pub(crate) collision_resolution_method: UMICollisionResolutionMethod,
    pub(crate) umi_bins: HashMap<UMIVec, HashSet<FastqPair>>,
    pub(crate) total_records: usize,
    pub(crate) good_records: usize,
    // ATCG
    pub(crate) quality_votes: Option<HashMap<UMIVec, (usize, usize, usize, usize)>>,
}

impl Default for PairHandler {
    fn default() -> Self {
        PairHandler {
            record_writers: (
                fastq::Writer::from_bufwriter(BufWriter::new(WriterMaybeGzip::NULL(io::sink()))),
                fastq::Writer::from_bufwriter(BufWriter::new(WriterMaybeGzip::NULL(io::sink())))),
            collision_resolution_method: UMICollisionResolutionMethod::KeepFirst,
            umi_bins: Default::default(),
            total_records: 0,
            good_records: 0,
            quality_votes: Default::default(),
        }
    }
}

impl PairHandler {
    pub(crate) fn write_pair(&mut self, pair: FastqPair) {
        // TODO: reimplement slicing etc
        self.record_writers.0.write(
            std::str::from_utf8(pair.0.name()).unwrap(),
            Option::from(pair.0.id()),
            pair.0.seq(),
            pair.0.qual(),
        )
            .expect("couldn't write out a forward record");
        self.record_writers.1.write(
            std::str::from_utf8(pair.1.name()).unwrap(),
            Option::from(pair.1.id()),
            pair.1.seq(),
            pair.1.qual(),
        ).expect("couldn't write out a reverse record");
    }

    pub(crate) fn insert_pair(&mut self, umi: &UMIVec, pair: &FastqPair) {
        match self.collision_resolution_method {
            UMICollisionResolutionMethod::None => {
                // special case: no comparison, etc., just go straight to disk
                self.good_records += 1;

                // write the record, add UMI
                let id_prefix = match umi.len() {
                    0 => "",
                    _ => std::str::from_utf8(umi).unwrap()
                };
                let pair_new = (
                    fastq::Record::with_attrs(
                        &*(id_prefix.to_owned() + " " + std::str::from_utf8(pair.0.name()).unwrap()),
                        pair.0.desc(),
                        &*pair.0.seq(),
                        &*pair.0.qual(),
                    ),
                    fastq::Record::with_attrs(
                        &*(id_prefix.to_owned() + " " + std::str::from_utf8(pair.1.name()).unwrap()),
                        pair.1.desc(),
                        &*pair.1.seq(),
                        &*pair.1.qual(),
                    )
                );
                self.write_pair(pair_new);
            }
            _ => {
                if !self.umi_bins.contains_key(umi) {
                    let mut set = HashSet::<FastqPair>::new();
                    match self.collision_resolution_method {
                        UMICollisionResolutionMethod::None => {
                            // handled above, but technically still a case here
                            unreachable!();
                        }
                        UMICollisionResolutionMethod::KeepFirst => {
                            // write the record immediately; save memory
                            self.write_pair(pair.clone());
                            // save an empty set so we don't come here again
                        }
                        UMICollisionResolutionMethod::QualityVote => {
                            // submit the "ballots" and save to disk later
                            // here, an empty set will also be used to take note of the UMI bin
                            todo!();
                        }
                        _ => {
                            // otherwise, we need to save this
                            set.insert(pair.clone());
                        }
                    }
                    // count the good record once for all these cases it could be overridden but the +1 will not change
                    self.good_records += 1;
                    // whatever we decided to put in the set (if anything), save it
                    self.umi_bins.insert(umi.clone(), set);
                } else {
                    let set = self.umi_bins.get_mut(umi).unwrap();

                    match self.collision_resolution_method {
                        UMICollisionResolutionMethod::None | UMICollisionResolutionMethod::QualityVote => {
                            // just handle it later somehow
                            set.insert(pair.clone());
                        }
                        UMICollisionResolutionMethod::KeepFirst => {}  // drop the new record
                        UMICollisionResolutionMethod::KeepLast => {
                            set.clear();
                            set.insert(pair.clone());
                        }
                        UMICollisionResolutionMethod::KeepLongestLeft | UMICollisionResolutionMethod::KeepLongestRight |
                        UMICollisionResolutionMethod::KeepLongestExtend => {
                            // slightly inefficient; but old must survive the clear
                            // todo: find a way to not copy this memory
                            let old = set.iter().next().unwrap().clone();

                            set.clear();
                            set.insert((
                                self.collision_resolution_method._compare_for_extension(&old.0, &pair.0),
                                self.collision_resolution_method._compare_for_extension(&old.1, &pair.1)
                            ));
                        }
                    }
                }
            }
        }
    }

    pub(crate) fn save_remaining(&mut self) {
        for (umi, pairs) in
        <HashMap<UMIVec, HashSet<(fastq::Record, fastq::Record)>> as Clone>::clone(&self.umi_bins).into_iter() {
            match self.collision_resolution_method {
                UMICollisionResolutionMethod::KeepFirst | UMICollisionResolutionMethod::None => {
                    // these records are already on disk
                }
                UMICollisionResolutionMethod::QualityVote => {
                    todo!()
                }
                _ => {
                    // conflict resolution has already selected a single read
                    self.write_pair(pairs.iter().next().unwrap().clone());
                }
            };
        }
    }
}
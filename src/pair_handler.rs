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

use crate::FastqPair;
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

    // TODO: respect these
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
    pub(crate) umi_bins: HashMap<Vec<u8>, HashSet<FastqPair>>,
    total_records: usize,
    good_records: usize,
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

    pub(crate) fn handle_pair(
        &mut self, umi: &Vec<u8>, new: &FastqPair) {
        if !self.umi_bins.contains_key(umi) {
            let mut set = HashSet::<FastqPair>::new();
            if self.collision_resolution_method == UMICollisionResolutionMethod::KeepFirst {
                // write the record immediately; save memory
                self.write_pair(new.clone());
                // save an empty set so we don't come here again
                self.umi_bins.insert(umi.clone(), set);
                return;
            } else {
                // otherwise, we need to save this
                set.insert(new.clone());
            }
            self.umi_bins.insert(umi.clone(), set);
            return;
        }

        let set = self.umi_bins.get_mut(umi).unwrap();

        match self.collision_resolution_method {
            UMICollisionResolutionMethod::None | UMICollisionResolutionMethod::QualityVote => {
                // just handle it later somehow
                set.insert(new.clone());
            }
            UMICollisionResolutionMethod::KeepFirst => {}  // drop the new record
            UMICollisionResolutionMethod::KeepLast => {
                set.clear();
                set.insert(new.clone());
            }
            UMICollisionResolutionMethod::KeepLongestLeft | UMICollisionResolutionMethod::KeepLongestRight |
            UMICollisionResolutionMethod::KeepLongestExtend => {
                // slightly inefficient; but old must survive the clear
                // todo: find a way to not copy this memory
                let old = set.iter().next().unwrap().clone();

                set.clear();
                set.insert((
                    self.collision_resolution_method._compare_for_extension(&old.0, &new.0),
                    self.collision_resolution_method._compare_for_extension(&old.1, &new.1)
                ));
            }
        }
    }

    pub(crate) fn save_all(&mut self) {
        for (umi, pairs) in
        <HashMap<Vec<u8>, HashSet<(fastq::Record, fastq::Record)>> as Clone>::clone(&self.umi_bins).into_iter() {
            match self.collision_resolution_method {
                UMICollisionResolutionMethod::KeepFirst => {
                    // these records are already on disk
                }
                UMICollisionResolutionMethod::None | UMICollisionResolutionMethod::QualityVote => {
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
use std::{io, iter};
use std::cmp::{max, Ordering};
use std::collections::{HashMap, HashSet};
use std::io::BufWriter;

use bio::bio_types::sequence::SequenceRead;
use bio::io::fastq;
use itertools::Itertools;
use strum_macros::VariantArray;

use crate::types::{BaseQualityVotes, FastqPair, OutputWriters, QualityVoteTotal, QualityVoteVec, UMIVec, WhichRead};
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

pub(crate) struct PairHandler {
    pub(crate) record_writers: OutputWriters,
    pub(crate) collision_resolution_method: UMICollisionResolutionMethod,
    pub(crate) umi_bins: HashMap<UMIVec, HashSet<FastqPair>>,
    pub(crate) records_total: usize,
    pub(crate) records_good: usize,
    pub(crate) records_written: usize,
    pub(crate) records_unpaired: (usize, usize),
    // ATCG order, only populated if --crm quality-vote
    pub(crate) quality_votes: HashMap<UMIVec, (QualityVoteVec, QualityVoteVec)>,
}

impl Default for PairHandler {
    fn default() -> Self {
        PairHandler {
            record_writers: OutputWriters {
                paired: (
                    fastq::Writer::from_bufwriter(BufWriter::new(WriterMaybeGzip::NULL(io::sink()))),
                    fastq::Writer::from_bufwriter(BufWriter::new(WriterMaybeGzip::NULL(io::sink())))
                ),
                unpaired: (
                    fastq::Writer::from_bufwriter(BufWriter::new(WriterMaybeGzip::NULL(io::sink()))),
                    fastq::Writer::from_bufwriter(BufWriter::new(WriterMaybeGzip::NULL(io::sink())))
                ),
            },
            collision_resolution_method: UMICollisionResolutionMethod::KeepFirst,
            umi_bins: Default::default(),
            records_total: 0,
            records_good: 0,
            records_written: 0,
            records_unpaired: (0, 0),
            quality_votes: Default::default(),
        }
    }
}

impl PairHandler {
    pub(crate) fn write_pair(&mut self, pair: FastqPair) {
        // TODO: reimplement slicing etc; increment some kind of dropped record counter
        self.records_written += 1;

        // if rec_fwr.seq().len() < (start_index_fwr + 1) as usize ||
        //     rec_rev.seq().len() < (start_index_rev + 1) as usize {
        //     continue;
        // }

        self.record_writers.paired.0.write(
            std::str::from_utf8(pair.0.name()).unwrap(),
            Option::from(pair.0.id()),
            pair.0.seq(),
            pair.0.qual(),
        )
            .expect("couldn't write out a forward record");
        self.record_writers.paired.0.write(
            std::str::from_utf8(pair.1.name()).unwrap(),
            Option::from(pair.1.id()),
            pair.1.seq(),
            pair.1.qual(),
        ).expect("couldn't write out a reverse record");
    }

    pub(crate) fn write_unpaired(&mut self, record: fastq::Record, which_read: WhichRead) {
        match which_read {
            WhichRead::FORWARD => {
                self.records_unpaired.0 += 1;
                self.record_writers.unpaired.0.write(
                    std::str::from_utf8(record.name()).unwrap(),
                    Option::from(record.id()),
                    record.seq(),
                    record.qual(),
                ).expect("couldn't write out an unpaired forward record")
            }
            WhichRead::REVERSE => {
                self.records_unpaired.1 += 1;
                self.record_writers.unpaired.1.write(
                    std::str::from_utf8(record.name()).unwrap(),
                    Option::from(record.id()),
                    record.seq(),
                    record.qual(),
                ).expect("couldn't write out an unpaired reverse record")
            }
        };
    }

    pub(crate) fn insert_pair(&mut self, umi: &UMIVec, pair: &FastqPair) {
        match self.collision_resolution_method {
            // special case: no comparison, etc., just go straight to disk
            UMICollisionResolutionMethod::None => {
                self.records_good += 1;

                let umi_space = [String::from_utf8(umi.clone()).unwrap(), " ".parse().unwrap()].concat();
                // write the record, add UMI
                let id_prefix = match umi.len() {
                    0 => "",
                    _ => &umi_space
                };
                let pair_new = (
                    fastq::Record::with_attrs(
                        &*(id_prefix.to_owned() + std::str::from_utf8(pair.0.name()).unwrap()),
                        pair.0.desc(),
                        &*pair.0.seq(),
                        &*pair.0.qual(),
                    ),
                    fastq::Record::with_attrs(
                        &*(id_prefix.to_owned() + std::str::from_utf8(pair.1.name()).unwrap()),
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
                            // handled above and returned to caller, nothing needs to be done down here
                            unreachable!();
                        }
                        // special cases: `umi_bins` is involved but only to indicate a UMI has been seen
                        UMICollisionResolutionMethod::QualityVote => {
                            // create the "ballots" and save to disk later
                            let mut votes = (
                                Vec::<BaseQualityVotes>::new(), Vec::<BaseQualityVotes>::new()
                            );
                            votes.0.extend(iter::repeat((0, 0, 0, 0)).take(pair.0.len() - umi.len()));
                            votes.1.extend(iter::repeat((0, 0, 0, 0)).take(pair.1.len()));

                            Self::update_vote_vec(&mut votes, pair, umi.len());
                            self.quality_votes.insert(umi.clone(), votes);
                        }
                        UMICollisionResolutionMethod::KeepFirst => {
                            // write the record immediately; save memory
                            self.write_pair(pair.clone());
                            // save an empty set so we don't come here again
                        }
                        // un-special cases: full comparison with the contents of `umi_bins` is necessary
                        _ => {
                            // otherwise, we need to save this
                            set.insert(pair.clone());
                        }
                    }
                    // count the good record once for all these cases it could be overridden but the +1 will not change
                    self.records_good += 1;
                    // whatever we decided to put in the set (if anything), save it
                    self.umi_bins.insert(umi.clone(), set);
                } else {
                    let set = self.umi_bins.get_mut(umi).unwrap();

                    match self.collision_resolution_method {
                        UMICollisionResolutionMethod::None | UMICollisionResolutionMethod::KeepFirst => {
                            // already handled above, no need for anything involving the set
                        }
                        // need to do a bit
                        UMICollisionResolutionMethod::QualityVote => {
                            // update the "ballots"
                            let mut votes = self.quality_votes.get_mut(umi).unwrap();
                            // stretch to size sufficient to fit data
                            votes.0.extend(
                                iter::repeat((0, 0, 0, 0))
                                    .take(max((pair.0.seq().len() - umi.len()).saturating_sub(votes.0.len()), 0)));
                            votes.1.extend(
                                iter::repeat((0, 0, 0, 0))
                                    .take(max(pair.1.seq().len().saturating_sub(votes.1.len()), 0)));

                            Self::update_vote_vec(&mut votes, pair, umi.len());
                        }
                        // un-special cases, again
                        UMICollisionResolutionMethod::KeepLast => {
                            // always replace, without any checks
                            set.clear();
                            set.insert(pair.clone());
                        }
                        UMICollisionResolutionMethod::KeepLongestLeft | UMICollisionResolutionMethod::KeepLongestRight |
                        UMICollisionResolutionMethod::KeepLongestExtend => {
                            // this clone isn't the best but meh
                            let old = set.iter().next().unwrap().clone();
                            set.remove(&old);

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

    fn update_vote_vec(votes: &mut (QualityVoteVec, QualityVoteVec), pair: &FastqPair, umi_len: usize) {
        for (ind, (base, qual)) in pair.0.seq().iter()
            .zip(pair.0.qual()).dropping(umi_len)
            .enumerate() {
            match base.to_ascii_uppercase() {
                b'A' => votes.0.get_mut(ind).unwrap().0 += *qual as QualityVoteTotal,
                b'T' => votes.0.get_mut(ind).unwrap().1 += *qual as QualityVoteTotal,
                b'C' => votes.0.get_mut(ind).unwrap().2 += *qual as QualityVoteTotal,
                b'G' => votes.0.get_mut(ind).unwrap().3 += *qual as QualityVoteTotal,
                b'N' => {}  // this read abstains for this base
                _ => unimplemented!()
            }
        }

        for (ind, (base, qual)) in pair.1.seq().iter()
            .zip(pair.1.qual())
            .enumerate() {
            match base.to_ascii_uppercase() {
                b'A' => votes.1.get_mut(ind).unwrap().0 += *qual as QualityVoteTotal,
                b'T' => votes.1.get_mut(ind).unwrap().1 += *qual as QualityVoteTotal,
                b'C' => votes.1.get_mut(ind).unwrap().2 += *qual as QualityVoteTotal,
                b'G' => votes.1.get_mut(ind).unwrap().3 += *qual as QualityVoteTotal,
                b'N' => {}  // this read abstains for this base
                _ => unimplemented!()
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
                    let votes = self.quality_votes.get(&umi).unwrap();

                    // for a tuple of vote totals:
                    let count_votes = |totals: &BaseQualityVotes| -> u8 {
                        // for each possible base (0..4), fetch the number of votes for that base
                        match (0..4).max_by_key(|i| match i {
                            0 => totals.0,
                            1 => totals.1,
                            2 => totals.2,
                            3 => totals.3,
                            _ => unimplemented!()
                        }).unwrap() {  // now convert the winning index to a base
                            0 => b'A',
                            1 => b'T',
                            2 => b'C',
                            3 => b'G',
                            _ => unimplemented!()
                        }
                    };
                    let resolved: (Vec<u8>, Vec<u8>) = (
                        votes.0.iter().map(count_votes).collect(), votes.1.iter().map(count_votes).collect());

                    // this quality score is entirely fake
                    self.write_pair((
                        fastq::Record::with_attrs(
                            std::str::from_utf8(&*umi).unwrap(),
                            Option::from("constructed by grebe from quality voting"),
                            &*resolved.0,
                            &*iter::repeat(b"~").take(resolved.0.len()).map(|x| x[0]).collect::<Vec<u8>>(),
                        ),
                        fastq::Record::with_attrs(
                            std::str::from_utf8(&*umi).unwrap(),
                            Option::from("constructed by grebe from quality voting"),
                            &*resolved.1,
                            &*iter::repeat(b"~").take(resolved.1.len()).map(|x| x[0]).collect::<Vec<u8>>(),
                        )
                    ));
                }
                _ => {
                    // conflict resolution has already selected a single read
                    self.write_pair(pairs.iter().next().unwrap().clone());
                }
            };
        }
    }
}
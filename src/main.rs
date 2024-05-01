use std::cmp::min;
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::process::exit;

use clap::{ArgGroup, ValueEnum, ValueHint};
use clap::builder::PossibleValue;
use clap::parser::ValueSource;
use editdistancek::edit_distance_bounded;
use indicatif::{ProgressBar, ProgressFinish};
use itertools::Itertools;
use pluralizer::pluralize;
use strum::VariantArray;

use pair_handler::UMICollisionResolutionMethod;
use types::FastqPair;

use crate::pair_handler::PairHandler;
use crate::reader::make_reader_pair;
use crate::types::UMIVec;

mod pair_handler;
mod reader;
mod writer;
mod types;

fn find_within_radius(umi_bins: &HashMap<UMIVec, HashSet<FastqPair>>, umi: &UMIVec, radius: usize)
                      -> Option<UMIVec> {
    umi_bins.keys().find(|proposed_umi| edit_distance_bounded(proposed_umi, umi, radius).is_some()).cloned()
}

impl ValueEnum for UMICollisionResolutionMethod {
    fn value_variants<'a>() -> &'a [Self] { &Self::VARIANTS }

    fn to_possible_value(&self) -> Option<PossibleValue> {
        Some(match self {
            Self::None => PossibleValue::new("none")
                .help("keep duplicates, prepend their assigned UMI to their sequence names"),
            UMICollisionResolutionMethod::KeepFirst => PossibleValue::new("first")
                .alias("keep-first")
                .alias("keep-left")
                .alias("kl")
                .help("keep the first sequence matched to this UMI, ignore any sequences that follow"),
            UMICollisionResolutionMethod::KeepLast => PossibleValue::new("last")
                .alias("keep-last")
                .alias("keep-right")
                .alias("kr")
                .help("keep the last sequence matched, ignore any sequences that come before"),
            UMICollisionResolutionMethod::KeepLongestLeft => PossibleValue::new("keep-longest-left")
                .alias("kl-left")
                .alias("kll")
                .help("keep the longest sequence matched, favor the earlier sequence when tied"),
            UMICollisionResolutionMethod::KeepLongestRight => PossibleValue::new("keep-longest-right")
                .alias("kl-right")
                .alias("klr")
                .help("keep the longest sequence matched, favor the later sequence when tied"),
            UMICollisionResolutionMethod::KeepLongestExtend => PossibleValue::new("keep-longest-extend")
                .alias("extend")
                .alias("kl-extend")
                .alias("kle")
                .help("keep the longest sequence, overwrite it if a read found later is longer and agrees completely \
                on base calls"),
            UMICollisionResolutionMethod::QualityVote => PossibleValue::new("quality-vote")
                .alias("quality-voting")
                .alias("vote")
                .alias("voting")
                .alias("qv")
                .help("create one final sequence by combining base calls and qualities from all matched reads"),
        })
    }
}


fn main() {
    let cmd = clap::command!("grebe")
        .about("Processing tool for Illumina sequencing data")
        .arg(clap::arg!(<"in-forward"> "forward (5'-3') reads to work with")
            .value_name("input forward .fastq")
            .value_parser(clap::value_parser!(PathBuf))
            .value_hint(ValueHint::FilePath))
        .arg(clap::arg!(<"in-reverse"> "reverse (3'-5') reads to work with")
            .value_name("input reverse .fastq")
            .value_parser(clap::value_parser!(PathBuf))
            .value_hint(ValueHint::FilePath))
        .arg(clap::arg!(--"phred64" "use the legacy phred64 encoding (over phred33) where score 0 \
        = \"@\" instead of \"!\"")
            .required(false)
            .default_value("false"))
        .arg(clap::arg!(-'u' <"UMI length"> "UMI length (strip this many bases off forward reads)")
            .id("umi-length")
            .visible_alias("umi-length")
            .value_parser(0..=15)
            .required(false)
            .default_value("0"))
        .arg(clap::arg!(--"collision-resolution-mode" <"mode"> "choose how to resolve UMI collisions")
            .alias("collision-resolution-method")
            .alias("conflict-resolution-mode")
            .alias("conflict-resolution-method")
            .visible_alias("crm")
            .value_parser(clap::value_parser!(UMICollisionResolutionMethod))
            .default_value("keep-first"))
        .arg(clap::arg!(-'l' <"levenshtein radius"> "bin UMIs together if at most this Levenshtein distance apart \
        (good for small library error tolerance, but very slow on genomic-scale data)")
            .id("levenshtein-radius")
            .visible_alias("levenshtein")
            .visible_alias("levenshtein-radius")
            .visible_alias("lr")
            .value_parser(0..=15)
            .required(false)
            .default_value("0"))
        .arg(clap::arg!(--"proactive-levenshtein" <"force-mode"> "(for advanced users) force guess-and-check approach \
        to levenshtein distance (default is true if -l is 1 or 2, false otherwise)")
            .visible_alias("pl")
            .value_parser(clap::value_parser!(bool))
            .required(false))
        .arg(clap::arg!(--"start-at" <"start index"> "start reads after this many base pairs (but process UMIs even if \
        they would be clipped); reads which become empty are dropped")
            .visible_alias("--start-index")
            .value_parser(0..=600)
            .required(false)
            .default_value("0"))
        .group(ArgGroup::new("left-slice")
            .arg("start-at"))
        .arg(clap::arg!(<"out-forward"> "where to place processed forward reads")  // TODO: output more sequence formats
            .value_name("output forward .fastq")
            .value_parser(clap::value_parser!(PathBuf))
            .value_hint(ValueHint::FilePath))
        .arg(clap::arg!(<"out-reverse"> "where to place processed reverse reads")
            .value_name("output reverse .fastq")
            .value_parser(clap::value_parser!(PathBuf))
            .value_hint(ValueHint::FilePath));
    // .arg(clap::arg!(<adapters>)
    //     .value_name("illumina adapters .fasta")
    //     .help("illumina adapters file (to prune them)")
    //     .value_parser(clap::value_parser!(std::path::PathBuf)));

    let args = cmd.get_matches();

    let input_paths = (
        args.get_one::<PathBuf>("in-forward").unwrap(),
        args.get_one::<PathBuf>("in-reverse").unwrap()
    );
    let record_readers = make_reader_pair(input_paths, true);
    // TODO: reverse reads here too
    let total_records = record_readers.0.records().count();
    let record_readers = make_reader_pair(input_paths, false);

    let record_writers = writer::make_writer_pair((
        args.get_one::<PathBuf>("out-forward").unwrap(),
        args.get_one::<PathBuf>("out-reverse").unwrap()
    ));

    let umi_length = *args.get_one::<i64>("umi-length").unwrap();

    let collision_resolution_method = args.get_one::<UMICollisionResolutionMethod>("collision-resolution-mode")
        .unwrap().to_owned();
    let mut pair_handler = PairHandler {
        record_writers,
        collision_resolution_method,
        total_records,
        ..Default::default()
    };

    // let start_index_arg = *args.get_one::<i64>("start-at").unwrap();
    // let start_index_rev = start_index_arg;
    // let start_index_fwr = max(start_index_arg, umi_length);

    let levenshtein_max = min(*args.get_one::<i64>("levenshtein-radius").unwrap(), umi_length);
    if levenshtein_max >= umi_length && args.value_source("levenshtein-radius") == Some(ValueSource::CommandLine) {
        eprintln!("warning: --levenshtein-max too high to be meaningful")
    }

    // TODO: debug print here
    let proactive_levenshtein = match args.get_one::<bool>("proactive-levenshtein") {
        Some(result) => {
            if umi_length == 0 {
                eprintln!("warning: --proactive-levenshtein is meaningless with no UMI")
            } else if levenshtein_max == 0 {
                eprintln!("warning: --proactive-levenshtein is meaningless with -l 0")
            }
            *result
        }
        None => levenshtein_max <= 2
    };

    let bar = ProgressBar::new(total_records as u64).with_finish(ProgressFinish::AndLeave);

    let pairs = record_readers.0.records().zip(record_readers.1.records());
    'pairs: for (rec_fwr, rec_rev) in pairs {
        bar.inc(1);

        let rec_fwr = match rec_fwr {
            Ok(result) => match result.check() {
                Ok(_) => result,
                Err(err) => {
                    eprintln!("forward record {} was invalid: {err}", pair_handler.total_records);
                    exit(1);
                }
            },
            Err(_) => {
                eprintln!("forward record {} was invalid", pair_handler.total_records);
                exit(1);
            }
        };
        let rec_rev = match rec_rev {
            Ok(result) => match result.check() {
                Ok(_) => result,
                Err(err) => {
                    eprintln!("reverse record {} was invalid: {err}", pair_handler.total_records);
                    exit(1);
                }
            },
            Err(_) => {
                eprintln!("reverse record {} was invalid", pair_handler.total_records);
                exit(1);
            }
        };

        if umi_length > 0 {
            let umi: UMIVec = rec_fwr.seq()[..umi_length as usize].iter().copied().collect();
            if levenshtein_max == 0 {
                pair_handler.insert_pair(&umi, &(rec_fwr, rec_rev));
            } else {
                if proactive_levenshtein {
                    // instead of checking the distance to elements of the set of known UMIs,
                    // generate UMIs within a certain distance and check them
                    // TODO: assumes no Ns outside of masked reads

                    // first, generate all options for <levenshtein_max> new base values
                    let new_bases = std::iter::repeat("ATCG".chars())
                        .take(levenshtein_max as usize)
                        .multi_cartesian_product();
                    // then, generate all options for <levenshtein_max> positions to replace at
                    for indices_to_replace in (0..umi_length).combinations(levenshtein_max as usize) {
                        // execute the replacement
                        for base_substitution in new_bases.clone() {
                            let mut umi_modified = umi.clone();
                            for (index, new_value) in (&indices_to_replace).iter().zip(base_substitution) {
                                umi_modified[*index as usize] = new_value as u8;
                            }
                            // if our result is already known, bail out
                            if pair_handler.umi_bins.contains_key(&umi_modified) {
                                pair_handler.insert_pair(&umi_modified, &(rec_fwr, rec_rev));
                                continue 'pairs;
                            }
                        }
                    }

                    // no proposed alternative was satisfactory; we have a new UMI
                    pair_handler.insert_pair(&umi, &(rec_fwr, rec_rev));
                } else {
                    // non-proactive mode; just check every known UMI and see if it's close enough

                    // check for perfect match first
                    if pair_handler.umi_bins.contains_key(&umi) {
                        pair_handler.insert_pair(&umi, &(rec_fwr, rec_rev));
                    } else {
                        match find_within_radius(&pair_handler.umi_bins, &umi, levenshtein_max as usize) {
                            None => pair_handler.insert_pair(&umi, &(rec_fwr, rec_rev)),
                            Some(found) => pair_handler.insert_pair(&found, &(rec_fwr, rec_rev))
                        }
                    }
                }
            }
        } else {
            pair_handler.collision_resolution_method = UMICollisionResolutionMethod::None;
            pair_handler.insert_pair(&vec![], &(rec_fwr, rec_rev));
        }
    }

    bar.finish_using_style();

    if umi_length > 0 {
        if pair_handler.records_written > 0 {
            println!("filtered {} down to {} via UMI, wrote {} after pair-level filtering",
                     pluralize("pair", pair_handler.total_records as isize, true),
                     pluralize("pair", pair_handler.good_records as isize, true),
                     pluralize("pair", pair_handler.records_written as isize, true))
        } else {
            println!("filtered {} down to {} via UMI; writing to disk...",
                     pluralize("pair", pair_handler.total_records as isize, true),
                     pluralize("pair", pair_handler.good_records as isize, true));
        }
    } else {
        println!("processed {}, wrote {} after pair-level filtering",
                 pluralize("pair", pair_handler.total_records as isize, true),
                 pluralize("pair", pair_handler.records_written as isize, true))
    }

    pair_handler.save_remaining();
    // TODO: verbose logging (masked reads, etc.)
    // TODO: exit codes
    // TODO: do things on quality scores
}

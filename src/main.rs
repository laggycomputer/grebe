use std::cmp::{max, min};
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::process::exit;

use bio::io::fastq;
use clap::{ArgGroup, ValueHint};
use clap::parser::ValueSource;
use editdistancek::edit_distance_bounded;
use itertools::Itertools;
use pluralizer::pluralize;

use crate::pair_handler::{PairHandler, UMICollisionResolutionMethod};
use crate::reader::reader_maybe_gzip;
use crate::writer::writer_maybe_gzip;

mod pair_handler;
mod reader;
mod writer;

type FastqPair = (fastq::Record, fastq::Record);

fn find_within_radius(umi_bins: &HashMap<Vec<u8>, HashSet<FastqPair>>, umi: &Vec<u8>, radius: usize)
                      -> Option<Vec<u8>> {
    umi_bins.keys().find(|proposed_umi| edit_distance_bounded(proposed_umi, umi, radius).is_some()).cloned()
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
        = \"@\" instead of \"1\"")
            .required(false)
            .default_value("false"))
        .arg(clap::arg!(-'u' <"UMI length"> "UMI length (strip this many bases off forward reads)")
            .id("umi-length")
            .visible_alias("umi-length")
            .value_parser(0..=15)
            .required(false)
            .default_value("0"))
        .arg(clap::arg!(--"collision-resolution-mode" <"mode"> "choose how to resolve UMI collisions")
            .visible_alias("collision-resolution-method")
            .visible_alias("conflict-resolution-mode")
            .visible_alias("conflict-resolution-method")
            .visible_alias("crm")
            .value_parser(clap::value_parser!(pair_handler::UMICollisionResolutionMethod))
            .default_value("keep-first"))
        .arg(clap::arg!(-'l' <"levenshtein radius"> "(0 to disable) bin UMIs together if at most this Levenshtein \
        distance apart (useful for small libraries to reduce error rates, but very slow on genomic-scale data)")
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
    let record_readers = (
        match reader_maybe_gzip(input_paths.0) {
            Ok((result, was_compressed)) => {
                if was_compressed { eprintln!("info: parsing {} as a gzip", input_paths.0.display()) }
                result
            }
            Err(_) => {
                eprintln!("couldn't open input forward .fastq");
                exit(1);
            }
        },
        match reader_maybe_gzip(input_paths.1) {
            Ok((result, was_compressed)) => {
                if was_compressed { eprintln!("info: parsing {} as a gzip", input_paths.1.display()) }
                result
            }
            Err(_) => {
                eprintln!("couldn't open input reverse .fastq");
                exit(1);
            }
        }
    );

    let output_paths = (
        args.get_one::<PathBuf>("out-forward").unwrap(),
        args.get_one::<PathBuf>("out-reverse").unwrap()
    );
    let record_writers = (
        match writer_maybe_gzip(output_paths.0) {
            Ok((result, was_compressed)) => {
                if was_compressed { eprintln!("info: writing {} as a gzip", output_paths.0.display()) }
                result
            }
            Err(_) => {
                eprintln!("couldn't open output forward .fastq; if it exists, remove it");
                exit(1);
            }
        },
        match writer_maybe_gzip(output_paths.1) {
            Ok((result, was_compressed)) => {
                if was_compressed { eprintln!("info: writing {} as a gzip", output_paths.1.display()) }
                result
            }
            Err(_) => {
                eprintln!("couldn't open output reverse .fastq; if it exists, remove it");
                exit(1);
            }
        }
    );

    let umi_length = *args.get_one::<i64>("umi-length").unwrap();

    let collision_resolution_method = args.get_one::<UMICollisionResolutionMethod>("collision-resolution-mode")
        .unwrap().to_owned();
    let mut pair_handler = PairHandler {
        record_writers,
        collision_resolution_method,
        umi_bins: HashMap::new(),
    };

    let start_index_arg = *args.get_one::<i64>("start-at").unwrap();
    let start_index_rev = start_index_arg;
    let start_index_fwr = max(start_index_arg, umi_length);

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

    // TODO: move all state to pair_handler where possible, will fix good_record count
    let mut total_records = 0;
    let mut good_records = 0;
    let pairs = record_readers.0.records().zip(record_readers.1.records());
    'pairs: for (rec_fwr, rec_rev) in pairs {
        total_records += 1;

        let rec_fwr = match rec_fwr {
            Ok(result) => result,
            Err(_) => {
                eprintln!("forward record {total_records} was invalid");
                exit(1);
            }
        };
        let rec_rev = match rec_rev {
            Ok(result) => result,
            Err(_) => {
                eprintln!("reverse record {total_records} was invalid");
                exit(1);
            }
        };

        if rec_fwr.seq().len() < (start_index_fwr + 1) as usize ||
            rec_rev.seq().len() < (start_index_rev + 1) as usize {
            continue;
        }

        let all_ns = |s: &u8| *s == ('N' as u8);
        if rec_fwr.seq().iter().all(all_ns) || rec_rev.seq().iter().all(all_ns) {
            continue;
        }

        if umi_length > 0 {
            let umi: Vec<u8> = rec_fwr.seq()[..umi_length as usize].iter().copied().collect();
            if levenshtein_max == 0 {
                pair_handler.handle_pair(&umi, &(rec_fwr, rec_rev));
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
                                pair_handler.handle_pair(&umi_modified, &(rec_fwr, rec_rev));
                                continue 'pairs;
                            }
                        }
                    }

                    // no proposed alternative was satisfactory; we have a new UMI
                    pair_handler.handle_pair(&umi, &(rec_fwr, rec_rev));
                } else {
                    if pair_handler.umi_bins.contains_key(&umi) {
                        pair_handler.handle_pair(&umi, &(rec_fwr, rec_rev));
                    } else {
                        match find_within_radius(&pair_handler.umi_bins, &umi, levenshtein_max as usize) {
                            None => pair_handler.handle_pair(&umi, &(rec_fwr, rec_rev)),
                            Some(found) => pair_handler.handle_pair(&found, &(rec_fwr, rec_rev))
                        }
                    }
                }
            }
        }

        good_records += 1;
    }

    println!("filtered {} down to {}; writing to disk...", pluralize("pair", total_records, true),
             pluralize("pair", good_records, true));

    pair_handler.save_all();
    // TODO: count records before starting and give progress indication
    // TODO: stop assuming forward and reverse reads appear in the proper order
    // TODO: verbose logging (masked reads, etc.)
    // TODO: exit codes
    // TODO: do things on quality scores
    // TODO: quality vote resolution
}

use std::cmp::{max, min};
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::process::exit;

use bio::alignment::distance::simd::hamming;
use bio::alphabets::dna;
use clap::{ArgGroup, ValueEnum, ValueHint};
use clap::builder::PossibleValue;
use clap::parser::ValueSource;
use indicatif::{ProgressBar, ProgressFinish};
use itertools::Itertools;
use pluralizer::pluralize;
use strum::VariantArray;

use pair_handling::UMICollisionResolutionMethod;
use types::FastqPair;

use crate::pair_handling::PairHandler;
use crate::reader::make_reader_pair;
use crate::types::{OutputWriters, UMIVec, WhichRead};
use crate::util::check_primer;

mod pair_handling;
mod reader;
mod writer;
mod types;
mod util;

fn find_within_radius(umi_bins: &HashMap<UMIVec, HashSet<FastqPair>>, umi: &UMIVec, radius: usize)
                      -> Option<UMIVec> {
    umi_bins.keys().find(|proposed_umi| hamming(proposed_umi, umi) <= radius as u64).cloned()
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
        .arg(clap::arg!(--"hr" <"hamming radius"> "bin UMIs together if at most this Hamming distance apart \
        (good for small library error tolerance, but very slow on genomic-scale data)")
            .id("hamming-radius")
            .visible_alias("hamming")
            .visible_alias("hamming-radius")
            .value_parser(0..=15)
            .required(false)
            .default_value("0"))
        .arg(clap::arg!(--"proactive-binning" <"force mode"> "(for advanced users, see docs; you shouldn't have to \
        set this)")
            .visible_alias("pb")
            .value_parser(clap::value_parser!(bool))
            .required(false))
        .arg(clap::arg!(--"forward-primer" <"forward primer"> "(IUPAC alphabet allowed) ensure forward reads begin \
        with this sequence; if -u is specified, forwards starting with this primer are considered failed UMI additions \
        and pair is discarded")
            .visible_alias("fp")
            .visible_alias("primer-forward")
            .visible_alias("pf")
            .required(false))
        .arg(clap::arg!(--"reverse-primer" <"reverse primer"> "(IUPAC alphabet allowed) ensure reverse reads end with \
        the reverse complement of this sequence; non-conforming pairs are discarded")
            .visible_alias("rp")
            .visible_alias("primer-reverse")
            .visible_alias("pr")
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
            .value_hint(ValueHint::FilePath))
        .arg(clap::arg!(["out-unpaired-forward"] "where to place unpaired forward reads")
            .value_name("output unpaired forward .fastq")
            .value_parser(clap::value_parser!(PathBuf))
            .value_hint(ValueHint::FilePath)
            .required(false))
        .arg(clap::arg!(["out-unpaired-reverse"] "where to place unpaired reverse reads")
            .value_name("output unpaired reverse .fastq")
            .value_parser(clap::value_parser!(PathBuf))
            .value_hint(ValueHint::FilePath)
            .required(false));
    // .arg(clap::arg!(<adapters>)
    //     .value_name("illumina adapters .fasta")
    //     .help("illumina adapters file (to prune them)")
    //     .value_parser(clap::value_parser!(std::path::PathBuf)));

    let args = cmd.get_matches();

    let umi_length = *args.get_one::<i64>("umi-length").unwrap() as u8;

    let collision_resolution_method;
    if umi_length == 0 {
        // silently override this; --crm is meaningless in this context
        collision_resolution_method = UMICollisionResolutionMethod::None;
    } else {
        collision_resolution_method = args.get_one::<UMICollisionResolutionMethod>("collision-resolution-mode")
            .unwrap().to_owned();
    }

    // let start_index_arg = *args.get_one::<i64>("start-at").unwrap();
    // let start_index_rev = start_index_arg;
    // let start_index_fwr = max(start_index_arg, umi_length);

    let hamming_radius = min(*args.get_one::<i64>("hamming-radius").unwrap() as u8, umi_length);
    if hamming_radius >= umi_length && args.value_source("hamming-radius") == Some(ValueSource::CommandLine) {
        eprintln!("warning: --hamming-max too high to be meaningful")
    }

    // TODO: debug print here
    let proactive_binning = match args.get_one::<bool>("proactive-hamming") {
        Some(result) => {
            if umi_length == 0 {
                eprintln!("warning: --proactive_binning is meaningless with no UMI")
            } else if hamming_radius == 0 {
                eprintln!("warning: --proactive_binning is meaningless with -l 0")
            }
            *result
        }
        // --pl true's intelligent binning makes it much slower for --crm none
        None => hamming_radius <= 3 && collision_resolution_method != UMICollisionResolutionMethod::None
    };


    let check_iupac_dna = |p: &String| !dna::iupac_alphabet().is_word(p.as_bytes());
    if args.get_one::<String>("forward-primer").is_some_and(check_iupac_dna) {
        eprintln!("forward primer not valid IUPAC DNA alphabet; refusing");
        exit(1);
    }
    if args.get_one::<String>("reverse-primer").is_some_and(check_iupac_dna) {
        eprintln!("reverse primer not valid IUPAC DNA alphabet; refusing");
        exit(1);
    }

    let enforce_primers = (
        args.get_one::<String>("forward-primer").map(|s| s.as_bytes()),
        args.get_one::<String>("reverse-primer").map(|s| s.as_bytes())
    );

    let input_paths = (
        args.get_one::<PathBuf>("in-forward"),
        args.get_one::<PathBuf>("in-reverse")
    );
    let record_readers = make_reader_pair(input_paths, true);
    let total_records = (record_readers.0.records().count(), record_readers.1.records().count());

    let record_readers = make_reader_pair(input_paths, false);

    let record_writers = OutputWriters {
        paired: writer::make_writer_pair((
            args.get_one::<PathBuf>("out-forward"),
            args.get_one::<PathBuf>("out-reverse")
        )),
        unpaired: writer::make_writer_pair((
            args.get_one::<PathBuf>("out-unpaired-forward"),
            args.get_one::<PathBuf>("out-unpaired-reverse")
        )),
    };

    let mut pair_handler = PairHandler {
        record_writers,
        collision_resolution_method,
        records_total: max(total_records.0, total_records.1),
        ..Default::default()
    };

    eprintln!("counted {}, working...", pluralize("pair", pair_handler.records_total as isize, true));
    let bar = ProgressBar::new(pair_handler.records_total as u64).with_finish(ProgressFinish::AndLeave);

    let pairs = record_readers.0.records().zip(record_readers.1.records());
    'pairs: for maybe_read_pair in pairs {
        bar.inc(1);

        // basic layout for this code:
        // 1. verify validity of reads
        // 2. verify validity of pair overall (paired/unpaired, sufficient length, primer, etc.)
        // 3. match UMI and allow handler struct to decide what to do from there

        // these checks permit me to go insane and unsafe every string parse
        let read_pair = (
            match maybe_read_pair.0 {
                Ok(result) => match result.check() {
                    Ok(_) => result,
                    Err(err) => {
                        eprintln!("forward record {} was invalid: {err}", pair_handler.records_total);
                        exit(1);
                    }
                },
                Err(_) => {
                    eprintln!("forward record {} was invalid", pair_handler.records_total);
                    exit(1);
                }
            },
            match maybe_read_pair.1 {
                Ok(result) => match result.check() {
                    Ok(_) => result,
                    Err(err) => {
                        eprintln!("reverse record {} was invalid: {err}", pair_handler.records_total);
                        exit(1);
                    }
                },
                Err(_) => {
                    eprintln!("reverse record {} was invalid", pair_handler.records_total);
                    exit(1);
                }
            }
        );

        let n_closure = |s: &u8| *s == b'N';
        match (read_pair.0.seq().iter().all(n_closure), read_pair.1.seq().iter().all(n_closure)) {
            (true, false) => {
                pair_handler.write_unpaired(read_pair.0, WhichRead::REVERSE);
                continue 'pairs;
            }
            (false, true) => {
                pair_handler.write_unpaired(read_pair.1, WhichRead::FORWARD);
                continue 'pairs;
            }
            (true, true) => {
                pair_handler.pair_drop_reason_count.both_masked += 1;
                continue 'pairs;
            }
            _ => {}
        }

        if enforce_primers.0.is_some() {
            if read_pair.0.seq().len() < umi_length as usize + enforce_primers.0.unwrap().len() {
                pair_handler.pair_drop_reason_count.no_forward_primer += 1;
                continue 'pairs;
            }

            let starts_with_primer = check_primer(enforce_primers.0.as_ref().unwrap(), &read_pair.0.seq())
                .unwrap_or_default();
            let starts_with_umi_then_primer = check_primer(
                enforce_primers.0.as_ref().unwrap(),
                &read_pair.0.seq()[umi_length as usize..],
            ).unwrap_or_default();

            if umi_length > 0 && starts_with_primer && !starts_with_umi_then_primer {
                // very unlikely the UMI then following seq is the primer; we will call this a bad UMI addition
                pair_handler.pair_drop_reason_count.umi_is_forward_primer += 1;
                continue 'pairs;
            } else if !starts_with_umi_then_primer {
                // primer not present where it should be
                pair_handler.pair_drop_reason_count.no_forward_primer += 1;
                continue 'pairs;
            }
        }
        if enforce_primers.1.is_some() {
            if read_pair.1.seq().len() < enforce_primers.1.unwrap().len() {
                pair_handler.pair_drop_reason_count.no_reverse_primer += 1;
                continue 'pairs;
            }

            let starts_with_primer = check_primer(enforce_primers.1.as_ref().unwrap(), read_pair.1.seq())
                .unwrap_or_default();

            if !starts_with_primer {
                pair_handler.pair_drop_reason_count.no_reverse_primer += 1;
                continue 'pairs;
            }
        }

        if umi_length > 0 {
            let umi: UMIVec = read_pair.0.seq()[..umi_length as usize].iter().copied().collect();
            if hamming_radius == 0 {
                pair_handler.insert_pair(&umi, &read_pair);
            } else {
                if pair_handler.umi_bins.contains_key(&umi) {
                    pair_handler.insert_pair(&umi, &read_pair);
                    continue 'pairs;
                }

                if proactive_binning {
                    // instead of checking the distance to elements of the set of known UMIs,
                    // generate UMIs within a certain distance and check them
                    // TODO: assumes no Ns outside of masked reads

                    let mut found_bins = HashSet::new();

                    // first, generate all options for <hamming_max> new base values
                    let new_bases = std::iter::repeat("ATCG".chars())
                        .take(hamming_radius as usize)
                        .multi_cartesian_product();
                    // then, generate all options for <hamming_radius> positions to replace at
                    for indices_to_replace in (0..umi_length).combinations(hamming_radius as usize) {
                        // execute the replacement
                        for base_substitution in new_bases.clone() {
                            let mut umi_modified = umi.clone();
                            for (index, new_value) in (&indices_to_replace).iter().zip(base_substitution) {
                                umi_modified[*index as usize] = new_value as u8;
                            }

                            if collision_resolution_method == UMICollisionResolutionMethod::None {
                                found_bins.insert(umi_modified);
                            } else {
                                if pair_handler.umi_bins.contains_key(&umi_modified) {
                                    pair_handler.insert_pair(&umi_modified, &read_pair);
                                    continue 'pairs;
                                }
                            }
                        }
                    }

                    // we reach this if:
                    // no known UMI was suitable
                    // or if --crm none and at least one known UMI was suitable
                    match found_bins.into_iter()
                        .max_by_key(|k| match pair_handler.umi_bins.get(k) {
                            None => 0,
                            Some(s) => s.len()
                        }) {
                        // if the first case above is true, the iterator stops immediately and we accept a new UMI
                        None => pair_handler.insert_pair(&umi, &read_pair),
                        // if the second case is true, we have found a "best" UMI (defined as the UMI with the biggest
                        // bin) and we use that one
                        Some(umi_modified) => pair_handler.insert_pair(&umi_modified, &read_pair)
                    }
                } else {
                    // non-proactive mode; just check every known UMI and see if it's close enough

                    match find_within_radius(&pair_handler.umi_bins, &umi, hamming_radius as usize) {
                        None => pair_handler.insert_pair(&umi, &read_pair),
                        Some(found) => pair_handler.insert_pair(&found, &read_pair)
                    }
                }
            }
        } else {
            pair_handler.insert_pair(&vec![], &read_pair);
        }
    }

    bar.finish_using_style();

    let saved_verb = "wrote";
    let dropped_verb = "no save path specified; dropped";
    let verbs = (
        match args.get_one::<PathBuf>("out-unpaired-forward").is_some() {
            true => saved_verb,
            false => dropped_verb,
        }, match args.get_one::<PathBuf>("out-unpaired-reverse").is_some() {
            true => saved_verb,
            false => dropped_verb,
        }
    );
    if pair_handler.records_unpaired.0 > 0 {
        println!("{} {}", verbs.0,
                 pluralize("unpaired forward read", pair_handler.records_unpaired.0 as isize, true));
    }
    if pair_handler.records_unpaired.1 > 0 {
        println!("{} {}", verbs.1,
                 pluralize("unpaired reverse read", pair_handler.records_unpaired.1 as isize, true));
    }

    let total_dropped = pair_handler.pair_drop_reason_count.total();
    println!("dropped {} for the following reasons:\n{}",
             pluralize("pair", total_dropped as isize, true),
             pair_handler.pair_drop_reason_count);

    if umi_length > 0 {
        if pair_handler.records_written > 0 {
            // assumption: records_written = records_good
            // this is valid in the current design where all good pairs are already written to disk at this point
            // said pairs have already been "forgotten", so there is no further information to report
            println!("filtered {} to {} via UMI, wrote to disk",
                     pluralize("remaining pair", (pair_handler.records_total - total_dropped) as isize, true),
                     pluralize("pair", pair_handler.records_written as isize, true));
        } else {
            // records_good is not populated, but number of bins is
            println!("filtered {} down to {} via UMI; writing to disk...",
                     pluralize("remaining pair", (pair_handler.records_total - total_dropped) as isize, true),
                     pluralize("pair", pair_handler.umi_bins.len() as isize, true));
        }
    } else {
        // everything that survived pair filtering went straight to disk (no UMI matching)
        println!("wrote {} after pair-level filtering",
                 pluralize("remaining pair", pair_handler.records_written as isize, true));
    }

    pair_handler.write_remaining();

    // TODO: verbose logging (masked reads, etc.)
    // TODO: exit codes
    // TODO: do things on quality scores
}

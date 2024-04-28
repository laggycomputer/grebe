use std::cmp::{max, min};
use std::collections::HashSet;
use std::fs::File;
use std::io;
use std::io::{BufReader, Read};
use std::path::PathBuf;
use std::process::exit;

use bio::bio_types::sequence::SequenceRead;
use bio::io::fastq;
use clap::{ArgGroup, ValueHint};
use clap::parser::ValueSource;
use editdistancek::edit_distance_bounded;
use flate2::bufread::MultiGzDecoder;
use itertools::Itertools;
use pluralizer::pluralize;

enum ReaderMaybeGzip {
    GZIP(MultiGzDecoder<BufReader<File>>),
    UNCOMPRESSED(BufReader<File>),
}

impl Read for ReaderMaybeGzip {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            ReaderMaybeGzip::GZIP(backer) => { backer.read(buf) }
            ReaderMaybeGzip::UNCOMPRESSED(backer) => { backer.read(buf) }
        }
    }
}

fn reader_maybe_gzip(path_buf: &PathBuf) -> Result<(fastq::Reader<BufReader<ReaderMaybeGzip>>, bool), io::Error> {
    let mut file = File::open(path_buf)?;
    let mut magic = [0; 2];
    file.read(&mut magic[..])?;

    let reopen = BufReader::new(File::open(path_buf)?);

    if magic.eq(&[0x1f, 0x8b]) {
        Ok((fastq::Reader::from_bufread(BufReader::new(ReaderMaybeGzip::GZIP(MultiGzDecoder::new(reopen)))), true))
    } else {
        Ok((fastq::Reader::from_bufread(BufReader::new(ReaderMaybeGzip::UNCOMPRESSED(reopen))), true))
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
        = \"@\" instead of \"1\"")
            .required(false)
            .default_value("false"))
        .arg(clap::arg!(-'u' <"UMI length"> "UMI length (strip this many bases off forward reads)")
            .id("umi-length")
            .alias("umi-length")
            .value_parser(0..=15)
            .required(false)
            .default_value("0"))
        .arg(clap::arg!(-'l' <"levenshtein radius"> "(0 to disable) bin UMIs together if at most this Levenshtein \
        distance apart (useful for small libraries to reduce error rates, but incurs time penalty quadratic in number \
        of reads)")
            .id("levenshtein-radius")
            .alias("levenshtein")
            .alias("levenshtein-radius")
            .alias("lr")
            .value_parser(0..=15)
            .required(false)
            .default_value("0"))
        .arg(clap::arg!(--"proactive-levenshtein" <"force-mode"> "(for advanced users) force guess-and-check approach \
        to levenshtein distance (default is true if -l is 1 or 2, false otherwise)")
            .alias("pl")
            .value_parser(clap::value_parser!(bool))
            .required(false))
        .arg(clap::arg!(--"start-at" <"start index"> "start reads after this many base pairs, including any UMI \
        stripping; reads of insufficient length are dropped")
            .alias("--start-index")
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

    let in_forward = args.get_one::<PathBuf>("in-forward").unwrap();
    let reader_fwr = match reader_maybe_gzip(in_forward) {
        Ok((result, was_compressed)) => {
            if was_compressed { eprintln!("info: parsing {} as a gzip", in_forward.display()) }
            result
        }
        Err(_) => {
            eprintln!("couldn't open input forward .fastq");
            exit(1);
        }
    };
    let in_reverse = args.get_one::<PathBuf>("in-reverse").unwrap();
    let reader_rev = match reader_maybe_gzip(in_reverse) {
        Ok((result, was_compressed)) => {
            if was_compressed { eprintln!("info: parsing {} as a gzip", in_reverse.display()) }
            result
        }
        Err(_) => {
            eprintln!("couldn't open input reverse .fastq");
            exit(1);
        }
    };

    let mut writer_fwr = match fastq::Writer::to_file(args.get_one::<PathBuf>("out-forward").unwrap()) {
        Ok(result) => result,
        Err(_) => {
            eprintln!("couldn't open output forward .fastq");
            exit(1);
        }
    };
    let mut writer_rev = match fastq::Writer::to_file(args.get_one::<PathBuf>("out-reverse").unwrap()) {
        Ok(result) => result,
        Err(_) => {
            eprintln!("couldn't open output reverse .fastq");
            exit(1);
        }
    };

    let umi_length = *args.get_one::<i64>("umi-length").unwrap();
    let mut seen_umis = HashSet::new();

    let start_index_arg = *args.get_one::<i64>("start-at").unwrap();
    let start_index_rev = start_index_arg;
    let start_index_fwr = max(start_index_rev, umi_length);

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

    let mut total_records = 0;
    let mut good_records = 0;
    let pairs = reader_fwr.records().zip(reader_rev.records());
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

        if rec_fwr.seq().len() < (start_index_fwr + 1) as usize || rec_rev.seq().len() < (start_index_rev + 1) as usize {
            continue;
        }

        let all_ns = |s: &u8| *s == ('N' as u8);
        if rec_fwr.seq().iter().all(all_ns) || rec_rev.seq().iter().all(all_ns) {
            continue;
        }

        if umi_length > 0 {
            let umi: Vec<u8> = rec_fwr.seq()[..umi_length as usize].iter().copied().collect();
            if levenshtein_max == 0 {
                if !seen_umis.insert(umi) {
                    continue;
                }
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
                            if seen_umis.contains(&umi_modified) {
                                continue 'pairs;
                            }
                        }
                    }
                } else {
                    if seen_umis.contains(&umi) || seen_umis.iter().any(|known_umi|
                        edit_distance_bounded(
                            known_umi.as_ref(),
                            umi.as_ref(),
                            levenshtein_max as usize,
                        ).is_some()
                    ) {
                        continue;
                    }
                }
                seen_umis.insert(umi);
            }
        }

        good_records += 1;

        writer_fwr.write(
            std::str::from_utf8(rec_fwr.name()).unwrap(),
            Option::from(rec_fwr.id()),
            &rec_fwr.seq()[start_index_fwr as usize..],
            rec_fwr.qual(),
        )
            .expect("couldn't write out a forward record");
        writer_rev.write(
            std::str::from_utf8(rec_rev.name()).unwrap(),
            Option::from(rec_rev.id()),
            &rec_rev.seq()[start_index_rev as usize..],
            rec_rev.qual(),
        )
            .expect("couldn't write out a reverse record");
    }

    println!("filtered {} down to {}", pluralize("pair", total_records, true), pluralize("pair", good_records, true));
    // TODO: stop assuming forward and reverse reads appear in the proper order
    // TODO: verbose logging (masked reads, etc.)
    // TODO: exit codes
    // TODO: support .f[ast]q.gz[ip] output
    // TODO: do things on quality scores
}

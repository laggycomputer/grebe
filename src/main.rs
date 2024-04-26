use std::collections::HashSet;
use std::path::PathBuf;
use std::process::exit;

use bio::bio_types::sequence::SequenceRead;
use bio::io::fastq;
use clap::ValueHint;
use pluralizer::pluralize;

fn main() {
    let cmd = clap::command!("cristatus")
        .bin_name("cristatus")
        .about("Processing tool for Illumina sequencing data")
        .arg(clap::arg!(<"in-forward"> "forward (5'-3') reads to work with")
            .value_name("input forward .fastq")
            .value_parser(clap::value_parser!(PathBuf))
            .value_hint(ValueHint::FilePath))
        .arg(clap::arg!(<"in-reverse"> "reverse (3'-5') reads to work with")
            .value_name("input reverse .fastq")
            .value_parser(clap::value_parser!(PathBuf))
            .value_hint(ValueHint::FilePath))
        .arg(clap::arg!(--"phred64" "use the legacy phred64 encoding (over phred33) where score 0 = \"@\"")
            .required(false))
        .arg(clap::arg!(-'u' <"umi-length"> "UMI length (strip this many bases off forward reads)")
            .alias("UMI-LENGTH")
            .value_parser(0..=15)
            .required(true))
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

    let reader_fwr = match fastq::Reader::from_file(args.get_one::<PathBuf>("in-forward").unwrap()) {
        Ok(result) => result,
        Err(_) => {
            eprintln!("couldn't open input forward .fastq");
            exit(1);
        }
    };
    let reader_rev = match fastq::Reader::from_file(args.get_one::<PathBuf>("in-reverse").unwrap()) {
        Ok(result) => result,
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

    let umi_length = args.get_one::<i64>("umi-length").unwrap();
    let mut seen_umis = HashSet::new();

    let mut total_records = 0;
    let mut good_records = 0;
    let pairs = reader_fwr.records().zip(reader_rev.records());
    for (rec_fwr, rec_rev) in pairs {
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

        let all_ns = |s: &u8| -> bool { *s == ('N' as u8) };
        if rec_fwr.seq().iter().all(all_ns) || rec_rev.seq().iter().all(all_ns) {
            continue;
        }

        let umi = String::from_utf8((&rec_fwr.seq()[..*umi_length as usize]).to_vec()).unwrap();
        if !seen_umis.insert(umi) {
            continue;
        }

        good_records += 1;

        let fwr_without_umi = &rec_fwr.seq()[*umi_length as usize..];
        writer_fwr.write(std::str::from_utf8(rec_fwr.name()).unwrap(), Option::from(rec_fwr.id()), fwr_without_umi, rec_fwr.qual())
            .expect("couldn't write out a forward record");
        writer_rev.write(std::str::from_utf8(rec_rev.name()).unwrap(), Option::from(rec_rev.id()), &rec_rev.seq(), rec_rev.qual())
            .expect("couldn't write out a reverse record");
    }

    println!("filtered {} down to {}", pluralize("pair", total_records, true), pluralize("pair", good_records, true));
    // TODO: stop assuming forward and reverse reads appear in the proper order
    // TODO: verbose logging (masked reads, etc.)
    // TODO: exit codes
}

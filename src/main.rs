use std::collections::HashSet;
use std::path::PathBuf;

use bio::bio_types::sequence::SequenceRead;
use bio::io::fastq;
use clap::ValueHint;

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
            .value_parser(0..15)
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

    let reader_fwr = fastq::Reader::from_file(args.get_one::<PathBuf>("in-forward").unwrap())
        .expect("couldn't open input forward .fastq");
    let reader_rev = fastq::Reader::from_file(args.get_one::<PathBuf>("in-reverse").unwrap())
        .expect("couldn't open input reverse .fastq");

    let mut writer_fwr = fastq::Writer::to_file(args.get_one::<PathBuf>("out-forward").unwrap())
        .expect("couldn't open output forward .fastq");
    let mut writer_rev = fastq::Writer::to_file(args.get_one::<PathBuf>("out-reverse").unwrap())
        .expect("couldn't open output reverse .fastq");

    // TODO: drop masked reads

    let umi_length = args.get_one::<i64>("umi-length").expect("need valid UMI length, even 0");
    let mut seen_umis = HashSet::new();

    let mut total = 0;
    let mut count = 0;
    let pairs = reader_fwr.records().zip(reader_rev.records());
    for (rec_fwr, rec_rev) in pairs {
        total += 1;

        let rec_fwr = rec_fwr.unwrap();
        let rec_rev = rec_rev.unwrap();

        let umi = String::from_utf8((&rec_fwr.seq()[..*umi_length as usize]).to_vec()).unwrap();
        if seen_umis.contains(&umi) {
            continue
        }
        seen_umis.insert(umi);
        count += 1;

        let fwr_without_umi = &rec_fwr.seq()[*umi_length as usize..];

        writer_fwr.write(std::str::from_utf8(rec_fwr.name()).unwrap(), Option::from(rec_fwr.id()), fwr_without_umi, rec_fwr.qual())
            .expect("couldn't write out a forward record");

        writer_rev.write(std::str::from_utf8(rec_rev.name()).unwrap(), Option::from(rec_rev.id()), &rec_rev.seq(), rec_rev.qual())
            .expect("couldn't write out a reverse record");
    }

    println!("{}", total);
    println!("{}", count);

    // TODO: stop assuming forward and reverse reads appear in the proper order
}

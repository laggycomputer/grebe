fn main() {
    let cmd = clap::Command::new("cristatus")
        .bin_name("cristatus")
        .about("Processing tool for Illumina sequencing data");

    let matches = cmd.get_matches();

    println!("Hello, world!");
}

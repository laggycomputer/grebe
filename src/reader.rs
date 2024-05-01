use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, Read};
use std::path::PathBuf;
use std::process::exit;

use bio::io::fastq;
use flate2::bufread::MultiGzDecoder;

pub(crate) enum ReaderMaybeGzip {
    GZIP(BufReader<MultiGzDecoder<BufReader<File>>>),
    UNCOMPRESSED(BufReader<File>),
    NULL(BufReader<io::Empty>),
}

impl Read for ReaderMaybeGzip {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            ReaderMaybeGzip::GZIP(backer) => backer.read(buf),
            ReaderMaybeGzip::UNCOMPRESSED(backer) => backer.read(buf),
            ReaderMaybeGzip::NULL(backer) => backer.read(buf),
        }
    }
}

impl BufRead for ReaderMaybeGzip {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match self {
            ReaderMaybeGzip::GZIP(backer) => backer.fill_buf(),
            ReaderMaybeGzip::UNCOMPRESSED(backer) => backer.fill_buf(),
            ReaderMaybeGzip::NULL(backer) => backer.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            ReaderMaybeGzip::GZIP(backer) => backer.consume(amt),
            ReaderMaybeGzip::UNCOMPRESSED(backer) => backer.consume(amt),
            ReaderMaybeGzip::NULL(backer) => backer.consume(amt),
        }
    }
}

pub(crate) fn reader_maybe_gzip(path_buf: &PathBuf) -> Result<(fastq::Reader<ReaderMaybeGzip>, bool), io::Error> {
    let mut file = File::open(path_buf)?;
    let mut magic = [0; 2];
    file.read(&mut magic[..])?;

    let reopen = BufReader::new(File::open(path_buf)?);

    if magic.eq(&[0x1f, 0x8b]) {
        Ok((fastq::Reader::from_bufread(ReaderMaybeGzip::GZIP(BufReader::new(MultiGzDecoder::new(reopen)))), true))
    } else {
        Ok((fastq::Reader::from_bufread(ReaderMaybeGzip::UNCOMPRESSED(reopen)), false))
    }
}

fn reader_from_path(maybe_path_buf: Option<&PathBuf>, silent: bool) -> fastq::Reader<ReaderMaybeGzip> {
    match maybe_path_buf {
        Some(path_buf) => match reader_maybe_gzip(path_buf) {
            Ok((result, was_compressed)) => {
                if was_compressed && !silent { eprintln!("info: parsing {} as a gzip", path_buf.display()) }
                result
            }
            Err(_) => {
                eprintln!("couldn't open input {} for reading", path_buf.display());
                exit(1);
            }
        }
        None => fastq::Reader::from_bufread(ReaderMaybeGzip::NULL(BufReader::new(io::empty())))
    }
}

pub(crate) fn make_reader_pair(input_paths: (Option<&PathBuf>, Option<&PathBuf>), silent: bool)
                               -> (fastq::Reader<ReaderMaybeGzip>, fastq::Reader<ReaderMaybeGzip>) {
    (reader_from_path(input_paths.0, silent), reader_from_path(input_paths.1, silent))
}

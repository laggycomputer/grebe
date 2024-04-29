use std::fs::File;
use std::io;
use std::io::{BufReader, Read};
use std::path::PathBuf;
use std::process::exit;

use bio::io::fastq;
use flate2::bufread::MultiGzDecoder;

pub(crate) enum ReaderMaybeGzip {
    GZIP(MultiGzDecoder<BufReader<File>>),
    UNCOMPRESSED(BufReader<File>),
    NULL(io::Empty),
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

pub(crate) fn reader_maybe_gzip(path_buf: &PathBuf) -> Result<(fastq::Reader<BufReader<ReaderMaybeGzip>>, bool), io::Error> {
    let mut file = File::open(path_buf)?;
    let mut magic = [0; 2];
    file.read(&mut magic[..])?;

    let reopen = BufReader::new(File::open(path_buf)?);

    if magic.eq(&[0x1f, 0x8b]) {
        Ok((fastq::Reader::from_bufread(BufReader::new(ReaderMaybeGzip::GZIP(MultiGzDecoder::new(reopen)))), true))
    } else {
        Ok((fastq::Reader::from_bufread(BufReader::new(ReaderMaybeGzip::UNCOMPRESSED(reopen))), false))
    }
}

fn reader_from_path(path_buf: &PathBuf, silent: bool) -> fastq::Reader<BufReader<ReaderMaybeGzip>> {
    match reader_maybe_gzip(path_buf) {
        Ok((result, was_compressed)) => {
            if was_compressed && !silent { eprintln!("info: parsing {} as a gzip", path_buf.display()) }
            result
        }
        Err(_) => {
            eprintln!("couldn't open {} for reading", path_buf.display());
            exit(1);
        }
    }
}

pub(crate) fn make_reader_pair(input_paths: (&PathBuf, &PathBuf), silent: bool)
                               -> (fastq::Reader<BufReader<ReaderMaybeGzip>>,
                                   fastq::Reader<BufReader<ReaderMaybeGzip>>) {
    (reader_from_path(input_paths.0, silent), reader_from_path(input_paths.1, silent))
}

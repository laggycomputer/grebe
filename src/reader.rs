use std::fs::File;
use std::io;
use std::io::{BufReader, Read};
use std::path::PathBuf;
use bio::io::fastq;
use flate2::bufread::MultiGzDecoder;

pub(crate) enum ReaderMaybeGzip {
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

use std::fs::{File, OpenOptions};
use std::io;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use bio::io::fastq;
use flate2::Compression;
use flate2::write::GzEncoder;

pub(crate) enum WriterMaybeGzip {
    GZIP(GzEncoder<File>),
    UNCOMPRESSED(File),
}

impl Write for WriterMaybeGzip {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            WriterMaybeGzip::GZIP(backer) => backer.write(buf),
            WriterMaybeGzip::UNCOMPRESSED(backer) => backer.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            WriterMaybeGzip::GZIP(backer) => backer.flush(),
            WriterMaybeGzip::UNCOMPRESSED(backer) => backer.flush(),
        }
    }
}

pub(crate) fn writer_maybe_gzip(path_buf: &PathBuf) -> Result<(fastq::Writer<WriterMaybeGzip>, bool), io::Error> {
    let file = OpenOptions::new().write(true).create_new(true).open(path_buf)?;

    if match path_buf.extension() {
        None => false,
        Some(ext) => ext == "gzip" || ext == "gz"
    } {
        Ok((fastq::Writer::from_bufwriter(BufWriter::new(
            WriterMaybeGzip::GZIP(GzEncoder::new(file, Compression::default())))), true))
    } else {
        Ok((fastq::Writer::from_bufwriter(BufWriter::new(WriterMaybeGzip::UNCOMPRESSED(file))), false))
    }
}

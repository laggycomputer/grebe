use std::{fs, io};
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, ErrorKind, Write};
use std::os::unix::fs::MetadataExt;
use std::path::PathBuf;
use std::process::exit;

use bio::io::fastq;
use flate2::Compression;
use flate2::write::GzEncoder;

pub(crate) enum WriterMaybeGzip {
    GZIP(GzEncoder<File>),
    UNCOMPRESSED(File),
    NULL(io::Sink),
}

impl Write for WriterMaybeGzip {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            WriterMaybeGzip::GZIP(backer) => backer.write(buf),
            WriterMaybeGzip::UNCOMPRESSED(backer) => backer.write(buf),
            WriterMaybeGzip::NULL(backer) => backer.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            WriterMaybeGzip::GZIP(backer) => backer.flush(),
            WriterMaybeGzip::UNCOMPRESSED(backer) => backer.flush(),
            WriterMaybeGzip::NULL(backer) => backer.flush(),
        }
    }
}

pub(crate) fn writer_maybe_gzip(path_buf: &PathBuf) -> Result<(fastq::Writer<WriterMaybeGzip>, bool), io::Error> {
    let meta = fs::metadata(path_buf);
    if meta.is_ok() && meta.unwrap().size() > 0 {
        return Err(io::Error::other(""));
    }

    let file = OpenOptions::new().write(true).create(true).open(path_buf)?;

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

pub(crate) fn writer_from_path(path_buf: &PathBuf) -> fastq::Writer<WriterMaybeGzip> {
    match writer_maybe_gzip(path_buf) {
        Ok((result, was_compressed)) => {
            if was_compressed { eprintln!("info: writing {} as a gzip", path_buf.display()) }
            result
        }
        Err(err) => {
            match err.kind() {
                ErrorKind::Other => eprintln!("refusing to overwrite nonempty file {}", path_buf.display()),
                _ => eprintln!("couldn't open output {} for writing", path_buf.display())
            }
            exit(1);
        }
    }
}

pub(crate) fn make_writer_pair(output_paths: (&PathBuf, &PathBuf))
                               -> (fastq::Writer<WriterMaybeGzip>, fastq::Writer<WriterMaybeGzip>) {
    (writer_from_path(output_paths.0), writer_from_path(output_paths.1))
}

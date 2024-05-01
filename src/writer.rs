use std::fs::{File, OpenOptions};
use std::io;
use std::io::{BufWriter, ErrorKind, Seek, SeekFrom, Write};
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
    let mut file = OpenOptions::new().write(true).create(true).open(path_buf)?;
    if file.seek(SeekFrom::End(0)).unwrap() > 0 {
        return Err(io::Error::other(""));
    }

    file.seek(SeekFrom::Start(0))?;

    if match path_buf.extension() {
        Some(ext) if ext == "gzip" || ext == "gz" => true,
        Some(_) => false,
        None => false,
    } {
        Ok((fastq::Writer::from_bufwriter(BufWriter::new(
            WriterMaybeGzip::GZIP(GzEncoder::new(file, Compression::default())))), true))
    } else {
        Ok((fastq::Writer::from_bufwriter(BufWriter::new(WriterMaybeGzip::UNCOMPRESSED(file))), false))
    }
}

pub(crate) fn writer_from_path(maybe_path_buf: Option<&PathBuf>) -> fastq::Writer<WriterMaybeGzip> {
    match maybe_path_buf {
        Some(path_buf) => match writer_maybe_gzip(path_buf) {
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
        },
        None => fastq::Writer::from_bufwriter(BufWriter::new(WriterMaybeGzip::NULL(io::sink())))
    }
}

pub(crate) fn make_writer_pair(output_paths: (Option<&PathBuf>, Option<&PathBuf>))
                               -> (fastq::Writer<WriterMaybeGzip>, fastq::Writer<WriterMaybeGzip>) {
    (writer_from_path(output_paths.0), writer_from_path(output_paths.1))
}

use std::fs::File;

use bio::bio_types::sequence::SequenceRead;
use bio::io::fastq;

use crate::FastqPair;

pub struct PairHandler {
    pub(crate) record_writers: (fastq::Writer<File>, fastq::Writer<File>),
}

impl PairHandler {
    pub fn write_pair(&mut self, pair: FastqPair) {
        // TODO: reimplement slicing etc
        self.record_writers.0.write(
            std::str::from_utf8(pair.0.name()).unwrap(),
            Option::from(pair.0.id()),
            pair.0.seq(),
            pair.0.qual(),
        )
            .expect("couldn't write out a forward record");
        self.record_writers.1.write(
            std::str::from_utf8(pair.1.name()).unwrap(),
            Option::from(pair.1.id()),
            pair.1.seq(),
            pair.1.qual(),
        ).expect("couldn't write out a reverse record");
    }
}
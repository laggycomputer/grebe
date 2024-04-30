use bio::io::fastq;

pub(crate) type FastqPair = (fastq::Record, fastq::Record);
pub(crate) type UMIVec = Vec<u8>;

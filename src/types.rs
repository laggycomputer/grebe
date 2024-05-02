use bio::io::fastq;

use crate::writer::WriterMaybeGzip;

pub(crate) type FastqPair = (fastq::Record, fastq::Record);
pub(crate) type UMIVec = Vec<u8>;
pub(crate) type QualityVoteTotal = u32;

pub(crate) struct OutputWriters {
    pub(crate) paired: (fastq::Writer<WriterMaybeGzip>, fastq::Writer<WriterMaybeGzip>),
    pub(crate) unpaired: (fastq::Writer<WriterMaybeGzip>, fastq::Writer<WriterMaybeGzip>),
}

pub(crate) enum WhichRead {
    FORWARD,
    REVERSE,
}

pub(crate) type BaseQualityVotes = (QualityVoteTotal, QualityVoteTotal, QualityVoteTotal, QualityVoteTotal);
pub(crate) type QualityVoteVec = Vec<BaseQualityVotes>;

pub(crate) struct PairDropReasonCount {
    pub(crate) both_masked: usize,
}

impl Default for PairDropReasonCount {
    fn default() -> Self {
        PairDropReasonCount {
            both_masked: Default::default(),
        }
    }
}
use bio::io::fastq;

pub(crate) type FastqPair = (fastq::Record, fastq::Record);
pub(crate) type UMIVec = Vec<u8>;
pub(crate) type QualityVoteTotal = u32;
pub(crate) type BaseQualityVotes = (QualityVoteTotal, QualityVoteTotal, QualityVoteTotal, QualityVoteTotal);
pub(crate) type QualityVoteVec = Vec<BaseQualityVotes>;

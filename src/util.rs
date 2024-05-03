use bio::alphabets::dna;
use bio::utils::TextSlice;

fn check_primer_base(bases: (&u8, &u8)) -> bool {
    let (primer_base, seq_base) = bases;
    match primer_base.to_ascii_uppercase() {
        b'A' | b'T' | b'C' | b'G' => primer_base.eq_ignore_ascii_case(seq_base),
        b'W' => "AT".contains(seq_base.to_ascii_uppercase() as char),
        b'M' => "AC".contains(seq_base.to_ascii_uppercase() as char),
        b'R' => "AG".contains(seq_base.to_ascii_uppercase() as char),
        b'Y' => "TC".contains(seq_base.to_ascii_uppercase() as char),
        b'K' => "TG".contains(seq_base.to_ascii_uppercase() as char),
        b'S' => "CG".contains(seq_base.to_ascii_uppercase() as char),
        // already verified this is valid fully specified DNA alphabet
        b'B' => seq_base.to_ascii_uppercase() != b'A',
        b'V' => seq_base.to_ascii_uppercase() != b'T',
        b'D' => seq_base.to_ascii_uppercase() != b'C',
        b'H' => seq_base.to_ascii_uppercase() != b'G',
        // why is this in a primer
        b'N' => true,
        _ => unimplemented!()
    }
}

pub(crate) fn check_primer(primer: TextSlice, seq: TextSlice) -> Result<bool, &'static str> {
    // willfully ignore IUPAC in sequence; if it has an N or anything besides a base call that's not something we want
    // anyway. also, trust the primer is valid already
    if !dna::alphabet().is_word(&seq[..primer.len()]) {
        return Err("seq invalid");
    }

    if seq.len() < primer.len() {
        return Ok(false);
    }

    return Ok(primer.iter().zip(seq.iter().take(primer.len())).all(check_primer_base));
}
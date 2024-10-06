""" Module with auxiliary functions for bio_tools.py, run_dna_rna_tools(),
filter_fastq() """


dna = {"a", "t", "g", "c"}
rna = {"a", "u", "g", "c"}
rna_complements = {
    "a": "u",
    "u": "a",
    "g": "c",
    "c": "g",
    "A": "U",
    "U": "A",
    "C": "G",
    "G": "C",
}
dna_complements = {
    "a": "t",
    "t": "a",
    "g": "c",
    "c": "g",
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
}


def transcribe(seq: str) -> str:
    for nucleotide in seq:
        if nucleotide == "T":
            seq = seq.replace("T", "U")
        if nucleotide == "t":
            seq = seq.replace("t", "u")
    return seq


def reverse(seq: str) -> str:
    seq = seq[::-1]
    return seq


def complement(seq: str) -> str:
    complement_seq = ""
    if set(seq.lower()).issubset(dna):
        for nucleotide in seq:
            complement_seq += dna_complements[nucleotide]
    else:
        for nucleotide in seq:
            complement_seq += rna_complements[nucleotide]
    return complement_seq


def reverse_complement(seq: str) -> str:
    reverse_complement_seq = ""
    seq = seq[::-1]
    if set(seq.lower()).issubset(dna):
        for nucleotide in seq:
            reverse_complement_seq += dna_complements[nucleotide]
    else:
        for nucleotide in seq:
            reverse_complement_seq += rna_complements[nucleotide]
    return reverse_complement_seq


def triplets_count(seq: str) -> int:
    triplets = len(seq) // 3
    return triplets


def gc_content(seq: str) -> float:
    gc = (seq.lower().count("g") + seq.lower().count("c")) / len(seq) * 100
    return gc


def is_palindrom(seq: str) -> bool:
    return seq == seq[::-1]

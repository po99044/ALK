from __future__ import annotations

from typing import Dict

_EXON_REF: Dict[str, str] = {
    "ALK_E20" : "accacccacctgcagtgtaccgccggaagcaccaggagctgcaagccatgcagatggagctgcagagccctgagtacaagctgagcaagctccgcacctcgaccatcatgaccgactacaaccccaactactgctttgctggcaagacctcctccatcagtgacctgaaggaggtgccgcggaaaaacatcaccctcattcggtgagcgccc".upper(),

    "ALK_E21" : "ttctcctttgcacaggggtctgggccatggcgcctttggggaggtgtatgaaggccaggtgtccggaatgcccaacgacccaagccccctgcaagtggctgtgaaggtaagaagtg".upper(),

    "ALK_E22" : "cccttctctgcccagacgctgcctgaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcaggtaaagccac".upper(),

    "ALK_E23" : "ctctctgctctgcagcaaattcaaccaccagaacattgttcgctgcattggggtgagcctgcaatccctgccccggttcatcctgctggagctcatggcggggggagacctcaagtccttcctccgagagacccgccctcgcccggtgagtgaga".upper(),

    "ALK_E24" : "tctgtctccccacagagccagccctcctccctggccatgctggaccttctgcacgtggctcgggacattgcctgtggctgtcagtatttggaggaaaaccacttcatccaccggtgagtcaaa".upper(),

    "ALK_E25" : "ttcctttcttcccagagacattgctgccagaaactgcctcttgacctgtccaggccctggaagagtggccaagattggagacttcgggatggcccgagacatctacaggtgagtaaag".upper(),

    "ALK_E26" : "tctccttccccacagggcgagctactatagaaagggaggctgtgccatgctgccagttaagtggatgcccccagaggccttcatggaaggaatattcacttctaaaacagacacatggtaagtcagc".upper(),
    "ALK_E27" : "ctgtcccatgcccaggtcctttggagtgctgctatgggaaatcttttctcttggatatatgccataccccagcaaaagcaaccaggaagttctggagtttgtcaccagtggaggccggatggacccacccaagaactgccctgggcctgtgtatgactct".upper(),
    "ALK_E28" : "tgcttcttcttttagataccggataatgactcagtgctggcaacatcagcctgaagacaggcccaactttgccatcattttggagaggattgaatactgcacccaggtaaaacatt".upper()

}



def exon(a: int) -> str:
    """
    Return uppercase reference sequence for ALK exon a (20–28).
    Used to define FP seed and expected target length for FASTQ trimming.
    """
    key = f"ALK_E{int(a)}"
    try:
        return _EXON_REF[key]
    except KeyError as e:
        raise ValueError(f"Unsupported exon: {a}. Must be one of 20–28.") from e

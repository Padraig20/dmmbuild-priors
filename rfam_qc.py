
# rfam_qc.py
"""
SYNOPSIS:

Compute quality control measures for each sequence in Rfam SEED
files.

DESCRIPTION:

Adapted from Paul P. Gardner's original Perl script.

Outputs one tab-separated line per family with columns:

1  FAMILY (GF AC if present, else GF ID)
2  MEAN_FRACTN_CANONICAL_BPs
3  COVARIATION
4  NO_SEQs
5  ALN_LENGTH
6  NO_BPs
7  NO_NUCs
8  mean_PID
9  max_PID
10 min_PID
11 mean_LEN
12 max_LEN
13 min_LEN
14 FRACTN_NUCs
15 FRAC_A
16 FRAC_C
17 FRAC_G
18 FRAC_U
19 MAX_DINUC   (one of R,Y,K,W,S,M with its value, e.g. "S:0.567")
20 CG_CONTENT  (C+G)

EXAMPLES:

Run it on the Rfam 10.0 SEED file to generate the corresponding Rfam.qc file:
python rfam_qc.py Rfam_10.0/Rfam.seed >> ./Rfam.qc

AUTHOR:

Patrick Styll

COPYRIGHT:

This program is free software. You may copy, modify, and redistribute
it under the same terms as Python itself.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple, Iterable, Optional

# IUPAC RNA nucleotide codes (Uracil replaces Thymine)
IUPAC_BASES = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "U": {"U"},
    "T": {"U"}, # T -> U (shouldnt appear in RNA though)
    "R": {"A", "G"},
    "Y": {"C", "U"},
    "S": {"C", "G"},
    "W": {"A", "U"},
    "K": {"G", "U"},
    "M": {"A", "C"},
    "B": {"C", "G", "U"},
    "D": {"A", "G", "U"},
    "H": {"A", "C", "U"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "U"},
}

def is_nucleotide(ch: str) -> bool:
    ch = ch.upper()
    return ch in IUPAC_BASES

def is_complementary(a: str, b: str) -> bool:
    a  = a.upper()
    b  = b.upper()
    ab = a + b
    return ab in {
        "AU", "AT",
        "UA", "TA",
        "CG", "GC",
        "GU", "GT",
        "UG", "TG",
        "RY", "YR",
        "MK", "KM",
        "SS", "WW",
    }

_OPEN_TO_CLOSE = {
    "(": ")",
    "[": "]",
    "{": "}",
    "<": ">",
}
# allow letters as well
for o, c in zip("ABCDEFGHIJKLMNOPQRSTUVWXYZ", "abcdefghijklmnopqrstuvwxyz"):
    _OPEN_TO_CLOSE[o] = c

_CLOSE_TO_OPEN = {v: k for k, v in _OPEN_TO_CLOSE.items()}

def parse_pairs_from_ss(ss: str) -> List[Tuple[int, int]]:
    stacks: Dict[str, List[int]] = defaultdict(list)
    pairs: List[Tuple[int, int]] = []
    for i, ch in enumerate(ss):
        if ch in _OPEN_TO_CLOSE:
            stacks[ch].append(i)
        elif ch in _CLOSE_TO_OPEN:
            o = _CLOSE_TO_OPEN[ch]
            if stacks[o]:
                j = stacks[o].pop()
                pairs.append((j, i))
            else:
                # unbalanced close => ignore
                pass
        else:
            # '.', gaps, etc. => ignore
            pass
    # unbalanced opens => ignore
    pairs.sort()
    return pairs

def comb2(n: int) -> int:
    return n * (n - 1) // 2

def pairwise_identity(seq1: str, seq2: str) -> float:
    """
    match = count positions where seq1 has a nucleotide and equals seq2
    len1  = count nucleotides in seq1
    len2  = count nucleotides in seq2
    pid   = match / max(len1, len2)
    """
    match = 0
    len1  = 0
    len2  = 0
    for a, b in zip(seq1, seq2):
        if is_nucleotide(a):
            len1 += 1
        if is_nucleotide(b):
            len2 += 1
        if is_nucleotide(a) and a.upper().replace("T", "U") == b.upper().replace("T", "U"):
            match += 1
    denom = max(len1, len2)
    return (match / denom) if denom else 0.0

def add_iupac_counts(ch: str, counts: Dict[str, float]) -> None:
    ch = ch.upper()
    if ch == "T":
        ch = "U"
    if ch not in IUPAC_BASES:
        return
    bases = IUPAC_BASES[ch]
    w = 1.0 / len(bases)
    for b in bases:
        counts[b] += w

@dataclass
class AlignmentBlock:
    fam: str
    seqs: Dict[str, str]  # name->aligned sequence
    ss_cons: str

def stockholm_blocks(handle: Iterable[str]) -> Iterable[AlignmentBlock]:
    fam_ac: Optional[str] = None
    fam_id: Optional[str] = None
    ss_parts: List[str]   = []
    seq_parts: Dict[str, List[str]] = defaultdict(list)

    def flush():
        nonlocal fam_ac, fam_id, ss_parts, seq_parts
        if not seq_parts:
            # empty block
            fam_ac    = fam_id = None
            ss_parts  = []
            seq_parts = defaultdict(list)
            return None

        fam     = fam_ac or fam_id or "UNKNOWN"
        seqs    = {k: "".join(v) for k, v in seq_parts.items()}
        ss_cons = "".join(ss_parts)
        blk     = AlignmentBlock(fam=fam, seqs=seqs, ss_cons=ss_cons)

        fam_ac    = fam_id = None
        ss_parts  = []
        seq_parts = defaultdict(list)
        return blk

    for raw in handle:
        line = raw.rstrip("\n")
        if not line:
            continue
        if line.strip() == "//":
            blk = flush()
            if blk is not None:
                yield blk
            continue
        if line.startswith("#=GF"):
            # e.g. "#=GF AC   RF00001"
            parts = line.split()
            if len(parts) >= 3:
                tag = parts[1]
                val = parts[2]
                if tag == "AC":
                    fam_ac = val
                elif tag == "ID": # only for fallback...
                    fam_id = val
            continue
        if line.startswith("#=GC"):
            parts = line.split()
            if len(parts) >= 3 and parts[1] == "SS_cons":
                ss_parts.append(parts[2])
            continue
        if line.startswith("#"):
            continue

        # sequence line: "<name> <aligned_seq_fragment>"
        parts = line.split()
        if len(parts) >= 2:
            name = parts[0]
            frag = parts[1]
            seq_parts[name].append(frag)

def compute_qc(block: AlignmentBlock) -> str:
    seqs  = block.seqs
    names = list(seqs.keys())
    nseq  = len(names)
    if nseq == 0:
        return ""

    # alignment length
    aln_len = len(next(iter(seqs.values())))

    ss_cons = block.ss_cons
    pairs   = parse_pairs_from_ss(ss_cons) if ss_cons else []
    n_bp    = len(pairs)

    paired_positions = 2 * n_bp # $paired

    # per-sequence canonical bp fraction + length stats + nucleotide counts
    mean_canon      = 0.0
    lens: List[int] = []
    total_nucs      = 0
    base_counts     = {"A": 0.0, "C": 0.0, "G": 0.0, "U": 0.0}

    seq_list = [seqs[n] for n in names]

    for s in seq_list:
        # ungapped length (count nucleotides including ambiguity codes)
        L = sum(1 for ch in s if is_nucleotide(ch))
        lens.append(L)
        total_nucs += L
        for ch in s:
            if is_nucleotide(ch):
                add_iupac_counts(ch, base_counts)

        if paired_positions == 0:
            mean_canon += 1.0
            continue

        consistent = 0
        for i, j in pairs:
            a = s[i] if i < len(s) else "-"
            b = s[j] if j < len(s) else "-"
            if is_complementary(a, b):
                consistent += 1
        mean_canon += (2.0 * consistent / paired_positions) if paired_positions else 1.0

    mean_canon /= nseq

    mean_len = sum(lens) / nseq
    max_len  = max(lens)
    min_len  = min(lens)

    frac_nucs = (total_nucs / (nseq * aln_len)) if (nseq and aln_len) else 0.0

    # A/C/G/U fractions over all nucleotides (ambiguity codes contribute fractionally)
    if total_nucs > 0:
        frac_a = base_counts["A"] / total_nucs
        frac_c = base_counts["C"] / total_nucs
        frac_g = base_counts["G"] / total_nucs
        frac_u = base_counts["U"] / total_nucs
    else:
        frac_a = frac_c = frac_g = frac_u = 0.0

    # MAX_DINUC
    # Which 2-letter nucleotide class (out of the six binary IUPAC classes)
    # has for the largest fraction of nucleotides in the alignment?
    dinuc_groups = {
        "R": frac_a + frac_g,  # purines
        "Y": frac_c + frac_u,  # pyrimidines
        "K": frac_g + frac_u,  # keto
        "W": frac_a + frac_u,  # weak (A/U)
        "S": frac_c + frac_g,  # strong (C/G)
        "M": frac_a + frac_c   # amino
    }
    max_code      = max(dinuc_groups, key=lambda k: dinuc_groups[k])
    max_dinuc_val = dinuc_groups[max_code]
    cg_content    = frac_c + frac_g

    # pairwise identity stats (PID; mean/max/min)
    if nseq < 2:
        mean_pid = max_pid = min_pid = 0.0
    else:
        total_pairs = 0
        pid_sum     = 0.0
        max_pid     = 0.0
        min_pid     = 1.0
        for i in range(nseq):
            si = seq_list[i]
            for j in range(i + 1, nseq):
                pj = pairwise_identity(si, seq_list[j])
                pid_sum     += pj
                total_pairs += 1
                if pj > max_pid:
                    max_pid = pj
                if pj < min_pid:
                    min_pid = pj
        mean_pid = pid_sum / total_pairs if total_pairs else 0.0
        if total_pairs == 0:
            max_pid = min_pid = 0.0

    # covariation metric (alicovar)
    alitot   = 0
    alicon   = 0
    aliincon = 0

    if pairs and nseq >= 2:
        for i, j in pairs:
            # iterate over all unordered sequence pairs (s1,s2)
            for a_idx in range(nseq - 1):
                s1 = seq_list[a_idx]
                c1 = s1[i]
                d1 = s1[j]
                for b_idx in range(a_idx + 1, nseq):
                    s2 = seq_list[b_idx]
                    c2 = s2[i]
                    d2 = s2[j]

                    if is_nucleotide(c1) and is_nucleotide(c2) and is_nucleotide(d1) and is_nucleotide(d2):
                        alitot += 2
                        if is_complementary(c1, d1) and is_complementary(c2, d2):
                            if c1 != c2:
                                alicon += 1
                            if d1 != d2:
                                alicon += 1
                        else:
                            if c1 != c2:
                                aliincon += 1
                            if d1 != d2:
                                aliincon += 1

    covariation = ((alicon - aliincon) / alitot) if alitot else 0.0

    out = [
        block.fam,
        f"{mean_canon:.5f}",
        f"{covariation:.5f}",
        str(nseq),
        str(aln_len),
        str(n_bp),
        str(total_nucs),
        f"{mean_pid:.3f}",
        f"{max_pid:.3f}",
        f"{min_pid:.3f}",
        f"{mean_len:.3f}",
        str(max_len),
        str(min_len),
        f"{frac_nucs:.3f}",
        f"{frac_a:.3f}",
        f"{frac_c:.3f}",
        f"{frac_g:.3f}",
        f"{frac_u:.3f}",
        f"{max_code}:{max_dinuc_val:.3f}",
        f"{cg_content:.3f}",
    ]
    return "\t".join(out)

def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description="Create an Rfam-style QC file from concatenated Stockholm seed alignments.")
    p.add_argument("stockholm", help="Input Stockholm file containing multiple alignments.")
    p.add_argument("--header", action="store_true", help="Write a header line.")
    args = p.parse_args(argv)

    inp = open(args.stockholm, "rt", encoding="utf-8", errors="replace")

    try:
        if args.header:
            print(
                "FAMILY\tMEAN_FRACTN_CANONICAL_BPs\tCOVARIATION\tNO_SEQs\tALN_LENGTH\tNO_BPs\tNO_NUCs\t"
                "mean_PID\tmax_PID\tmin_PID\tmean_LEN\tmax_LEN\tmin_LEN\tFRACTN_NUCs\t"
                "FRAC_A\tFRAC_C\tFRAC_G\tFRAC_U\tMAX_DINUC\tCG_CONTENT\n"
            )

        for blk in stockholm_blocks(inp):
            line = compute_qc(blk)
            if line:
                print(line)
    finally:
        inp.close()

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
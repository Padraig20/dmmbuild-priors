import os
import ctypes
import numpy as np
import matplotlib.pyplot as plt

eps = 1e-10

# Setting up C-Types for calling the dummer-build shared library
lib = ctypes.CDLL(os.path.abspath("dummer-build.so"))

# Struct mirrors
class GapPriors(ctypes.Structure):
    _fields_ = [
        ("match",     ctypes.c_double),
        ("insStart",  ctypes.c_double),
        ("delStart",  ctypes.c_double),
        ("insEnd",    ctypes.c_double),
        ("insExtend", ctypes.c_double),
        ("delEnd",    ctypes.c_double),
        ("delExtend", ctypes.c_double),
    ]
    
class ArrayDouble(ctypes.Structure):
    _fields_ = [
        ("data", ctypes.POINTER(ctypes.c_double)),
        ("len", ctypes.c_size_t),
    ]

# Signatures
lib.get_relative_entropies_c.argtypes = [
    ctypes.c_char_p,         # filename
    ctypes.c_double,         # symfrac
    ctypes.c_double,         # ere
    ctypes.c_double,         # bwMaxiter
    ctypes.c_double,         # bwMaxDiff
    ctypes.c_bool,           # countOnly  (matches C _Bool)
    GapPriors,               # by value
]
lib.get_relative_entropies_c.restype = ArrayDouble

# Setting up the initial gap priors

# These are the initial gap priors for protein sequences,
# originally trained by Graeme Mitchison in the mid-1990's
gp = GapPriors(
    match=0.7939, insStart=0.0278, delStart=0.0135,
    insEnd=0.1551, insExtend=0.1331,
    delEnd=0.9002, delExtend=0.5630
)

# For nucleotide sequences, trained on a portion of
# the rmark dataset
#gp = GapPriors(
#    match=2.0, insStart=0.1, delStart=0.1,
#    insEnd=0.12, insExtend=0.4,
#    delEnd=0.5, delExtend=1.0
#)

# All *.sto or *.stk files in the current directory
def get_msa_filenames(path: str) -> list[str]:
    return [os.path.join(path, f) for f in os.listdir(path) if f.endswith(".sto") or f.endswith(".stk")]

def get_relative_entropies(msas: list[str], priors: GapPriors) -> list[float]:
    results = []
    for msa in msas:
        res = lib.get_relative_entropies_c(
            msa.encode('utf-8'),
            0.50,          # symfrac
            0.59,          # ere (protein)
            1.0, 1e-6,     # bwMaxiter, bwMaxDiff
            True,          # countOnly -> true for testing!
            priors
        )
        results += res.data[:res.len]
    return results

def do_plot(entropies):
    plt.figure(figsize=(8, 5))
    plt.hist(entropies, bins=150, alpha=0.8, label="Relative Entropy")
    plt.xlabel("Relative Entropy")
    plt.ylabel("Count")
    plt.title("Relative Entropies (Count-Based)")
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    plt.figure(figsize=(8, 3))
    box = plt.boxplot(entropies, vert=False, patch_artist=True, widths=0.8,
                      boxprops=dict(facecolor='#87CEEB', color='#1E90FF'),
                      medianprops=dict(color='#FF4500', linewidth=2),
                      whiskerprops=dict(color='#1E90FF'),
                      capprops=dict(color='#1E90FF'), # <- obvious ChatGPT code ;)
                      flierprops=dict(marker='o', markerfacecolor='#FFD700', markersize=6, alpha=0.5))
    plt.xlabel("Relative Entropy (Count-Based)", fontsize=12)
    plt.title("Relative Entropies", fontsize=14, fontweight='bold')
    plt.grid(axis='x', linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    msas = get_msa_filenames("alignments")
    print(f"Found {len(msas)} MSA files: {msas}")

    entropies = np.array(get_relative_entropies(msas, gp))
    do_plot(entropies)

    entropies = entropies[entropies < 1000]
    do_plot(entropies)
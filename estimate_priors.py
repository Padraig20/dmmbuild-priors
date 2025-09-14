import os
import ctypes
import dirichlet
import numpy as np

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
lib.build_hmm_c.argtypes = [
    ctypes.c_char_p,         # filename
    ctypes.c_double,         # symfrac
    ctypes.c_double,         # ere
    ctypes.c_double,         # esigma
    ctypes.c_double,         # bwMaxiter
    ctypes.c_double,         # bwMaxDiff
    ctypes.c_bool,           # countOnly  (matches C _Bool)
    GapPriors,               # by value
]
lib.build_hmm_c.restype = ArrayDouble

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

def get_transition_probabilities(msas: list[str], priors: GapPriors) -> list[list[float]]:
    results = []
    for msa in msas:
        res = lib.build_hmm_c(
            msa.encode('utf-8'),
            0.50,          # symfrac
            0.59, 45.0,    # ere (protein), esigma
            100.0, 1e-6,   # bwMaxiter, bwMaxDiff
            False,         # countOnly -> true for testing!
            priors
        )
        for i in range(int(res.len/7)):
            probs = []
            for j in range(7):
                probs += [res.data[i*7 + j]]
            results.append(probs)
    return results

def get_new_gap_priors(probs: list[list[float]]) -> GapPriors:

    probs = np.array(probs) + eps # avoid zeros for dirichlet MLE

    a_d_g = np.array([prob[0:3] for prob in probs])
    b_bp  = np.array([prob[3:5] for prob in probs])
    e_ep  = np.array([prob[5:7] for prob in probs])
    
    try:
        a = dirichlet.mle(a_d_g)
        b = dirichlet.mle(b_bp)
        e = dirichlet.mle(e_ep)
    except Exception as ex:
        print("-"*100)
        print(f"Exception during Dirichlet MLE: {ex}")
        print(f"Most likely we're simply done: Returning previous gap priors...")
        print("-"*100)
        return gp
    
    return GapPriors(
        match=a[0], insStart=a[1], delStart=a[2],
        insEnd=b[0], insExtend=b[1],
        delEnd=e[0], delExtend=e[1]
    )
    
def get_max_diff(gp1: GapPriors, gp2: GapPriors) -> float:
    diffs = [
        abs(gp1.match     - gp2.match),
        abs(gp1.insStart  - gp2.insStart),
        abs(gp1.delStart  - gp2.delStart),
        abs(gp1.insEnd    - gp2.insEnd),
        abs(gp1.insExtend - gp2.insExtend),
        abs(gp1.delEnd    - gp2.delEnd),
        abs(gp1.delExtend - gp2.delExtend),
    ]
    return max(diffs)

def print_gap_priors(gp: GapPriors):
    print(f" match:     {gp.match:.6f}")
    print(f" insStart:  {gp.insStart:.6f}")
    print(f" delStart:  {gp.delStart:.6f}")
    print(f" insEnd:    {gp.insEnd:.6f}")
    print(f" insExtend: {gp.insExtend:.6f}")
    print(f" delEnd:    {gp.delEnd:.6f}")
    print(f" delExtend: {gp.delExtend:.6f}")

if __name__ == "__main__":
    msas = get_msa_filenames("alignments")
    print(f"Found {len(msas)} MSA files: {msas}")
    
    max_iterations = 100
    
    for iteration in range(max_iterations):
        
        probs = get_transition_probabilities(msas, gp)
        new_gp = get_new_gap_priors(probs)
        
        max_diff = get_max_diff(gp, new_gp)
        
        gp = new_gp
        
        print(f"Iteration {iteration+1}:")
        print_gap_priors(gp)
        
        if max_diff < eps:
            print("Converged.")
            break
    
    print()
    print("="*100)
    print()
    print(f"Final gap priors after {iteration+1} iterations:")
    print_gap_priors(gp)
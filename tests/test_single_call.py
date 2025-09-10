import ctypes
import os

lib = ctypes.CDLL(os.path.abspath("dummer-build.so"))

# ---- Struct mirrors ----
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

# ---- Signatures ----
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

# ---- MAIN ----
gp = GapPriors(
    match=0.7939, insStart=0.0278, delStart=0.0135,
    insEnd=0.1551, insExtend=0.1331,
    delEnd=0.9002, delExtend=0.5630
)

res = lib.build_hmm_c(
    b"alignments/fn3.sto",
    0.5,           # symfrac
    0.59, 45.0,    # ere, esigma
    1000.0, 1e-6,  # bwMaxiter, bwMaxDiff
    False,         # countOnly
    gp
)

for i in range(res.len):
    if i%7 == 6:
        print(f" {res.data[i]:.6f}\n")
    else:
        print(f" {res.data[i]:.6f}", end="")
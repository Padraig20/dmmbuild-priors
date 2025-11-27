import os
import dirichlet
import numpy as np
import subprocess
import tempfile

eps = 1e-3

class GapPriors:
    def __init__(self, match: float, insStart: float, delStart: float,
                 insEnd: float, insExtend: float,
                 delEnd: float, delExtend: float):
        self.match     = match
        self.insStart  = insStart
        self.delStart  = delStart
        self.insEnd    = insEnd
        self.insExtend = insExtend
        self.delEnd    = delEnd
        self.delExtend = delExtend

def get_hmm_data(msa: str, method: str, use_counts_only: bool) -> list[list[float]]:
    # write HMM to a temporary file using the priors passed in
    out_hmm = tempfile.NamedTemporaryFile(delete=False, suffix=".hmm")
    out_hmm.close()

    cmd = [
        "./bin/dummer-build",
        "--maxiter", "100",
        "--pnone",
        msa
    ]
    
    if method == "easel":
        cmd.insert(1, "--counts")
    
    if use_counts_only:
        cmd.insert(1, "--countonly")
          # end of one HMM; keep going to read the ne
    try:
        # run the external build and write its stdout into the temporary HMM file
        with open(out_hmm.name, "w") as fh_out:
            subprocess.run(cmd, check=True, stdout=fh_out, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"hmmbuild failed: {e}")
    
    probs = []
    with open(out_hmm.name, "r") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s == "//":          # end of one HMM; keep going to read the next
                continue
            parts = s.split()
            if len(parts) != 7:
                continue
            try:
                vals = np.array([float(v) for v in parts], dtype=float)
            except ValueError:
                # skip non-numeric lines
                continue
            # file stores negative log-probabilities => convert to normal probs
            if method == "dirichlet":
                p = np.exp(-vals)
            else:  # method == "easel"
                p = np.array([vals[0], vals[3], vals[2], vals[3], vals[4], vals[5], vals[6]])
            probs.append(p.tolist())

    # clean up temporary HMM file
    try:
        os.unlink(out_hmm.name)
    except OSError:
        pass

    print(f"Read {len(probs)} transition probability rows from HMM.")

    return probs

def get_new_gap_priors(probs: list[list[float]]) -> GapPriors:

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
        print(f"No idea what happened. Terminating...")
        print("-"*100)
        exit(1)
    
    return GapPriors(
        match=a[0], insStart=a[1], delStart=a[2],
        insEnd=b[0], insExtend=b[1],
        delEnd=e[0], delExtend=e[1]
    )

def get_new_gap_priors_easel(probs: list[list[float]]) -> GapPriors:

    # write probs to temporary CSV files
    tmp_m = tempfile.NamedTemporaryFile(delete=False, suffix="_m.csv")
    tmp_m.close()
    tmp_i = tempfile.NamedTemporaryFile(delete=False, suffix="_i.csv")
    tmp_i.close()
    tmp_d = tempfile.NamedTemporaryFile(delete=False, suffix="_d.csv")
    tmp_d.close()

    np.savetxt(tmp_m.name, np.array(probs)[:,:3], delimiter=" ", fmt="%.8f")
    np.savetxt(tmp_i.name, np.array(probs)[:,3:5], delimiter=" ", fmt="%.8f")
    np.savetxt(tmp_d.name, np.array(probs)[:,5:7], delimiter=" ", fmt="%.8f")

    tmp_out_m = tempfile.NamedTemporaryFile(delete=False, suffix="_m_out")
    tmp_out_m.close()
    tmp_out_i = tempfile.NamedTemporaryFile(delete=False, suffix="_i_out")
    tmp_out_i.close()
    tmp_out_d = tempfile.NamedTemporaryFile(delete=False, suffix="_d_out")
    tmp_out_d.close()

    # run Easel's dirichlet tool on each file
    cmd_m = ["./easel/miniapps/esl-mixdchlet", "fit", "-s", "7", "1", "3", tmp_m.name, tmp_out_m.name]
    cmd_i = ["./easel/miniapps/esl-mixdchlet", "fit", "-s", "7", "1", "2", tmp_i.name, tmp_out_i.name]
    cmd_d = ["./easel/miniapps/esl-mixdchlet", "fit", "-s", "7", "1", "2", tmp_d.name, tmp_out_d.name]

    try:
        subprocess.run(cmd_m, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run(cmd_i, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run(cmd_d, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Easel dirichlet tool failed: {e}")
    print("Successfully fit Dirichlet mixtures using Easel.")

    def _parse_alphas(path: str, expected_dim: int) -> np.ndarray:
        with open(path, "r") as fh:
            lines = [ln.strip() for ln in fh if ln.strip()]
        if len(lines) < 2:
            raise RuntimeError(f"Unexpected output format in {path}: need 2 lines")
        parts = [float(x) for x in lines[1].split()]
        if len(parts) != expected_dim + 1:
            raise RuntimeError(
                f"Unexpected number of fields in {path}: got {len(parts)}, expected {expected_dim + 1}"
            )
        # drop the first field (mixture weight), keep the alphas
        return np.array(parts[1:], dtype=float)

    a = _parse_alphas(tmp_out_m.name, 3)  # match, insStart, delStart
    b = _parse_alphas(tmp_out_i.name, 2)  # insEnd, insExtend
    e = _parse_alphas(tmp_out_d.name, 2)  # delEnd, delExtend

    # cleanup temporary files
    for p in (
        tmp_m.name, tmp_i.name, tmp_d.name,
        tmp_out_m.name, tmp_out_i.name, tmp_out_d.name
    ):
        try:
            os.unlink(p)
        except OSError:
            pass

    return GapPriors(
        match=a[0], insStart=a[1], delStart=a[2],
        insEnd=b[0], insExtend=b[1],
        delEnd=e[0], delExtend=e[1]
    )

def print_gap_priors(gp: GapPriors):
    print(f" match:     {gp.match:.6f}")
    print(f" insStart:  {gp.insStart:.6f}")
    print(f" delStart:  {gp.delStart:.6f}")
    print(f" insEnd:    {gp.insEnd:.6f}")
    print(f" insExtend: {gp.insExtend:.6f}")
    print(f" delEnd:    {gp.delEnd:.6f}")
    print(f" delExtend: {gp.delExtend:.6f}")

if __name__ == "__main__":
    msa = "../Pfam-A.filtered_lt50.stk"
    print(f"Looking at file: {msa}")

    # set this to "dirichlet" to use Dirichlet MLE
    # set this to "easel" to use Easel's built-in method
    method = "easel"

    # set this to true to run without Baum-Welch
    counts_only = False

    # set to true if mean counts should be displayed
    verbose = False
                
    probs = get_hmm_data(msa, method, counts_only)

    if verbose:
        # compute and print mean of each of the 7 probability columns
        arr = np.array(probs, dtype=float)
        if arr.size == 0:
            print("No probability rows to average.")
        else:
            mean_probs = arr.mean(axis=0)
            names = ["match", "insStart", "delStart", "insEnd", "insExtend", "delEnd", "delExtend"]
            print()
            print("\n".join(f"{n}: {v:.6f}" for n, v in zip(names, mean_probs)))
            print()

    if method == "dirichlet":
        gp = get_new_gap_priors(probs)
    elif method == "easel":
        gp = get_new_gap_priors_easel(probs)
    else:
        raise ValueError(f"Unknown method: {method}")
    
    print()
    print("="*100)
    print()
        
    print_gap_priors(gp)
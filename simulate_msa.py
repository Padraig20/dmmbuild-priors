import random
import numpy as np
import matplotlib.pyplot as plt

def simulate_msa(
    n_seqs=10,
    length=50,
    deletion_start_prob=0.05,
    deletion_extend_prob=0.7,
    insertion_start_prob=0.05,
    insertion_extend_prob=0.5,
    alphabet="ACDEFGHIKLMNPQRSTVWY",
    seed=None,
    plot=False,
):
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    base_seq = [random.choice(alphabet) for _ in range(length)]
    msa = []

    for s in range(n_seqs):
        seq = base_seq.copy()
        new_seq = []
        i = 0
        while i < len(seq):
            # 'deletion' block
            if random.random() < deletion_start_prob:
                new_seq.append("-")
                i += 1
                while i < len(seq) and random.random() < deletion_extend_prob:
                    new_seq.append("-")
                    i += 1
            # 'insertion' block
            elif random.random() < insertion_start_prob:
                ins_run = [random.choice(alphabet)]
                while random.random() < insertion_extend_prob:
                    ins_run.append(random.choice(alphabet))
                new_seq.extend(ins_run)
            # 'match' block
            else:
                new_seq.append(seq[i])
                i += 1
        msa.append("".join(new_seq))

    # equalize alignment lengths (...pad shorter ones, correct?)
    max_len = max(len(s) for s in msa)
    msa = [s.ljust(max_len, "-") for s in msa]

    print("#=GF ID Simulated_MSA")
    for i, s in enumerate(msa):
        print(f"seq{i+1:02d}  {s}")
    print("//")

    if plot:
        mat = np.array([[1 if c == '-' else 0 for c in s] for s in msa])
        plt.imshow(mat, cmap="Reds", interpolation="nearest")
        plt.xlabel("Alignment position")
        plt.ylabel("Sequence index")
        plt.title("Gap map (red = gap)")
        plt.show()

    return msa

if __name__ == "__main__":

    #simulate_msa(
    #    n_seqs=20,
    #    length=50,
    #    deletion_start_prob=0.05,
    #    deletion_extend_prob=0.9,  # long gaps
    #    insertion_start_prob=0.02,
    #    insertion_extend_prob=0.5,
    #    seed=42,
    #    plot=True,
    #)

    #for i in range(50):
    #    simulate_msa(
    #        n_seqs=40,
    #        length=300,
    #        deletion_start_prob=0.1,
    #        deletion_extend_prob=0.9,  # long gaps
    #        insertion_start_prob=0.02,
    #        insertion_extend_prob=0.5,
    #        seed=i,
    #    )


    #for i in range(50):
    #    simulate_msa(
    #        n_seqs=40,
    #        length=300,
    #        deletion_start_prob=0.1,
    #        deletion_extend_prob=0.2,  # short gaps
    #        insertion_start_prob=0.02,
    #        insertion_extend_prob=0.5,
    #        seed=i,
    #    )

    #for i in range(50):
    #    simulate_msa(
    #        n_seqs=40,
    #        length=300,
    #        deletion_start_prob=0.1,
    #        deletion_extend_prob=0.9,  # long gaps
    #        insertion_start_prob=0.0,
    #        insertion_extend_prob=0.0,
    #        seed=i,
    #    )


    for i in range(50):
        simulate_msa(
            n_seqs=40,
            length=300,
            deletion_start_prob=0.0,
            deletion_extend_prob=0.0,  # short gaps
            insertion_start_prob=0.7,
            insertion_extend_prob=0.2,
            seed=i,
        )
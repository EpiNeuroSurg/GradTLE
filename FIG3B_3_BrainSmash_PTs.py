import numpy as np
from scipy.stats import pearsonr
from brainsmash.mapgen.memmap import txt2memmap
from brainsmash.mapgen.sampled import Sampled
import os

DATA_DIR = "/host/verges/tank/data/mo/demo/BrainSMASH"

# load PET
y = np.loadtxt(f"{DATA_DIR}/pet_val_pt_3B.txt")
print("Loaded PET:", y.shape)

# distance
distmat = np.load(f"{DATA_DIR}/distmat.npy", mmap_mode='r')
index = np.load(f"{DATA_DIR}/index.npy")

# ===== FIXED SURROGATE GENERATION =====
n_perm = 10000   # start smaller for speed

gen = Sampled(y, distmat, index)

# ✅ ALWAYS generate fresh (no loading broken file)
surrogates = np.zeros((n_perm, len(y)))

for i in range(n_perm):
    if i % 50 == 0:
        print(f"Progress: {i}/{n_perm}")

    surrogates[i] = gen(1)[0]

# ===== ANALYSIS =====
np.random.seed(42)
idx = np.random.choice(len(y), size=min(8000, len(y)), replace=False)

y_sub = y[idx]

for xf in [
    f"{DATA_DIR}/PTs_CorThick.txt",
    f"{DATA_DIR}/PTs_CorCurv.txt",
    f"{DATA_DIR}/PTs_MidFLAIR.txt",
    f"{DATA_DIR}/PTs_Myelin.txt"
]:
    x = np.loadtxt(xf)
    x_sub = x[idx]

    r_emp, _ = pearsonr(x_sub, y_sub)

    s_sub = surrogates[:, idx]

    x_z = (x_sub - x_sub.mean()) / x_sub.std()
    s_z = (s_sub - s_sub.mean(axis=1, keepdims=True)) / s_sub.std(axis=1, keepdims=True)

    r_null = (s_z @ x_z) / len(x_z)

    p = (np.sum(np.abs(r_null) >= abs(r_emp)) + 1) / (len(r_null) + 1)

    print("====", os.path.basename(xf))
    print("r =", round(r_emp, 4), "p =", p)

print("Done.")

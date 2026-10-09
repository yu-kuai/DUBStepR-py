# Generate stress-test inputs: subsampled real data + synthetic clustered data, varied params.
import numpy as np, scipy.io, scipy.sparse as sp, json
from pathlib import Path
import sys
REF = str(Path(__file__).resolve().parents[1] / "data") + "/"
OUT = sys.argv[1].rstrip("/") + "/"
X = sp.load_npz(REF + "vignette_pbmc1k_lognorm.npz").tocsc()
G = np.array(open(REF + "vignette_pbmc1k_genes.txt").read().split("\n"))
jobs = []
def add(name, m, genes, **p):
    scipy.io.mmwrite(OUT + name + ".mtx", sp.csr_matrix(m), precision=17)
    open(OUT + name + ".genes", "w").write("\n".join(genes))
    jobs.append(dict(name=name, min_cells=p.get("min_cells", "NA"), species=p.get("species", "human"),
                     optimise=p.get("optimise", True), k=p.get("k", 10), num_pcs=p.get("num_pcs", 20)))
rng = np.random.default_rng(0)
# real-data subsamples of cells
for i, n in enumerate([150, 300, 500, 750]):
    cells = np.sort(rng.choice(X.shape[1], n, replace=False))
    add(f"sub{n}", X[:, cells], G)
# random gene subsets (changes GGC / bins)
for i, frac in enumerate([0.3, 0.6]):
    genes = np.sort(rng.choice(X.shape[0], int(frac * X.shape[0]), replace=False))
    add(f"genes{int(frac*100)}", X[genes], G[genes])
# parameter variations on full data
add("k5_pcs10", X, G, k=5, num_pcs=10)
add("k20_pcs30", X, G, k=20, num_pcs=30)
add("mincells20", X, G, min_cells=20)
add("mincells150", X, G, min_cells=150)
add("mouse", X, G, species="mouse")
add("rat", X, G, species="rat")
add("noopt", X, G, optimise=False)
# synthetic: NB counts with cluster-specific marker programs, log-normalised
for s, (ncell, ngene, nclust) in enumerate([(400, 2000, 4), (1200, 3000, 8), (2500, 4000, 12)]):
    r = np.random.default_rng(100 + s)
    lab = r.integers(0, nclust, ncell)
    base = r.lognormal(-1.0, 1.2, ngene)
    mu = np.tile(base, (ncell, 1))
    for c in range(nclust):
        markers = r.choice(ngene, 40, replace=False)
        mu[np.ix_(lab == c, markers)] *= r.uniform(3, 15, 40)
    mu *= r.lognormal(0, 0.3, ncell)[:, None]
    counts = r.negative_binomial(2, 2 / (2 + mu))
    norm = np.log1p(counts / counts.sum(1, keepdims=True) * 1e4)
    names = np.array([f"G{j}" for j in range(ngene)])
    add(f"synth{ncell}x{ngene}", norm.T, names)
json.dump(jobs, open(OUT + "jobs.json", "w"))
print(len(jobs), "jobs")

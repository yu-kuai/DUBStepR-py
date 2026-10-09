import json, time, warnings, numpy as np, pandas as pd, scipy.io, scipy.sparse as sp, os, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
import dubstepr
OUT = sys.argv[1].rstrip("/") + "/"
rows = []
for j in json.load(open(OUT + "jobs.json")):
    n = j["name"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = scipy.io.mmread(OUT + n + ".mtx")
    genes = open(OUT + n + ".genes").read().split("\n")
    mc = None if j["min_cells"] == "NA" else j["min_cells"]
    t0 = time.time()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = dubstepr.dubstepr(sp.csr_matrix(m), gene_names=genes, min_cells=mc, species=j["species"],
                                    optimise_features=j["optimise"], k=j["k"], num_pcs=j["num_pcs"])
        py_err = None
    except Exception as e:
        py_err = repr(e)
    tp = time.time() - t0
    r_err = open(OUT + n + ".R.error").read().strip() if os.path.exists(OUT + n + ".R.error") else None
    if r_err or py_err:
        rows.append(dict(job=n, status=f"R error: {r_err} | Py error: {py_err}")); continue
    ci = pd.read_csv(OUT + n + ".R.corrinfo", sep="\t", keep_default_na=False)
    row = dict(job=n, ggc=f"{len(res.feature_genes)}/{len(ci)}",
               order=list(res.feature_genes) == list(ci.g),
               zdiff=float(np.max(np.abs(res.corr_info["corr.range"].to_numpy() - ci.z.to_numpy()))) if len(ci) == len(res.corr_info) else np.nan,
               elbow=f"{res.elbow_pt}/{open(OUT + n + '.R.elbow').read().strip()}", py_s=round(tp, 1))
    if j["optimise"]:
        opt = open(OUT + n + ".R.opt").read().split(); di = pd.read_csv(OUT + n + ".R.di", sep="\t")
        row.update(opt=f"{len(res.optimal_feature_genes)}/{len(opt)}", opt_same=list(res.optimal_feature_genes) == opt,
                   di_steps=len(di), di_reldiff=float(np.max(np.abs(res.density_index.to_numpy() / di.di.to_numpy() - 1))) if len(di) == len(res.density_index) else np.nan)
    rows.append(row)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 20)
print(pd.DataFrame(rows).to_string(index=False))

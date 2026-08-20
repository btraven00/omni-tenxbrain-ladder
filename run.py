#!/usr/bin/env python3
"""One rung of the TENxBrainData PCA-scaling ladder, as a DATA-contract dir.

Input is the raw 10x HDF5 (EH1039, 1.3M E18 mouse brain cells) fetched by the
omni-data stage; output is one rung:

    {name}.h5ad                    DATA contract: counts in layers/counts
    {name}.clusters_truth.tsv      (CSR, f32), no X, flat string obs columns
    {name}.clusters_truth_num.txt
    {name}_properties.yaml

Draws are NESTED: a rung of n cells is a strict subset of any larger rung at
the same seed, so a change in a downstream number is "more cells", not
"different cells". Composition is held fixed across rungs by drawing
proportionally within --group_col (the 133 sequencing libraries). Both
properties come from omni-scsampler; only the *selection* is delegated to it --
the row gather below is streaming, because the top rungs (2.6e9 nonzeros at
1.3M cells) do not fit in memory as scipy CSR, and neither scsampler's h5ad
source nor anndata's backed mode can materialize them.

There are no cell-type labels for this dataset. The truth files are emitted
header-only on purpose -- enough for the stage contract. Metric stages that
consume them must be excluded for these datasets.

Self-check (no data needed):  ./run.py --selftest
"""
import argparse
import os
import sys

import h5py
import numpy as np
import pandas as pd
from scsampler import sample

PROPS = {"batch_var": "Mouse", "sample_var": "Library", "labels_var": None}
GZIP = 4  # counts are small integers stored f32; ~3x for a few MB/s of CPU

# Streaming block limits. The source slice covering a block is read in one
# contiguous go, so SPAN -- not BLOCK -- is what bounds peak RAM: a sparse rung's
# 5000 selected cells are spread over the WHOLE source, so capping cells alone
# reads all 2.6e9 nonzeros at once (31GB decompressed, int64 indices). Capping
# the span costs one sequential pass over the source per sparse rung (~1 min).
BLOCK = 5_000        # cells
SPAN = 1 << 24       # source nonzeros: ~16.8M x 12B = ~200MB of data+indices

SD = h5py.string_dtype(encoding="utf-8")


def dict_group(group):
    """Mark a group as an anndata mapping. Empty ones need this as much as full
    ones: anndataR switch()es on encoding-type, and a missing attribute reads
    as NULL rather than "", which aborts the read with "EXPR must be a length 1
    vector". Python's anndata tolerates the omission, so it only shows up in R."""
    group.attrs["encoding-type"] = "dict"
    group.attrs["encoding-version"] = "0.1.0"
    return group


def strings(group, name, values):
    ds = group.create_dataset(name, data=list(values), dtype=SD)
    ds.attrs["encoding-type"] = "string-array"
    ds.attrs["encoding-version"] = "0.2.0"


def frame(group, columns):
    group.attrs["encoding-type"] = "dataframe"
    group.attrs["encoding-version"] = "0.2.0"
    group.attrs["_index"] = "_index"
    group.attrs.create("column-order", columns, dtype=SD)


def read_obs(f, mm, group_col):
    """Barcodes + the library/mouse split, as a positionally-indexed frame.
    10x carries the library index in the barcode's "-N" suffix; the 133
    libraries are the only batch axis this dataset has."""
    barcodes = f[f"{mm}/barcodes"].asstr()[:]
    lib = np.array([b.rsplit("-", 1)[-1] for b in barcodes])
    levels = sorted(set(lib), key=int)
    if len(levels) < 2:
        raise SystemExit(f"barcode suffixes give {len(levels)} libraries; the "
                         "-N convention does not hold for this file")
    # TENxBrainData's colData splits the libraries into the two E18 animals at
    # the midpoint of the library range; reproduced here so we don't have to
    # pull the R-only colData .rds just for one 2-level column.
    second = set(levels[len(levels) // 2:])
    return pd.DataFrame({"_index": barcodes, group_col: lib,
                         "Mouse": np.where(np.isin(lib, list(second)), "B", "A")},
                        index=np.arange(len(barcodes)))


def select(obs, n, group_col, seed):
    """Nested, composition-preserving draw. Delegated to omni-scsampler so the
    ladder shares its RNG contract (per-group keys => strict nesting)."""
    idx = sample.proportional(obs, group_col=group_col, n_target=n, seed=seed,
                              nested=True, allow_undersized=True)
    return np.sort(np.asarray(idx))


def blocks(starts, ends, span=SPAN, block=BLOCK):
    """Split the selected rows into contiguous source reads of <=span nonzeros
    and <=block cells. A single row wider than span still gets its own block."""
    b, n = 0, len(starts)
    while b < n:
        e = min(int(np.searchsorted(ends, starts[b] + span, "right")),
                b + block, n)
        yield slice(b, max(e, b + 1))
        b = max(e, b + 1)


def write_rung(f, mm, obs, idx, path, group_col):
    """Stream the selected rows into a DATA-contract h5ad.

    The 10x matrix is genes x cells in CSC, so its column pointers are already
    per-cell: CSC(genes x cells) is bit-for-bit CSR(cells x genes). The gather
    is therefore a pure copy of row segments -- no transpose, no scipy, no
    decompressed copy of the whole matrix anywhere."""
    src_indptr = f[f"{mm}/indptr"][:]
    n_genes = int(f[f"{mm}/shape"][0])
    src_data, src_indices = f[f"{mm}/data"], f[f"{mm}/indices"]

    starts, ends = src_indptr[idx], src_indptr[idx + 1]
    lengths = ends - starts
    out_indptr = np.zeros(len(idx) + 1, dtype=np.int64)
    np.cumsum(lengths, out=out_indptr[1:])
    nnz = int(out_indptr[-1])

    with h5py.File(path, "w") as g:
        g.attrs["encoding-type"] = "anndata"
        g.attrs["encoding-version"] = "0.1.0"
        counts = g.create_group("layers/counts")
        counts.attrs["encoding-type"] = "csr_matrix"
        counts.attrs["encoding-version"] = "0.1.0"
        counts.attrs.create("shape", [len(idx), n_genes], dtype=np.int64)
        dict_group(g["layers"])
        counts.create_dataset("indptr", data=out_indptr)
        data = counts.create_dataset("data", (nnz,), dtype=np.float32,
                                     chunks=(min(nnz, 1 << 20),), compression="gzip",
                                     compression_opts=GZIP)
        indices = counts.create_dataset("indices", (nnz,), dtype=np.int32,
                                        chunks=(min(nnz, 1 << 20),), compression="gzip",
                                        compression_opts=GZIP)

        for sl in blocks(starts, ends):
            lo, hi = int(starts[sl][0]), int(ends[sl][-1])
            # One contiguous read covering the block, then keep only the rows
            # we selected out of it. Every source element is touched at most
            # once per rung, sequentially.
            buf_d, buf_i = src_data[lo:hi], src_indices[lo:hi]
            take = np.concatenate([np.arange(s - lo, e - lo)
                                   for s, e in zip(starts[sl], ends[sl])])
            o0, o1 = int(out_indptr[sl.start]), int(out_indptr[sl.stop])
            data[o0:o1] = buf_d[take]
            indices[o0:o1] = buf_i[take]

        sub = obs.iloc[idx]
        o = g.create_group("obs")
        strings(o, "_index", sub["_index"])
        for c in (group_col, "Mouse"):
            strings(o, c, sub[c])
        frame(o, [group_col, "Mouse"])

        v = g.create_group("var")
        ids = f[f"{mm}/genes"].asstr()[:]
        strings(v, "_index", ids)
        strings(v, "ID", ids)
        strings(v, "Symbol", f[f"{mm}/gene_names"].asstr()[:])
        strings(v, "Type", ["Gene Expression"] * n_genes)
        frame(v, ["ID", "Symbol", "Type"])

        for empty in ("obsm", "varm", "obsp", "varp", "uns"):
            dict_group(g.create_group(empty))
    return nnz


def companions(stem, n_cells, group_col):
    # No cell-type labels exist for this dataset. Header only: valid for the
    # stage contract, useless as a label source -- which is the honest state.
    with open(f"{stem}.clusters_truth.tsv", "w") as fh:
        fh.write("cell_id\ttruths\n")
    with open(f"{stem}.clusters_truth_num.txt", "w") as fh:
        fh.write("0\n")
    props = {**PROPS, "sample_var": group_col}
    with open(f"{stem}_properties.yaml", "w") as fh:
        for k, v in props.items():
            fh.write(f"{k}: {'null' if v is None else v}\n")
        fh.write(f"n_cells: {n_cells}\n")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output_dir", required=True)
    p.add_argument("--name", required=True,
                   help="stem for every emitted file (api 0.5: the module id)")
    p.add_argument("--tenx_raw_h5", required=True,
                   help="raw 10x HDF5 from the omni-data stage")
    p.add_argument("--n_cells", type=int, required=True,
                   help="rung size; >= the source size means 'every cell'")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--group_col", default="Library",
                   help="obs column the proportional draw preserves")
    args = p.parse_args(argv)

    stem = f"{args.output_dir}/{args.name}"
    with h5py.File(args.tenx_raw_h5, "r") as f:
        mm = list(f)[0]
        obs = read_obs(f, mm, args.group_col)
        print(f"source: {len(obs)} cells x {int(f[f'{mm}/shape'][0])} genes, "
              f"{obs[args.group_col].nunique()} libraries, "
              f"{len(f[f'{mm}/data'])} nonzeros", file=sys.stderr)
        idx = select(obs, args.n_cells, args.group_col, args.seed)
        nnz = write_rung(f, mm, obs, idx, f"{stem}.h5ad.part", args.group_col)
    os.replace(f"{stem}.h5ad.part", f"{stem}.h5ad")
    companions(stem, len(idx), args.group_col)
    print(f"{args.name}: {len(idx)} cells, {nnz} nonzeros", file=sys.stderr)


def selftest():
    """Two properties the ladder rests on, neither needing the 4GB source:
    blocks() partitions the rows without ever asking for more than span at
    once (all that stands between a sparse rung and a 31GB read), and
    scsampler's draw is actually nested (a rung must be a subset of the next)."""
    rng = np.random.default_rng(0)
    src = np.zeros(20_001, dtype=np.int64)
    np.cumsum(rng.integers(1, 200, 20_000), out=src[1:])
    for n_sel in (20_000, 40, 1):  # dense draw, very sparse draw, single row
        idx = np.sort(rng.choice(20_000, n_sel, replace=False))
        starts, ends = src[idx], src[idx + 1]
        got = list(blocks(starts, ends, span=1000, block=7))
        assert [s.start for s in got] == [0] + [s.stop for s in got[:-1]]
        assert got[-1].stop == len(idx), got[-1]
        for s in got:
            assert 0 < s.stop - s.start <= 7
            assert ends[s][-1] - starts[s][0] <= 1000 or s.stop - s.start == 1
    print("blocks() ok")

    obs = pd.DataFrame({"_index": [f"c{i}" for i in range(1000)],
                        "Library": rng.integers(1, 8, 1000).astype(str)})
    prev = None
    for n in (50, 100, 200, 400):
        cur = set(select(obs, n, "Library", seed=42))
        assert prev is None or prev <= cur, f"rung {n} is not a superset"
        prev = cur
    print("nesting ok")


if __name__ == "__main__":
    selftest() if sys.argv[1:2] == ["--selftest"] else main()

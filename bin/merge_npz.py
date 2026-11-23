#!/usr/bin/env python3
"""Concatenate the region axis of several multiBigwigSummary .npz files.

score_peaks splits its reference BED into chunks and runs one
multiBigwigSummary per chunk (BED-file `-p` only parallelises by chromosome,
so a single process otherwise ends up doing all of chr1 while the rest idle).
Each chunk produces a deeptools .npz holding

    matrix  : ndarray, (regions_in_chunk, n_samples)
    labels  : ndarray of sample names

Stacking the matrices in chunk order reproduces exactly the array a single
unchunked run would have written (multiBigwigSummary keeps BED input order),
so plotPCA / plotCorrelation read the result unchanged.

    merge_npz.py '<glob for the part .npz files>' <output .npz>
"""
import glob
import sys

import numpy as np

pattern, out = sys.argv[1], sys.argv[2]
parts = sorted(glob.glob(pattern))
if not parts:
    sys.exit(f"merge_npz.py: no files match {pattern!r}")

matrices = []
labels = None
for f in parts:
    d = np.load(f, allow_pickle=True)
    matrices.append(d["matrix"])
    if labels is None:
        labels = d["labels"]

np.savez_compressed(out, matrix=np.vstack(matrices), labels=labels)

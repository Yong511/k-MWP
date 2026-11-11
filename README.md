# k-MWP

k-MWP is a C++ implementation for finding a maximum-weight $k$-plex in a vertex-weighted graph. It maximizes the sum of the selected vertex weights subject to the $k$-plex constraint and the minimum size $|S| \ge 2k-1$.

## Build and run

Use GNU Make and a C++17-compatible GCC toolchain on Linux. Run these commands from the repository root:

```bash
cd src
make -j2
./MWPlex -g 620 -k 3 -a MWPA -w 0
```

Run `MWPlex` from `src/`: it reads input from `../dataset/<id>/`. Set `-g` to `500`, `620`, or `889`, and use `-k` to choose the plex parameter. The solver sets the minimum size to `2*k-1` automatically. `-a MWPA` selects the full solver; `-a MWPA-` and `-a Base` select the two internal variants. `-w 0` loads the bundled `weight.txt`.

The final output includes `k-MWP Size` and `Total Weight`. The legacy `kplexes.txt` output can use internal vertex IDs after graph reduction; do not treat those IDs as original graph labels.

## Included datasets

The repository contains the three smallest distinct benchmark graphs, using the same preprocessed inputs as the experiment runs. A separate toy graph and duplicate benchmark folders are not included.

| Folder | Dataset | Vertices | Undirected edges | Expected weight at `-k 3 -a MWPA -w 0` |
| --- | --- | ---: | ---: | ---: |
| `500` | [US Airport Network](https://toreopsahl.com/datasets/#usairport) | 500 | 2,980 | 542,443,079 |
| `620` | [fb-pages-food](https://networkrepository.com/fb-pages-food.php) | 620 | 2,091 | 1,416 |
| `889` | [wiki-Vote](https://networkrepository.com/soc-wiki-Vote.php) | 889 | 2,914 | 1,307 |

Each folder contains `b_degree.bin` (the integer-width marker, vertex count, adjacency-entry count, and degree array), `b_adj.bin` (the adjacency lists), and `weight.txt` (one vertex weight per line, in zero-based vertex order). The binary adjacency-entry count includes both directions of each undirected edge. The included `edges.txt` and `b_weight.bin` are auxiliary files; the command above uses `weight.txt` rather than `b_weight.bin`.

Only the solver source, build files, and these three datasets are included. The external baselines, full benchmark suite, and experiment logs are not part of this repository.

## Third-party notices

The bundled command-line parser [`popl`](https://github.com/badaix/popl) is distributed under its [MIT license](LICENSES/popl-MIT.txt). Dataset links above identify the upstream sources; the parser license does not apply to those datasets or to the rest of this repository.

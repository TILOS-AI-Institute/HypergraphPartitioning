# Hint partitions

K-SpecPart is a **refinement** tool: it improves an existing partition ("hint")
rather than partitioning from scratch. This folder holds the hint partitions
used by the benchmark sweep ([`../titan_sweep.jl`](../titan_sweep.jl)) and is a
convenient place to keep your own hints.

## Layout

```
hints/
└── <suite>/                        e.g. Titan23
    └── ub_factor_<imb>/            e.g. ub_factor_2  (imbalance used to make the hint)
        ├── 2_way/
        │   └── <design>.hgr.specpart.ubfactor.<imb>.part.2
        ├── 3_way/
        │   └── <design>.hgr.specpart.ubfactor.<imb>.part.3
        └── 4_way/
            └── <design>.hgr.specpart.ubfactor.<imb>.part.4
```

Each hint file is an hMETIS-style partition: **one 0-indexed block id per line**,
in vertex order (the format `titan_sweep.jl` and `run_kspecpart.jl` expect via
`--hint`).

## Naming convention

`titan_sweep.jl` locates a hint as:

```
<HINT_DIR>/<k>_way/<design>.hgr.specpart.ubfactor.<IMB>.part.<k>
```

where `HINT_DIR` defaults to `hints/Titan23/ub_factor_2` (override with
`SWEEP_HINT_DIR`), `IMB` is the sweep imbalance (default `2`), `design` is the
benchmark name (matching `<design>.hgr`), and `k` is the block count.

Example: `hints/Titan23/ub_factor_2/4_way/des90.hgr.specpart.ubfactor.2.part.4`

## Notes

- Hint files are **tracked in git** (an explicit exception in the top-level
  `.gitignore`, which otherwise ignores `*.part.*` run artifacts).
- The `.gitkeep` files just keep the empty `*_way/` directories in git until the
  hints are uploaded.

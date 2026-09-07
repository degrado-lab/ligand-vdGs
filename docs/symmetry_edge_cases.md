# Symmetry-determination edge cases in `generate_vdgs`

---

## Summary of edge-case classes

| # | Class | Representative fragments | Ligands affected | Severity |
|---|---|---|---|---|
| 1 | `_has_deliberate_charge_assignment` carve-out's true branch never fires in production | `[O-]S([O-])([O-])O` (3 ligands) | 0 | Low (unreachable in practice, not unreferenced) |
| 2 | 24-permutation CGs multiply every pairwise RMSD by 24 | `O=P(O)(O)O`, `O=S(=O)(O)O`, `C[N+](C)(C)C` | 3.8k | Medium (cost) |
| 3 | Ring-conformer / cis-trans cases | `CC1OC1C`, `CN1CCC1` | ~125 | Low |

Settled classes are no longer listed here: ring positions merged with exocyclic
ones (`_automorphism_graph` carries no ring-membership label) is in
`docs/pitfalls.md` under Known issues, and the per-fragment symmetry-override
inconsistency was fixed (every group now comes from the derived rules).

---

## 1. Carve-out whose true branch is unreachable in production

`_has_deliberate_charge_assignment` (four or more single-bonded terminal atoms of
one element with non-uniform charges) fires on exactly one fragment in the whole
dictionary — `[O-]S([O-])([O-])O` — which has **3 ligands**. The function itself
is not dead code: it is called from `_find_resonance_terminal_groups` on every
fragment's automorphism computation. But its `True` branch never fires for a
fragment that actually reaches a production library, because
`select_fragments()` (`extract_fragment_smiles.py` /
`make_sge_scripts_for_frags.py`) filters the fragment dict to the
`--min-instances` (default 250 estimated CG occurrences in the parent database) *before* any per-fragment vdG
generation runs, and hit finding reads the stored group from `cg_symmetry.npz`
rather than re-deriving it. So at default settings, no production run has ever
exercised the behavior its docstring example describes, even though the code
path executes on every call.

---

## 2. 24-permutation CGs

| Fragment | perms | ligands |
|---|---|---|
| `O=P(O)(O)O` | 24 | 3420 |
| `C[N+](C)(C)C` | 24 | 222 |
| `O=S(=O)(O)O` | 24 | 166 |

Every stage-1 pairwise RMSD for these buckets is 24 Kabsch fits instead of 1.
Phosphate is the fragment the combinatorial-explosion note in `CLAUDE.md` refers
to. All 24 are needed, so speed has to come from the prefilters -- the
CA-CG-COM fingerprint and the backbone lower bound.

---

## 3. Ring / cis-trans cases

| Fragment | perms | ligands | parent-asym% | note |
|---|---|---|---|---|
| `CC1OC1C` | 2 | 55 | **85%** | 2,3-dimethyloxirane: the swap merges cis and trans, and is improper |
| `CN1CCC1` | 2 | 69 | 33% | |
| `C1CCN1` | 2 | 76 | 7% | azetidine, fine |
| `C1COC1` | 2 | 83 | 6% | oxetane, fine |
| `NC1CC1`, `CC1(N)CC1`, `CNC1CC1`, `NCC1CC1` | 2 | 268/83/122/84 | 0-1% | fine |

---

The two standing prohibitions that apply to every class above — do not drop
"improper" permutations, do not replace a consumer's enumeration with a reliance
on the stored CG atom order — are in `docs/pitfalls.md` under Caveats, "CG atom
correspondence comes from the automorphism group, not from stored order".

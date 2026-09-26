"""Cheap necessary-condition filter for SMARTS-vs-ligand matching: skip fragments a
ligand's atom counts fall short of, without calling OpenBabel's matcher. Necessary,
not sufficient; unclassifiable SMARTS get an all-zero row."""
import re
import numpy as np

# Element symbols in the order cg.ring_size_constraints scans them: two-letter
# first, so 'Cl' is not read as 'C' followed by a ring-closure digit.
_ORGANIC_SUBSET = ('Cl', 'Br', 'B', 'C', 'N', 'O', 'P', 'S', 'F', 'I',
                   'b', 'c', 'n', 'o', 'p', 's')

# Atom classes counted; every key atom and every OBMol atom lands in exactly
# one, making a per-class count a necessary condition for a match.
_CLS_AROMATIC, _CLS_ACYCLIC, _CLS_RING, _CLS_ANY = 0, 1, 2, 3
_N_CLS = 4

def parse_pattern_atoms(smarts):
    """[(element, is_aromatic, primitives)] in pattern-atom order, or None if the
    key uses grammar this module doesn't parse (disables the prefilter for that
    pattern rather than risking a wrong skip)."""
    atoms, i, n = [], 0, len(smarts)
    while i < n:
        char = smarts[i]
        if char == '[':
            depth, j = 1, i + 1
            while j < n and depth:
                if smarts[j] == '[':
                    depth += 1
                elif smarts[j] == ']':
                    depth -= 1
                j += 1
            primitives = smarts[i + 1:j - 1]
            if ',' in primitives:  # e.g. `[c,n]`: no single element required.
                return None
            match = re.match(r'([A-Z][a-z]?|[a-z])', primitives)
            if not match:
                return None
            symbol = match.group(1)
            atoms.append((symbol.capitalize(), symbol.islower(), primitives))
            i = j
            continue
        if char == '*':
            return None
        for symbol in _ORGANIC_SUBSET:
            if smarts.startswith(symbol, i):
                atoms.append((symbol.capitalize(), symbol.islower(), ''))
                i += len(symbol)
                break
        else:
            # Bond, branch, ring-closure digit, %nn, or dot: not an atom.
            i += 2 if char == '%' else 1
    return atoms

def build_prefilter(fragments, ring_constraints):
    """(requirement matrix, element column index, ring-size column index).

    Row i holds the atom count fragment i needs per (element, class) and per
    ring size; `prefilter_candidates` skips any fragment a ligand falls short
    of, without invoking OpenBabel.
    """
    parsed = [parse_pattern_atoms(f) for f in fragments]
    symbols = sorted({e for atoms in parsed if atoms for e, _, _ in atoms})
    sym_col = {s: i * _N_CLS for i, s in enumerate(symbols)}
    ring_offset = _N_CLS * len(symbols)
    sizes = sorted({size for constraints in ring_constraints
                    for size in constraints if size})
    ring_col = {size: ring_offset + i for i, size in enumerate(sizes)}
    req = np.zeros((len(fragments), ring_offset + len(sizes)), dtype=np.int16)
    for row, (atoms, constraints) in enumerate(zip(parsed, ring_constraints)):
        if atoms is None or len(atoms) != len(constraints):
            continue
        for (element, aromatic, primitives), ring_size in zip(atoms, constraints):
            base = sym_col[element]
            req[row, base + _CLS_ANY] += 1
            if aromatic:
                req[row, base + _CLS_AROMATIC] += 1
            elif ring_size is not None:
                req[row, base + _CLS_RING] += 1
                req[row, ring_col[ring_size]] += 1
            elif '!R' in primitives:
                req[row, base + _CLS_ACYCLIC] += 1
            # Otherwise the atom constrains only its element, already counted.
    return req, sym_col, ring_col

def prefilter_candidates(mol, prefilter):
    """Indices of fragments whose atom counts this molecule could supply."""
    from openbabel import openbabel as ob
    req, sym_col, ring_col = prefilter
    counts = np.zeros(req.shape[1], dtype=np.int16)
    for atom in ob.OBMolAtomIter(mol):
        base = sym_col.get(ob.GetSymbol(atom.GetAtomicNum()))
        if base is None:
            continue  # No fragment needs this element.
        counts[base + _CLS_ANY] += 1
        if atom.IsAromatic():
            counts[base + _CLS_AROMATIC] += 1
        elif not atom.IsInRing():
            counts[base + _CLS_ACYCLIC] += 1
        else:
            counts[base + _CLS_RING] += 1
            column = ring_col.get(atom.MemberOfRingSize())
            if column is not None:
                counts[column] += 1
    return np.flatnonzero((req <= counts).all(axis=1))

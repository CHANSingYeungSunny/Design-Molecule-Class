"""Representation adapters: reference molecule -> the exact bytes the model reads.

Operational definitions (protocol section 5)
-------------------------------------------
smiles    RDKit canonical SMILES, stereo-free pool. For Identity positives a
          *rooted* non-canonical traversal is used, so the two strings differ
          while describing the same molecule; the pair is never a character
          comparison of one string with itself.
selfies   SELFIES encoded from the same reference structure, with a decode
          round-trip check back to the identical canonical SMILES.
our_code  The professor prototype's own constructor text (``to_code()``).
          Identity positives use the same graph with the atom/bond *records*
          reordered -- the model is given the code, never a SMILES.

Not in round 1: img, graph2d, graph3d, moljson. The registry marks them
``implemented: false`` with the reason; they are listed so the comparison table
has a place to live, but no stub pretends to produce them.

None of these functions may leak the answer: no names, no reference SMILES
alongside the image/code, no standard values, no internal example id.
"""
from __future__ import annotations

import dataclasses

import selfies as sf
from rdkit import Chem

from .core_loader import load_core


class RepresentationError(RuntimeError):
    """The reference structure cannot be expressed in this representation."""


# --------------------------------------------------------------------------
# reference structure helpers
# --------------------------------------------------------------------------

def reference_mol(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RepresentationError('RDKit cannot parse %r' % smiles)
    return mol


def canonical_smiles(mol: Chem.Mol) -> str:
    """Canonical, stereo-free reference SMILES.

    Explicit hydrogens are removed first. A few ChEBI-20 rows are written with
    every hydrogen as a real ``[H]`` vertex; without RemoveHs, RDKit's canonical
    output for those is not stable under a re-parse, so the same molecule would
    compare unequal to itself. Removing H here normalises both spellings to the
    identical string without changing the molecule (the formula is unchanged).
    """
    return Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=False, canonical=True)


def rooted_smiles(mol: Chem.Mol, root: int) -> str:
    """A valid, deliberately non-canonical traversal starting at ``root``."""
    return Chem.MolToSmiles(mol, isomericSmiles=False, canonical=False,
                            rootedAtAtom=root)


# --------------------------------------------------------------------------
# smiles
# --------------------------------------------------------------------------

def smiles_text(smiles: str) -> str:
    mol = reference_mol(smiles)
    out = canonical_smiles(mol)
    # the pool filter requires Our Code to read the canonical form too, so the
    # same molecule list is used for every representation
    load_core().Molecule(out)
    return out


def smiles_identity_variant(smiles: str, variant: int) -> str:
    mol = reference_mol(smiles)
    if variant == 0:
        return canonical_smiles(mol)
    atoms = list(mol.GetAtoms())
    root = (variant * 7 + 1) % len(atoms)
    out = rooted_smiles(mol, root)
    if out == canonical_smiles(mol) and len(atoms) > 1:
        # extremely rare for a small ring-only molecule; pick the last atom
        out = rooted_smiles(mol, len(atoms) - 1)
    return out


# --------------------------------------------------------------------------
# selfies
# --------------------------------------------------------------------------

def selfies_text(smiles: str) -> str:
    mol = reference_mol(smiles)
    canonical = canonical_smiles(mol)
    encoded = sf.encoder(canonical)
    if encoded is None:
        raise RepresentationError('SELFIES cannot encode %r' % canonical)
    # decode round-trip check: the encoding must come back as the same molecule
    back = sf.decoder(encoded)
    if back is None or canonical_smiles(reference_mol(back)) != canonical:
        raise RepresentationError(
            'SELFIES round-trip changed the molecule for %r' % canonical)
    return encoded


def selfies_identity_variant(smiles: str, variant: int) -> str:
    return sf.encoder(smiles_identity_variant(smiles, variant))


# --------------------------------------------------------------------------
# our code
# --------------------------------------------------------------------------

def _core_molecule(smiles: str):
    core = load_core()
    return core.Molecule(smiles)


def our_code_text(smiles: str) -> str:
    """The professor prototype's constructor text for this molecule."""
    try:
        return _core_molecule(smiles).to_code()
    except Exception as exc:  # noqa: BLE001 - surfaced as a conversion failure
        raise RepresentationError('Our Code cannot express %r: %s' % (smiles, exc)) from exc


def our_code_identity_variant(smiles: str, variant: int) -> str:
    """Same graph, records reordered -- a different legal arrangement.

    The atom ids are preserved; only the order of the (atoms, bonds) tuples
    changes. validate_graph() accepts this, and the resulting text differs
    from variant 0 while describing the identical molecule.
    """
    mol = _core_molecule(smiles)
    if variant == 0:
        return mol.to_code()
    atoms = list(mol.atoms)
    bonds = list(mol.bonds)
    if len(atoms) > 1:
        shift = (variant * 5 + 1) % len(atoms)
        atoms = atoms[shift:] + atoms[:shift]
    if len(bonds) > 1:
        shift = (variant * 3 + 1) % len(bonds)
        bonds = bonds[shift:] + bonds[:shift]
    reordered = dataclasses.replace(mol, atoms=tuple(atoms), bonds=tuple(bonds))
    return reordered.to_code()


# --------------------------------------------------------------------------
# registry
# --------------------------------------------------------------------------

@dataclasses.dataclass(frozen=True)
class Representation:
    name: str
    kind: str                      # 'text' or 'image'
    implemented: bool
    render: object                 # callable(smiles) -> payload
    render_variant: object         # callable(smiles, variant) -> payload
    note: str = ''


REPRESENTATIONS = {
    'smiles': Representation(
        'smiles', 'text', True, smiles_text, smiles_identity_variant,
        'RDKit canonical SMILES; Identity positives use a rooted non-canonical traversal.'),
    'selfies': Representation(
        'selfies', 'text', True, selfies_text, selfies_identity_variant,
        'SELFIES from the same reference structure, with a decode round-trip check.'),
    'our_code': Representation(
        'our_code', 'text', True, our_code_text, our_code_identity_variant,
        "The professor prototype's Molecule(...) constructor text; Identity positives reorder records."),
    }

# The five ADDED representations (img, graph2d, graph3d, cml, moljson) are
# deliberately NOT listed above. They have real renderers in ``reps_ext``, and
# keeping a second, stale copy here is what caused two separate defects: a
# report layer that could not find 'cml' at all, and a `implemented: False`
# placeholder that made four representations with hundreds of real responses
# print as "model does not support this" / "not implemented". ``render`` below
# is the single entry point and reaches both registries, so there is one place
# to look and no second source of truth to drift.
ROUND1_REPRESENTATIONS = ('our_code', 'smiles', 'selfies')
ADDED_REPRESENTATIONS = ('img', 'graph2d', 'graph3d', 'cml', 'moljson')
ALL_REPRESENTATION_NAMES = ('smiles', 'selfies', 'our_code') + ADDED_REPRESENTATIONS


def render(name: str, smiles: str, variant: int = 0):
    """Render a representation payload.

    The five ADDED representations live in reps_ext and are reached from here,
    so every caller keeps using one entry point and the item builder needs no
    special case.
    """
    if name in REPRESENTATIONS and REPRESENTATIONS[name].implemented:
        rep = REPRESENTATIONS[name]
        return rep.render_variant(smiles, variant) if variant else rep.render(smiles)
    from . import reps_ext
    if name in reps_ext.RENDERERS:
        return reps_ext.render(name, smiles, variant)
    rep = REPRESENTATIONS.get(name)
    raise RepresentationError('%s is not implemented: %s'
                              % (name, rep.note if rep else 'unknown representation'))
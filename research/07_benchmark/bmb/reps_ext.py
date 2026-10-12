"""The five ADDED representations: img, Graph2D, Graph3D, CML, MolJSON.

Operational definitions follow the extension brief. Two things are stated here
once so no report has to re-derive them:

  * Graph2D / Graph3D are MOL V2000 CONNECTION TABLES AS TEXT, not pictures and
    not graph neural networks. "Graph" in the task list means the atom/bond
    relation, and this project serialises it as a MOL block.
  * The two differ only in where the coordinates come from: Graph2D uses the
    2D layout (z = 0), Graph3D uses a computed conformer. Graph3D is therefore
    an EXTRA-INFORMATION condition, not a different task.

Nothing here ever emits a name, SMILES, database id, descriptor, functional
group label or answer. Every function is pure: same input -> same bytes.
"""
from __future__ import annotations

import hashlib
import io
import json
import pathlib
import sys
import xml.etree.ElementTree as ET

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

RDLogger.DisableLog('rdApp.*')

MOLJSON_SRC = pathlib.Path(__file__).resolve().parents[3] / 'reference' / 'MolJSON' / 'src'
if str(MOLJSON_SRC) not in sys.path:
    sys.path.insert(0, str(MOLJSON_SRC))

# --------------------------------------------------------------------------
# frozen conversion parameters (fixed before any new main request)
# --------------------------------------------------------------------------
CONVERSION = {
    'version': 'ext-v2',
    'img': {
        'size_px': 1024,
        'background': 'white',
        'drawer': 'rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DCairo',
        'metadata': 'stripped by re-encoding the PNG so only pixel data remains',
        'no_title_no_atom_numbers_no_highlights': True,
        'identity_views': 'two renders sharing the rule; neither depends on the label',
    },
    'graph2d': {
        'format': 'MOL V2000',
        'coords': 'RDKit 2D layout (AllChem.Compute2DCoords), z = 0',
        'title_line': 'blank',
        'kekulize': True,
        'roundtrip_check': 'Chem.MolFromMolBlock',
    },
    'graph3d': {
        'format': 'MOL V2000',
        'embed': 'AllChem.EmbedMultipleConfs, ETKDGv3, useSmallRingTorsions=True',
        'seed': 20261005,
        'n_candidates': 10,
        'optimize': 'MMFFOptimizeMoleculeConfs, MMFF94s',
        'max_iters': 500,
        'max_iters_extended': 2000,
        'fallback': 'UFF if MMFF has no parameters for the molecule, recorded per molecule',
        'select': 'lowest energy among candidates that CONVERGED (converged flag 0)',
        'on_failure': ('no candidate converged even at max_iters_extended, or embedding '
                       'produced nothing under both the declared rule and the random-coordinate '
                       'fallback, or neither MMFF nor UFF covers the molecule -> the molecule is '
                       'recorded as a conversion failure. It is never replaced, and no score of '
                       '0 is written for it.'),
        'fallback_order': ('v1 rule first and unchanged; random-coordinate embedding only if '
                           'nothing embedded; the extended iteration cap only if nothing '
                           'converged. A molecule the v1 rule already expressed never reaches a '
                           'fallback, so frozen v1 payloads are unaffected.'),
        'keep_hydrogens': False,
        'stereo': ('includeStereo=False, so the chiral flag stays 0 and no wedge/hash bond '
                   'is written: the export cannot add a stereochemical label the source '
                   'structure did not carry. Reading a 3D block back makes RDKit PERCEIVE '
                   'stereocentres, which is why the round-trip check compares the '
                   'stereo-free constitution.'),
        '_project_choice': ('candidate count, seed, iteration cap, optimiser and fallback '
                            'order are THIS project\'s choices, not a published standard'),
    },
    'cml': {
        'schema': 'CML atomArray / bondArray',
        'coords': 'none (structure only, so it stays a text-format comparison)',
        'aromaticity': 'explicit Kekule bond orders plus an aromatic flag',
    },
    'moljson': {
        'source': 'github.com/oxpig/MolJSON',
        'commit': '047b135e628b091ec57a65d6da06989554f3a229',
        'license': 'BSD-3-Clause',
        'schema_file': 'reference/MolJSON/schemas/moljson.schema.json',
        'required_keys': ['atoms', 'bonds', 'charges', 'aromatic_n_h'],
        'atom_id_style': 'a',
        '_atom_id_style_reason': (
            "MolToJSON's default style names an atom after its element and a running "
            "counter (C1, C2, ... O4). That spells out how many atoms of each element "
            "the molecule has, which IS the answer to the atom-count task, so it would "
            "hand one representation a free win. The sequential style (a1, a2, ...) is "
            "the same converter, the same schema and the same document shape, with ids "
            "that carry no element count. This is a fairness adjustment, not a change "
            "of schema."),
        'include_empty_fields': True,
        '_include_empty_fields_reason': (
            "the schema lists charges and aromatic_n_h as required; the converter's "
            "default omits them when empty, so the payload would fail its own schema. "
            "Emitting them as null is the schema's declared representation of 'none'."),
    },
}

IDENTITY_VIEWS = (0, 1)

# Fixed, label-independent view angles. Slot A always uses view 0 and slot B
# always uses view 1, for a positive pair and a negative pair alike; nothing
# here is chosen from the answer.
VIEW_2D_DEGREES = 37.0
VIEW_3D_EULER_DEGREES = (37.0, -23.0, 11.0)


class ConversionError(RuntimeError):
    """This molecule cannot be expressed in this representation."""


# --------------------------------------------------------------------------
# memoisation
# --------------------------------------------------------------------------
# The renderers are pure, but a conformer search is expensive and the same
# molecule is rendered once per task that draws it. Caching the RESULT keeps
# every payload byte-identical to an uncached run while removing the repeats.
# Nothing here is keyed on the answer, and the cache never survives a process,
# so a frozen payload is still the one written to disk by extension-prep.
_CACHE = {}
CACHE_STATS = {'memory_hits': 0, 'disk_hits': 0, 'misses': 0}


def _cached(key, build):
    if key not in _CACHE:
        _CACHE[key] = build()
    return _CACHE[key]


# A conformer search is the one expensive step, and it is also the one the
# protocol requires to be FROZEN: the same molecule must serialise to the same
# coordinates on every run, and a coordinate must never be regenerated between
# the freeze and the request. Keeping the finished blocks on disk turns that
# requirement into something the file system enforces rather than something the
# code merely intends. The key carries the conversion version and every
# parameter that shapes the geometry, so changing any of them invalidates the
# entry instead of silently reusing an old conformer.
CONFORMER_CACHE_DIR = pathlib.Path(__file__).resolve().parents[1] / 'data' / 'prepared' / 'ext_cache'


def _conformer_cache_key(smiles):
    spec = CONVERSION['graph3d']
    # Explicit values, not the prose descriptions: clarifying a docstring must
    # not invalidate a cache, but changing any number that shapes the geometry
    # must. Editing the text fields would otherwise silently be treated as a
    # parameter change, and a real parameter change would be missed if only the
    # prose were keyed on.
    material = '|'.join(str(x) for x in [
        CONVERSION['version'], smiles, spec['format'], spec['seed'],
        spec['n_candidates'], spec['max_iters'], spec['max_iters_extended'],
        spec['optimize'], spec['fallback'],
        spec['keep_hydrogens'], VIEW_3D_EULER_DEGREES,
    ])
    return hashlib.sha256(material.encode('utf-8')).hexdigest()


def cache_entry(smiles):
    path = CONFORMER_CACHE_DIR / ('%s.json' % _conformer_cache_key(smiles))
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding='utf-8'))
    except Exception:  # noqa: BLE001 - a corrupt entry is a miss, not a crash
        return None


def store_entry(smiles, entry):
    CONFORMER_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    path = CONFORMER_CACHE_DIR / ('%s.json' % _conformer_cache_key(smiles))
    path.write_text(json.dumps(entry, ensure_ascii=False), encoding='utf-8')


def _mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ConversionError('RDKit cannot parse %r' % smiles)
    return mol


def _clean_mol(smiles, add_hs=False):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ConversionError('RDKit cannot parse %r' % smiles)
    if add_hs:
        mol = Chem.AddHs(mol)
    return mol


# --------------------------------------------------------------------------
# img
# --------------------------------------------------------------------------

def img_png(smiles, variant=0):
    """A 2D structure PNG with no title, no atom numbers and no metadata."""
    return _cached(('img', smiles, variant), lambda: _img_png_uncached(smiles, variant))


def _img_png_uncached(smiles, variant):
    mol = _mol(smiles)
    # a deterministic layout; variant rotates the depiction by re-running the
    # coordinate generation, which is label-independent
    AllChem.Compute2DCoords(mol)
    if variant:
        conf = mol.GetConformer()
        angle = 37.0 * variant
        import math
        cos_a, sin_a = math.cos(math.radians(angle)), math.sin(math.radians(angle))
        for i in range(mol.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (p.x * cos_a - p.y * sin_a,
                                     p.x * sin_a + p.y * cos_a, 0.0))
    size = CONVERSION['img']['size_px']
    drawer = rdMolDraw2D.MolDraw2DCairo(size, size)
    options = drawer.drawOptions()
    options.addStereoAnnotation = False
    options.addAtomIndices = False
    options.includeMetadata = False      # RDKit would otherwise embed the molecule
    options.clearBackground = True
    options.padding = 0.08
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    raw = drawer.GetDrawingText()
    # re-encode through Pillow: this drops every PNG text chunk RDKit may have
    # written, so only pixels survive
    from PIL import Image
    with Image.open(io.BytesIO(raw)) as image:
        buffer = io.BytesIO()
        image.convert('RGB').save(buffer, format='PNG', optimize=False)
        return buffer.getvalue()


def img_identity_variant(smiles, variant):
    return img_png(smiles, variant=1 if variant else 0)


def png_has_text_metadata(blob):
    """True if the PNG carries any tEXt/iTXt/zTXt chunk (metadata)."""
    import struct
    if blob[:8] != b'\x89PNG\r\n\x1a\n':
        return True
    pos, end = 8, len(blob)
    while pos + 8 <= end:
        length = struct.unpack('>I', blob[pos:pos + 4])[0]
        kind = blob[pos + 4:pos + 8]
        if kind in (b'tEXt', b'iTXt', b'zTXt', b'eXIf'):
            return True
        if kind == b'IEND':
            break
        pos += 12 + length
    return False


# --------------------------------------------------------------------------
# graph2d
# --------------------------------------------------------------------------

def _mol_block(mol, conf_id=-1, three_d=False):
    mol = Chem.Mol(mol)
    mol.SetProp('_Name', '')          # blank title: no name, id or SMILES
    # includeStereo=False keeps the chiral flag at 0 and writes no wedge/hash
    # bond, so the export cannot ADD a stereochemical label the reference
    # structure did not carry. The reference pool is stereo-free by construction.
    return Chem.MolToMolBlock(mol, confId=conf_id, kekulize=True,
                              includeStereo=False)


def _constitution(smiles_or_mol):
    """Canonical SMILES with stereochemistry OFF.

    A MOL block written from 3D coordinates makes RDKit PERCEIVE stereocentres
    when it is read back, so the isomeric SMILES gains @ and / marks even though
    not one atom or bond changed. Comparing with stereochemistry ON therefore
    reports a difference for every molecule that happens to contain a
    stereocentre -- which is a defect in the check, not in the conversion.
    Constitution (atoms, bonds, bond orders, charges) is what must be preserved.
    """
    mol = (smiles_or_mol if isinstance(smiles_or_mol, Chem.Mol)
           else Chem.MolFromSmiles(smiles_or_mol))
    return Chem.MolToSmiles(mol, isomericSmiles=False)


def _block_has_no_stereo_markers(block):
    """True when the MOL block itself carries no stereo label.

    Field 12-14 of the counts line is the chiral flag; the bond block's field
    4 is the bond stereo care box (1 = wedge up, 6 = wedge down). Both must be
    clear, otherwise the export would be adding a stereo claim of its own.
    """
    lines = block.splitlines()
    if len(lines) < 4:
        return False
    counts = lines[3]
    if counts[12:15].strip() not in ('0', ''):
        return False
    try:
        n_atoms = int(counts[0:3])
        n_bonds = int(counts[3:6])
    except ValueError:
        return False
    if len(lines) < 4 + n_atoms + n_bonds:
        return False
    # ONLY the bond block. Reading past it would run into the property lines,
    # where "M  CHG  1   4  -1" has a 1 in the same column as a wedge flag.
    for line in lines[4 + n_atoms: 4 + n_atoms + n_bonds]:
        if len(line) >= 12 and line[9:12].strip() in ('1', '6'):
            return False
    return True


def _verify_block(block, smiles, label):
    back = Chem.MolFromMolBlock(block, sanitize=True, removeHs=False)
    if back is None:
        raise ConversionError('%s MOL block does not round-trip for %r' % (label, smiles))
    if _constitution(back) != _constitution(smiles):
        raise ConversionError('%s round-trip changed the constitution for %r'
                              % (label, smiles))
    if not _block_has_no_stereo_markers(block):
        raise ConversionError('%s added a stereochemical marker for %r' % (label, smiles))
    return back


def _rotate_2d(conformer, degrees):
    import math
    cos_a, sin_a = math.cos(math.radians(degrees)), math.sin(math.radians(degrees))
    for i in range(conformer.GetNumAtoms()):
        p = conformer.GetAtomPosition(i)
        conformer.SetAtomPosition(i, (p.x * cos_a - p.y * sin_a,
                                      p.x * sin_a + p.y * cos_a, 0.0))


def _renumber_order(n, variant):
    """A fixed atom permutation. Same rule for every molecule and every label."""
    if n <= 1:
        return list(range(n))
    shift = (n // 3) or 1
    return list(range(shift, n)) + list(range(shift))


def graph2d_text(smiles, variant=0):
    def build():
        mol = _clean_mol(smiles)
        AllChem.Compute2DCoords(mol)
        if variant:
            mol = Chem.RenumberAtoms(mol, _renumber_order(mol.GetNumAtoms(), variant))
            _rotate_2d(mol.GetConformer(0), VIEW_2D_DEGREES)
        block = _mol_block(mol, conf_id=0)
        _verify_block(block, smiles, 'Graph2D')
        return block

    return _cached(('graph2d', smiles, variant), build)


def graph2d_identity_variant(smiles, variant):
    return graph2d_text(smiles, variant)


# --------------------------------------------------------------------------
# graph3d
# --------------------------------------------------------------------------

def _embed(mol3d_heavy, use_random_coords=False, n_conformers=None):
    params = AllChem.ETKDGv3()
    params.randomSeed = CONVERSION["graph3d"]["seed"]
    params.useSmallRingTorsions = True
    if use_random_coords:
        params.useRandomCoords = True
    return list(AllChem.EmbedMultipleConfs(
        mol3d_heavy,
        numConfs=n_conformers or CONVERSION['graph3d']['n_candidates'],
        params=params))


def _optimise(mol3d_heavy, method, max_iters):
    if method == 'MMFF94s':
        return AllChem.MMFFOptimizeMoleculeConfs(mol3d_heavy, maxIters=max_iters,
                                                 mmffVariant='MMFF94s')
    return AllChem.UFFOptimizeMoleculeConfs(mol3d_heavy, maxIters=max_iters)


def _best_converged(ids, results):
    best, best_energy = None, None
    for conf_id, (converged, energy) in zip(ids, results):
        if converged != 0:
            continue
        if best_energy is None or energy < best_energy:
            best, best_energy = conf_id, energy
    return best, best_energy


def _conformer(mol3d_heavy, smiles):
    """Embed, optimise and pick the lowest-energy converged conformer.

    FALLBACK DISCIPLINE (revision ext-v2)
    -------------------------------------
    The v1 rule is tried FIRST and unchanged: ETKDGv3 with a fixed seed,
    ``useSmallRingTorsions``, N candidates; MMFF94s (UFF if MMFF has no
    parameters) optimised to convergence within the declared iteration cap; the
    lowest-energy CONVERGED candidate wins.

    Two fallbacks are reached ONLY when the v1 rule produces no usable
    conformer at all. That ordering is the whole point: a molecule that v1 could
    already express never reaches a fallback, so its payload cannot move, and
    the frozen v1 results stay reproducible. A fallback changes coverage, never
    the geometry of something already covered.

      * embedding produced nothing -> retry ETKDGv3 with ``useRandomCoords``,
        same seed, same candidate count. This is a documented RDKit option for
        difficult embeddings, not a relaxation of the chemistry.
      * nothing converged at the declared cap -> re-optimise the SAME
        candidates with a higher cap. A conformer that is still unconverged is
        never accepted, and an unoptimised geometry is never written.

    Neither fallback touches bonds, charges, aromaticity or the metal atoms: a
    molecule whose force field has no parameters is still a hard failure.
    ``path`` in the returned metadata records which route was taken, so a
    reader can tell a v1 payload from a recovered one.
    """
    spec = CONVERSION['graph3d']
    path = []
    ids = _embed(mol3d_heavy)
    if ids:
        path.append('etkdgv3')
    else:
        ids = _embed(mol3d_heavy, use_random_coords=True)
        if not ids:
            raise ConversionError('embedding failed for %r under both the declared rule '
                                  'and the random-coordinate fallback' % smiles)
        path.append('etkdgv3+random_coords')

    if AllChem.MMFFHasAllMoleculeParams(mol3d_heavy):
        method = 'MMFF94s'
    elif AllChem.UFFHasAllMoleculeParams(mol3d_heavy):
        method = 'UFF'
    else:
        raise ConversionError('no MMFF or UFF parameters for %r' % smiles)

    results = _optimise(mol3d_heavy, method, spec['max_iters'])
    best, best_energy = _best_converged(ids, results)
    if best is not None:
        path.append('converged@%d' % spec['max_iters'])
    else:
        results = _optimise(mol3d_heavy, method, spec['max_iters_extended'])
        best, best_energy = _best_converged(ids, results)
        if best is None:
            raise ConversionError(
                'no conformer converged for %r even at %d iterations'
                % (smiles, spec['max_iters_extended']))
        path.append('converged@%d' % spec['max_iters_extended'])
    return mol3d_heavy, best, best_energy, method, '+'.join(path)


def _rotate_3d(conformer, euler_degrees):
    """A rigid rotation. Distances and angles are untouched, so the molecule is
    the same one; only the frame it is written in changes."""
    import math
    ax, ay, az = (math.radians(a) for a in euler_degrees)
    cx, sx = math.cos(ax), math.sin(ax)
    cy, sy = math.cos(ay), math.sin(ay)
    cz, sz = math.cos(az), math.sin(az)
    for i in range(conformer.GetNumAtoms()):
        p = conformer.GetAtomPosition(i)
        # Rx then Ry then Rz
        y, z = p.y * cx - p.z * sx, p.y * sx + p.z * cx
        x = p.x * cy + z * sy
        z = -p.x * sy + z * cy
        x, y = x * cz - y * sz, x * sz + y * cz
        conformer.SetAtomPosition(i, (x, y, z))


def _graph3d_mol_uncached(smiles):
    """Embed, optimise, drop the hydrogens added only to drive the geometry.

    The heavy-atom set is therefore the same one Graph2D carries, which is what
    makes the Graph2D/Graph3D pair isolate its coordinates.
    """
    base = _clean_mol(smiles)
    with_h = Chem.AddHs(base)
    mol, conf_id, energy, method, path = _conformer(with_h, smiles)
    if not CONVERSION['graph3d'].get('keep_hydrogens', False):
        mol = Chem.RemoveHs(mol)
        if not mol.GetNumConformers():
            raise ConversionError('conformer lost when removing hydrogens for %r' % smiles)
        conf_id = 0
    return mol, conf_id, energy, method, path


def _graph3d_mol(smiles):
    """A private copy of the cached conformer, so a caller may rotate it freely."""
    mol, conf_id, energy, method, path = _cached(('graph3d', smiles),
                                           lambda: _graph3d_mol_uncached(smiles))
    return Chem.Mol(mol), conf_id, energy, method, path


def graph3d_text(smiles, variant=0):
    """MOL V2000 with computed 3D coordinates.

    Variant 0 is the conformer as it comes out of the optimiser; variant 1 is the
    same conformer rigidly rotated. The conformer itself is never re-chosen per
    call and never depends on the answer.
    """
    key = ('graph3d_block', smiles, variant)
    if key in _CACHE:
        CACHE_STATS['memory_hits'] += 1
        return _CACHE[key]
    entry = cache_entry(smiles)
    if entry is not None:
        CACHE_STATS['disk_hits'] += 1
        block = entry['block_v1' if variant else 'block_v0']
        _CACHE[key] = block
        return block
    CACHE_STATS['misses'] += 1

    def emit():
        mol, conf_id, energy, method, path = _graph3d_mol(smiles)
        out = {}
        for index, use_variant in ((0, 0), (1, 1)):
            work = Chem.Mol(mol)
            if use_variant:
                _rotate_3d(work.GetConformer(conf_id), VIEW_3D_EULER_DEGREES)
            block = _mol_block(work, conf_id=conf_id, three_d=True)
            _verify_block(block, smiles, 'Graph3D')
            if '3D' not in block.splitlines()[1]:
                raise ConversionError('Graph3D block is not marked 3D for %r' % smiles)
            out[index] = block
        conf = mol.GetConformer(conf_id)
        zs = [abs(conf.GetAtomPosition(i).z) for i in range(mol.GetNumAtoms())]
        store_entry(smiles, {
            'smiles': smiles,
            'conversion_version': CONVERSION['version'],
            'block_v0': out[0], 'block_v1': out[1],
            'method': method, 'energy_kcal': energy, 'conformer_path': path,
            'n_candidates': CONVERSION['graph3d']['n_candidates'],
            'is_3d_flag': bool(conf.Is3D()),
            'max_abs_z': max(zs) if zs else 0.0,
            'nearly_planar': (max(zs) < 0.1) if zs else True,
        })
        return out[variant]

    block = emit()
    _CACHE[key] = block
    return block


def graph3d_identity_variant(smiles, variant):
    return graph3d_text(smiles, variant)


def graph3d_geometry(smiles):
    """Conformer provenance and a planarity flag, for the manifest and coverage
    report. Never sent to the model."""
    entry = cache_entry(smiles)
    if entry is None:
        graph3d_text(smiles, 0)          # builds and stores the entry
        entry = cache_entry(smiles)
    if entry is None:
        raise ConversionError('no conformer entry could be produced for %r' % smiles)
    # a genuinely planar molecule is expected to be nearly flat; the flag
    # RECORDS that, it does not fail the molecule
    return {k: entry[k] for k in ('method', 'energy_kcal', 'conformer_path', 'n_candidates',
                                  'is_3d_flag', 'max_abs_z', 'nearly_planar')}


# --------------------------------------------------------------------------
# cml
# --------------------------------------------------------------------------

def cml_text(smiles, variant=0):
    return _cached(('cml', smiles, variant), lambda: _cml_text_uncached(smiles, variant))


def _cml_text_uncached(smiles, variant):
    mol = _clean_mol(smiles)
    if variant:
        # the same molecule written with a different atom labelling and a
        # different array order; the graph, the charges and the aromatic flags
        # are all unchanged
        mol = Chem.RenumberAtoms(mol, _renumber_order(mol.GetNumAtoms(), variant))
    Chem.Kekulize(mol, clearAromaticFlags=False)
    root = ET.Element('molecule')
    atoms = ET.SubElement(root, 'atomArray')
    for atom in mol.GetAtoms():
        attrs = {'id': 'a%d' % atom.GetIdx(), 'elementType': atom.GetSymbol()}
        if atom.GetFormalCharge():
            attrs['formalCharge'] = str(atom.GetFormalCharge())
        if atom.GetIsAromatic():
            attrs['aromatic'] = 'true'
        if atom.GetTotalNumHs() == 0 and atom.GetSymbol() in ('N', 'O'):
            attrs['hydrogenCount'] = '0'
        ET.SubElement(atoms, 'atom', attrs)
    bonds = ET.SubElement(root, 'bondArray')
    for bond in mol.GetBonds():
        order = {'SINGLE': 'S', 'DOUBLE': 'D', 'TRIPLE': 'T', 'AROMATIC': 'A'}.get(
            str(bond.GetBondType()), 'S')
        ET.SubElement(bonds, 'bond', {
            'id': 'b%d' % bond.GetIdx(),
            'atomRefs2': 'a%d a%d' % (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
            'order': order,
            'aromatic': 'true' if bond.GetIsAromatic() else 'false'})
    text = ET.tostring(root, encoding='unicode')
    if not text.startswith('<?xml'):
        text = '<?xml version="1.0" encoding="UTF-8"?>\n' + text
    ET.fromstring(text)          # XML syntax must be valid
    return text


def cml_roundtrip_mol(smiles):
    """Read our own CML back to a molecule, for the lossless check."""
    text = cml_text(smiles)
    root = ET.fromstring(text)
    rw = Chem.RWMol()
    index = {}
    for atom in root.find('atomArray'):
        a = Chem.Atom(atom.get('elementType'))
        if atom.get('formalCharge'):
            a.SetFormalCharge(int(atom.get('formalCharge')))
        if atom.get('aromatic') == 'true':
            a.SetIsAromatic(True)
        if atom.get('hydrogenCount') == '0':
            a.SetNoImplicit(True)
        index[atom.get('id')] = rw.AddAtom(a)
    order = {'S': Chem.BondType.SINGLE, 'D': Chem.BondType.DOUBLE,
             'T': Chem.BondType.TRIPLE, 'A': Chem.BondType.AROMATIC}
    for bond in root.find('bondArray'):
        left, right = bond.get('atomRefs2').split()
        rw.AddBond(index[left], index[right], order[bond.get('order')])
    return rw.GetMol()


def cml_identity_variant(smiles, variant):
    return cml_text(smiles, variant)


# --------------------------------------------------------------------------
# moljson
# --------------------------------------------------------------------------

def _moljson_payload(smiles, variant=0):
    from moljson import MolToJSON
    mol = _mol(smiles)
    if variant:
        mol = Chem.RenumberAtoms(mol, _renumber_order(mol.GetNumAtoms(), variant))
    spec = CONVERSION['moljson']
    return MolToJSON(mol,
                     atom_id_style=spec['atom_id_style'],
                     include_empty_fields=spec['include_empty_fields'])


def moljson_text(smiles, variant=0):
    return _cached(('moljson', smiles, variant),
                   lambda: _moljson_text_uncached(smiles, variant))


def _moljson_text_uncached(smiles, variant):
    from moljson import MolFromJSON, CheckRoundTrip
    mol = _mol(smiles)
    payload = _moljson_payload(smiles, variant)
    ok, _a, _b, _rt = CheckRoundTrip(mol)
    if not ok:
        raise ConversionError('MolJSON round-trip failed for %r' % smiles)
    required = set(moljson_schema_keys()['required'])
    missing = required - set(payload)
    if missing:
        raise ConversionError('MolJSON payload is missing schema key(s) %s for %r'
                              % (sorted(missing), smiles))
    try:
        back = MolFromJSON(payload)
    except Exception as exc:  # noqa: BLE001
        raise ConversionError('MolFromJSON rejected our payload for %r: %s' % (smiles, exc))
    if back is None or Chem.MolToSmiles(back) != Chem.MolToSmiles(mol):
        raise ConversionError('MolJSON round-trip gives a different molecule for %r' % smiles)
    return json.dumps(payload, ensure_ascii=False, separators=(',', ':'))


def moljson_schema_keys():
    path = MOLJSON_SRC.parent / 'schemas' / 'moljson.schema.json'
    data = json.loads(path.read_text(encoding='utf-8'))
    return {'required': data.get('required'), 'properties': list(data.get('properties', {}))}


def moljson_identity_variant(smiles, variant):
    return moljson_text(smiles, variant)


# --------------------------------------------------------------------------
# registry
# --------------------------------------------------------------------------

REPRESENTATION_KIND = {
    'smiles': 'text', 'selfies': 'text', 'our_code': 'text',
    'img': 'image', 'graph2d': 'text', 'graph3d': 'text',
    'cml': 'text', 'moljson': 'text',
}

EXTENSION_NAMES = ('img', 'graph2d', 'graph3d', 'cml', 'moljson')

RENDERERS = {
    'img': (img_png, img_identity_variant),
    'graph2d': (graph2d_text, graph2d_identity_variant),
    'graph3d': (graph3d_text, graph3d_identity_variant),
    'cml': (cml_text, cml_identity_variant),
    'moljson': (moljson_text, moljson_identity_variant),
}

PILOT_DISPLAY_NAME = {
    'smiles': 'SMILES', 'selfies': 'SELFIES', 'our_code': 'Our Code',
    'img': 'img', 'graph2d': 'Graph2D\nMOL 2D文本', 'graph3d': 'Graph3D\nMOL 3D文本',
    'cml': 'CML', 'moljson': 'MolJSON',
}


def render(name, smiles, variant=0):
    if name not in RENDERERS:
        raise KeyError('no extension renderer for %r' % name)
    fn = RENDERERS[name][1 if variant else 0]
    return fn(smiles, variant)


def payload_sha256(payload):
    blob = payload if isinstance(payload, bytes) else str(payload).encode('utf-8')
    return hashlib.sha256(blob).hexdigest()


def strip_metadata_png(blob):
    """Re-encode so the PNG contains pixels only."""
    from PIL import Image
    with Image.open(io.BytesIO(blob)) as image:
        buffer = io.BytesIO()
        image.convert('RGB').save(buffer, format='PNG', optimize=False)
        return buffer.getvalue()
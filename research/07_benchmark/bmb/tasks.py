"""Task definitions: what the model is asked, and how its answer is read.

The prompt is byte-identical across representations apart from two things:
the representation note (one sentence saying what notation this is) and the
payload itself. No prompt contains the answer, a molecule name, a reference
value, or the internal example id. Every prompt asks for a compact JSON object
so the endpoint's json_object mode is legal (the keyword 'JSON' is present).

The model never calls a tool, never executes Our Code, never runs RDKit.
"""
from __future__ import annotations

import hashlib
import json
import random

from . import labels, reps

TASK_ORDER = ('t1_identity', 't2_functional_groups', 't3_atom_counts', 't4_qed', 't5_hydrogen_bond_donors')

TASK_TITLES_CN = {
    't1_identity': 'Identity 同一性判断',
    't2_functional_groups': '子结构识别（官能团）',
    't3_atom_counts': '结构理解（重原子计数）',
    't4_qed': 'QED 估计（本项目新增任务）',
    't5_hydrogen_bond_donors': '基础描述符（HBD 计数）',
}

REP_NOTES = {
    'smiles': ('The molecules are given as SMILES strings. A SMILES string is a text '
               'walk over the molecular graph; branches are in parentheses, ring closures '
               'are digits, and lowercase letters mark aromatic atoms.'),
    'selfies': ('The molecules are given as SELFIES strings. SELFIES is a text notation '
                'in which every bracket group is an atom or branch symbol; any well-formed '
                'SELFIES string decodes to a valid molecule.'),
    'our_code': ('Each molecule is given as Python source that builds explicit records. '
                 'Atom(id=..., element=..., formal_charge=..., isotope=..., aromatic=..., '
                 'explicit_hs=..., atom_map=...) describes one atom; Bond(begin=..., end=..., '
                 'bond_type=...) describes one bond between two atom ids; '
                 'Molecule(atoms=(...), bonds=(...)) collects them. Atom.id is only a local '
                 'label inside that one expression and carries no chemical meaning. '
                 'explicit_hs=None means the element\'s normal implicit-hydrogen rule applies; '
                 'an integer is an exact hydrogen count.'),
    'moljson': ('The molecule is given as a MolJSON document: a JSON object with an "atoms" '
                'list (each entry has an id and an element), a "bonds" list (source, target, '
                'order), and an "aromatic_n_h" list naming aromatic nitrogens carrying a '
                'hydrogen.'),
    'cml': ('The molecule is given as a CML document, which is an XML format for chemical '
            'information. The <atomArray> holds one <atom> element per atom, with an id, an '
            'elementType, and optional formalCharge, hydrogenCount and aromatic attributes. '
            'The <bondArray> holds one <bond> element per bond, giving its atomRefs2 pair and '
            'its order (S=single, D=double, T=triple, A=aromatic). Atom ids are local labels '
            'inside this one document and carry no chemical meaning.'),
    'img': ('The molecule is given as a 2D structure image. Read the atoms and bonds from the '
            'picture; no text encoding of the molecule is provided.'),
    'graph2d': ('The molecule is given as a MOL V2000 connection table. The header line is '
                'blank; the counts line gives the number of atoms and bonds; the atom block '
                'lists one atom per line with its x, y, z coordinates and element symbol; the '
                'bond block lists one bond per line as atom1, atom2 and bond order (1=single, '
                '2=double, 3=triple, 4=aromatic). This is a 2D table: the z coordinate is 0. '
                'Atom indices are 1-based and only identify positions inside this table.'),
    'graph3d': ('The molecule is given as a MOL V2000 connection table with the same layout as '
                'a 2D one, except the atom block now carries real x, y, z coordinates in '
                'angstrom from a computed 3D conformer, and the header marks the table as 3D. '
                'The atoms are the same heavy atoms as in the 2D table; only the coordinates '
                'differ.'),
}

_IMAGE_PLACEHOLDER = '[the structure image is attached]'


def image_data_uri(blob):
    import base64
    return 'data:image/png;base64,' + base64.b64encode(blob).decode('ascii')


def user_message(rep, text, images=()):
    """A user turn: plain text, or a multimodal block list for image input.

    The image travels as a real image_url content block. It is never inlined as
    base64 in the text and never replaced by a text description of the picture,
    which would turn the image condition into a text one.
    """
    if not images:
        return {'role': 'user', 'content': text}
    blocks = [{'type': 'text', 'text': text}]
    for blob in images:
        blocks.append({'type': 'image_url', 'image_url': {'url': image_data_uri(blob)}})
    return {'role': 'user', 'content': blocks}


def payload_text(rep, payload):
    """(text shown to the model, image blobs). Text reps inline the payload."""
    if rep == 'img':
        return _IMAGE_PLACEHOLDER, (payload,)
    return payload, ()


class Task:
    name = ''
    title_cn = ''
    expects_pair = False

    def build_items(self, pool, count, seed, rep, split):
        raise NotImplementedError

    def messages(self, rep, item, run_config):
        raise NotImplementedError

    def parse(self, raw):
        raise NotImplementedError


# --------------------------------------------------------------------------
# JSON plumbing
# --------------------------------------------------------------------------

def extract_json(raw):
    """Pull the first JSON object out of a raw response. Returns (obj, error)."""
    if not isinstance(raw, str):
        return None, 'response was not a string'
    text = raw.strip()
    if text.startswith('```'):
        text = text.split('```')[1] if text.count('```') >= 2 else text.strip('`')
        if text.lstrip().lower().startswith('json'):
            text = text.lstrip()[4:]
    start = text.find('{')
    end = text.rfind('}')
    if start < 0 or end <= start:
        return None, 'no JSON object found'
    try:
        return json.loads(text[start:end + 1]), None
    except json.JSONDecodeError as exc:
        return None, 'invalid JSON: %s' % exc


def system_message():
    return {
        'role': 'system',
        'content': ('You are answering chemistry questions. You cannot run code, call tools, '
                    'or look anything up. Answer only from the structure you are given. '
                    'Reply with a single compact JSON object and nothing else.'),
    }


def _hash_messages(messages):
    blob = json.dumps(messages, ensure_ascii=False, sort_keys=True)
    return hashlib.sha256(blob.encode('utf-8')).hexdigest()


# --------------------------------------------------------------------------
# T1 Identity
# --------------------------------------------------------------------------

T1_INSTRUCTION = (
    'Structure A and structure B are two encodings of possibly the same molecule.\n'
    'Decide whether they describe the SAME molecule.\n'
    'Two encodings describe the same molecule when they have the same atoms, the same '
    'bonds between those atoms, the same bond orders and the same formal charges, '
    'regardless of which atom the text happens to start from or how the records are ordered.\n'
    'Reply with JSON only: {{"same_molecule": true}} or {{"same_molecule": false}}.'
)


class IdentityTask(Task):
    name = 't1_identity'
    title_cn = TASK_TITLES_CN['t1_identity']
    expects_pair = True

    def build_pairs(self, pool, count, seed):
        """The Identity pair list, decided WITHOUT rendering anything.

        This is the fix for silent substitution. The list of pairs is a property
        of the molecule pool and the seed alone, so every representation is
        asked exactly the same questions. Previously the pair list was built
        while rendering, and a molecule that could not be rendered was skipped
        in favour of the next candidate -- so Graph3D ended up with 200 pairs of
        which only 196 matched the other seven, without that being visible as a
        failure. A representation that cannot render an endpoint now simply does
        not answer that item; the item is reported as a conversion failure and
        the endpoints are never swapped for easier ones.
        """
        rng = random.Random(seed)
        n_pos = count // 2
        n_neg = count - n_pos
        by_formula, by_heavy = {}, {}
        for item in pool:
            by_formula.setdefault(item['formula'], []).append(item)
            by_heavy.setdefault(item['heavy_atoms'], []).append(item)

        pairs, used = [], set()
        candidates = pool[:]
        rng.shuffle(candidates)
        for entry in candidates:
            if len(pairs) >= n_pos:
                break
            if entry['cluster'] in used:
                continue
            used.add(entry['cluster'])
            pairs.append({
                'kind': 'pos',
                'example_id': '%s-pos' % entry['cid'],
                'cluster': entry['cluster'],
                'clusters': (entry['cluster'], entry['cluster']),
                'smiles_a': entry['smiles'], 'smiles_b': entry['smiles'],
                'label': {'same_molecule': True},
            })

        ordered = sorted(pool, key=lambda item: item['cluster'])
        rng.shuffle(ordered)
        for entry in ordered:
            if len(pairs) >= n_pos + n_neg:
                break
            partner = None
            for group in (by_formula.get(entry['formula'], []),
                          by_heavy.get(entry['heavy_atoms'], [])):
                pool_ids = [c for c in group if c['cluster'] != entry['cluster']]
                if pool_ids:
                    partner = rng.choice(pool_ids)
                    break
            if partner is None:
                continue
            key = tuple(sorted((entry['cluster'], partner['cluster'])))
            if key in used:
                continue
            used.add(key)
            pairs.append({
                'kind': 'neg',
                'example_id': '%s-neg' % entry['cid'],
                'cluster': entry['cluster'],
                'clusters': (entry['cluster'], partner['cluster']),
                'smiles_a': entry['smiles'], 'smiles_b': partner['smiles'],
                'label': {'same_molecule': False},
            })
        return pairs

    def build_items(self, pool, count, seed, rep, split, pair_views=None):
        """``pair_views=(va, vb)`` applies the SAME two views to slot A and slot B
        for a positive pair and a negative pair alike, so the view a molecule is
        drawn in carries no information about the answer.

        Left as None, the earlier three-representation construction is kept
        byte-for-byte (positives 0/1, negatives 0/0). That construction is the
        one the frozen main-02 responses were produced under, so changing it here
        would silently invalidate them instead of reusing them.
        """
        view_a, view_b = pair_views if pair_views else (0, 1)
        neg_view_b = view_b if pair_views else 0
        items = []
        # the PAIR LIST comes from the representation-independent builder above;
        # this loop only renders it. A pair that cannot be rendered is dropped
        # and reported, never replaced by the next candidate.
        for pair in self.build_pairs(pool, count, seed):
            variant_b = view_b if pair['label']['same_molecule'] else neg_view_b
            try:
                a = reps.render(rep, pair['smiles_a'], variant=view_a)
                b = reps.render(rep, pair['smiles_b'], variant=variant_b)
            except Exception:  # noqa: BLE001
                continue
            if pair['label']['same_molecule'] and a == b:
                # a positive whose two views are identical would be answerable
                # by string equality; it is dropped rather than substituted
                continue
            items.append({
                'example_id': 't1-%s-%s' % (rep, pair['example_id']),
                'task': self.name, 'rep': rep, 'split': split,
                'cluster': pair['cluster'],
                # both endpoints, so cluster resampling cannot ignore one side
                'clusters': pair['clusters'],
                'payload': {'a': a, 'b': b},
                'label': pair['label'],
            })
        return items

    def messages(self, rep, item, run_config):
        text_a, images_a = payload_text(rep, item['payload']['a'])
        text_b, images_b = payload_text(rep, item['payload']['b'])
        content = (
            'Representation: %s\n\nA:\n%s\n\nB:\n%s\n\n%s'
            % (REP_NOTES[rep], text_a, text_b, T1_INSTRUCTION.format())
        )
        return [system_message(),
                user_message(rep, content, images_a + images_b)]

    def parse(self, raw):
        obj, error = extract_json(raw)
        if error:
            return None, error
        if not isinstance(obj, dict) or 'same_molecule' not in obj:
            return None, 'missing key same_molecule'
        value = obj['same_molecule']
        if isinstance(value, str):
            low = value.strip().lower()
            if low in ('true', 'yes'):
                value = True
            elif low in ('false', 'no'):
                value = False
        if not isinstance(value, bool):
            return None, 'same_molecule is not a boolean'
        return {'same_molecule': value}, None


# --------------------------------------------------------------------------
# T2 functional groups
# --------------------------------------------------------------------------

class FunctionalGroupTask(Task):
    name = 't2_functional_groups'
    title_cn = TASK_TITLES_CN['t2_functional_groups']

    def label_names(self):
        return [item['name'] for item in labels.enabled_labels()]

    def build_items(self, pool, count, seed, rep, split):
        rng = random.Random(seed)
        chosen = pool[:]
        rng.shuffle(chosen)
        items = []
        for entry in chosen[:count]:
            try:
                payload = reps.render(rep, entry['smiles'], variant=0)
            except Exception:  # noqa: BLE001
                continue
            items.append({
                'example_id': 't2-%s-%s' % (rep, entry['cid']),
                'task': self.name, 'rep': rep, 'split': split,
                'cluster': entry['cluster'],
                'payload': {'molecule': payload},
                'label': labels.functional_group_labels(entry['smiles'], only_enabled=True),
            })
        return items

    def messages(self, rep, item, run_config):
        names = self.label_names()
        listing = ', '.join('"%s"' % n for n in names)
        instruction = (
            'For the molecule below, decide for EACH functional group in the list whether it '
            'is present. Judge only from the structure given.\n'
            'Groups: %s\n'
            'Reply with JSON only, one boolean per group, for example '
            '{{"%s": true}}.' % (listing, names[0])
        )
        text, images = payload_text(rep, item['payload']['molecule'])
        content = ('Representation: %s\n\nMolecule:\n%s\n\n%s'
                   % (REP_NOTES[rep], text, instruction))
        return [system_message(), user_message(rep, content, images)]

    def parse(self, raw):
        obj, error = extract_json(raw)
        if error:
            return None, error
        if not isinstance(obj, dict):
            return None, 'not a JSON object'
        out = {}
        for name in self.label_names():
            if name not in obj:
                return None, 'missing group %s' % name
            value = obj[name]
            if isinstance(value, str):
                low = value.strip().lower()
                value = {'true': True, 'yes': True, 'false': False, 'no': False}.get(low, value)
            if not isinstance(value, bool):
                return None, 'group %s is not a boolean' % name
            out[name] = value
        return out, None


# --------------------------------------------------------------------------
# T3 atom counts
# --------------------------------------------------------------------------

class AtomCountTask(Task):
    name = 't3_atom_counts'
    title_cn = TASK_TITLES_CN['t3_atom_counts']

    def elements(self):
        return list(labels.COUNT_ELEMENTS)

    def build_items(self, pool, count, seed, rep, split):
        rng = random.Random(seed)
        chosen = pool[:]
        rng.shuffle(chosen)
        items = []
        for entry in chosen[:count]:
            try:
                payload = reps.render(rep, entry['smiles'], variant=0)
            except Exception:  # noqa: BLE001
                continue
            items.append({
                'example_id': 't3-%s-%s' % (rep, entry['cid']),
                'task': self.name, 'rep': rep, 'split': split,
                'cluster': entry['cluster'],
                'payload': {'molecule': payload},
                'label': labels.atom_counts(entry['smiles']),
            })
        return items

    def messages(self, rep, item, run_config):
        elements = self.elements()
        instruction = (
            'Count how many atoms of each listed element are in the molecule below. '
            'Count heavy atoms only; do not count hydrogen atoms.\n'
            'Elements: %s\n'
            'Reply with JSON only, one non-negative integer per element, for example '
            '{{"%s": 0}}.' % (', '.join(elements), elements[0])
        )
        text, images = payload_text(rep, item['payload']['molecule'])
        content = ('Representation: %s\n\nMolecule:\n%s\n\n%s'
                   % (REP_NOTES[rep], text, instruction))
        return [system_message(), user_message(rep, content, images)]

    def parse(self, raw):
        obj, error = extract_json(raw)
        if error:
            return None, error
        if not isinstance(obj, dict):
            return None, 'not a JSON object'
        out = {}
        for element in self.elements():
            if element not in obj:
                return None, 'missing element %s' % element
            value = obj[element]
            if isinstance(value, bool) or not isinstance(value, int):
                if isinstance(value, str) and value.strip().lstrip('-').isdigit():
                    value = int(value.strip())
                else:
                    return None, 'element %s is not an integer' % element
            if value < 0:
                return None, 'element %s is negative' % element
            out[element] = value
        return out, None


# --------------------------------------------------------------------------
# T4 QED
# --------------------------------------------------------------------------

class QedTask(Task):
    name = 't4_qed'
    title_cn = TASK_TITLES_CN['t4_qed']

    def build_items(self, pool, count, seed, rep, split):
        rng = random.Random(seed)
        chosen = pool[:]
        rng.shuffle(chosen)
        items = []
        for entry in chosen[:count]:
            try:
                payload = reps.render(rep, entry['smiles'], variant=0)
            except Exception:  # noqa: BLE001
                continue
            items.append({
                'example_id': 't4-%s-%s' % (rep, entry['cid']),
                'task': self.name, 'rep': rep, 'split': split,
                'cluster': entry['cluster'],
                'payload': {'molecule': payload},
                'label': {'qed': labels.qed_reference(entry['smiles'])},
            })
        return items

    def messages(self, rep, item, run_config):
        instruction = (
            'Estimate the drug-likeness score QED (Quantitative Estimate of Drug-likeness) '
            'for the molecule below, on its usual scale from 0 (least drug-like) to 1 '
            '(most drug-like). Answer from the structure alone; do not run any code and do '
            'not attempt to compute descriptors externally.\n'
            'Reply with JSON only: {"qed": <a number between 0 and 1>}.'
        )
        text, images = payload_text(rep, item['payload']['molecule'])
        content = ('Representation: %s\n\nMolecule:\n%s\n\n%s'
                   % (REP_NOTES[rep], text, instruction))
        return [system_message(), user_message(rep, content, images)]

    def parse(self, raw):
        obj, error = extract_json(raw)
        if error:
            return None, error
        if not isinstance(obj, dict) or 'qed' not in obj:
            return None, 'missing key qed'
        value = obj['qed']
        if isinstance(value, bool):
            return None, 'qed is a boolean'
        if isinstance(value, str):
            try:
                value = float(value.strip())
            except ValueError:
                return None, 'qed is not a number'
        if not isinstance(value, (int, float)):
            return None, 'qed is not a number'
        value = float(value)
        if not 0.0 <= value <= 1.0:
            return None, 'qed outside [0, 1]: %r' % value
        return {'qed': value}, None


# --------------------------------------------------------------------------
# T5 descriptors
# --------------------------------------------------------------------------

class DescriptorTask(Task):
    name = 't5_hydrogen_bond_donors'
    title_cn = TASK_TITLES_CN['t5_hydrogen_bond_donors']

    def build_items(self, pool, count, seed, rep, split):
        rng = random.Random(seed)
        chosen = pool[:]
        rng.shuffle(chosen)
        items = []
        for entry in chosen[:count]:
            try:
                payload = reps.render(rep, entry['smiles'], variant=0)
            except Exception:  # noqa: BLE001
                continue
            reference = labels.descriptor_reference(entry['smiles'])
            items.append({
                'example_id': 't5-%s-%s' % (rep, entry['cid']),
                'task': self.name, 'rep': rep, 'split': split,
                'cluster': entry['cluster'],
                'payload': {'molecule': payload},
                'label': {'hbd': reference['hbd']},
            })
        return items

    def messages(self, rep, item, run_config):
        instruction = (
            'Count the hydrogen-bond donors (HBD) in the molecule below, using the standard '
            'definition: the number of nitrogen or oxygen atoms that carry at least one '
            'hydrogen atom.\n'
            'Reply with JSON only: {"hbd": <a non-negative integer>}.'
        )
        text, images = payload_text(rep, item['payload']['molecule'])
        content = ('Representation: %s\n\nMolecule:\n%s\n\n%s'
                   % (REP_NOTES[rep], text, instruction))
        return [system_message(), user_message(rep, content, images)]

    def parse(self, raw):
        obj, error = extract_json(raw)
        if error:
            return None, error
        if not isinstance(obj, dict) or 'hbd' not in obj:
            return None, 'missing key hbd'
        value = obj['hbd']
        if isinstance(value, bool):
            return None, 'hbd is a boolean'
        if isinstance(value, str) and value.strip().lstrip('-').isdigit():
            value = int(value.strip())
        if not isinstance(value, int):
            return None, 'hbd is not an integer'
        if value < 0:
            return None, 'hbd is negative'
        return {'hbd': value}, None


TASKS = {
    't1_identity': IdentityTask(),
    't2_functional_groups': FunctionalGroupTask(),
    't3_atom_counts': AtomCountTask(),
    't4_qed': QedTask(),
    't5_hydrogen_bond_donors': DescriptorTask(),
}


def get(name: str) -> Task:
    if name not in TASKS:
        raise KeyError('unknown task %r; known: %s' % (name, ', '.join(TASK_ORDER)))
    return TASKS[name]


def prompt_hash(rep: str, item, run_config) -> str:
    return _hash_messages(get(item['task']).messages(rep, item, run_config))
**Yes. Define your own `parse_smiles()` and `write_smiles()` utilities, and make `Molecule` depend on those utilities rather than RDKit.** The utilities themselves can be implemented in ordinary Python.

There are two different outcomes: moving RDKit calls into `utils.py` merely hides the dependency; implementing the utilities independently removes it. **I built the second version.**

molecule_standalone.py · molecule_no_rdkit.zip

**Important:** this is a working, connectivity-focused prototype, not a complete RDKit replacement. It rejects specified stereochemistry and does not perform full chemical validation or canonicalization.

## 1. What should replace the RDKit imports?

I would separate three responsibilities:

| Component                          | Responsibility                                                              | Implementation in the new prototype            |
| ---------------------------------- | --------------------------------------------------------------------------- | ---------------------------------------------- |
| `Atom`, `Bond`, `Molecule`         | Store the molecular structure and expose operations                         | Plain Python classes                           |
| `parse_smiles()`, `write_smiles()` | Convert between SMILES and explicit graph records                           | Local Python utilities                         |
| Chemical validation                | Check valence, aromaticity, stereochemistry, and other chemical constraints | Separate; not comprehensively implemented here |

The dependency structure becomes:

```
SMILES
   │
   │ parse_smiles() — your own Python function
   ▼
Molecule(atoms=..., bonds=...)
   │
   │ write_smiles() — your own Python function
   ▼
SMILES
```

In the supplied implementation, all three classes and their utilities live in **one source file**. Its imports are only from Python’s standard library: `dataclasses`, `collections`, `typing`, `ast`, and `re`.

This requires new implementations rather than copying the RDKit Python wrappers. RDKit’s underlying SMILES parser is C++ and depends on its graph objects, stereochemistry components, and generated parsing machinery. Its writer likewise operates on RDKit’s molecular infrastructure. [GitHub](https://raw.githubusercontent.com/rdkit/rdkit/master/Code/GraphMol/SmilesParse/SmilesParse.cpp)

Borrowing the graph design and implementing the required serialization algorithms is therefore a cleaner approach than trying to extract a few RDKit functions.

## 2. How the class works without RDKit

The following runs with the downloaded module and requires no third-party installation:

```
from molecule_standalone import Atom, Bond, Molecule

mol = Molecule(smiles="CCO")

assert mol.decode_to_smiles() == "CCO"
assert mol.neighbors("a2") == ("a1",)
assert mol.shortest_path("a0", "a2") == ("a0", "a1", "a2")

# Change the explicit molecular structure.
edited = mol.replace_atom("a2", element="S")

assert edited.decode_to_smiles() == "CCS"
assert mol.decode_to_smiles() == "CCO"
```

The import above loads **your local source file**, not RDKit.

Calling `mol.to_code()` produces:

```
Molecule(
    atoms=(
        Atom(id="a0", element="C"),
        Atom(id="a1", element="C"),
        Atom(id="a2", element="O"),
    ),
    bonds=(
        Bond(begin="a0", end="a1"),
        Bond(begin="a1", end="a2"),
    ),
)
```

**This constructor expression can be the LLM-facing molecular representation.** Provide the shared class and utility definitions once; individual molecular examples then contain no imports.

The conversion also works in the other direction:

```
code = mol.to_code()
reconstructed = Molecule.from_code(code)

assert reconstructed.decode_to_smiles() == "CCO"
```

`from_code()` reads a restricted constructor expression through an abstract syntax tree. It does not use `eval()` or execute arbitrary Python statements.

Most importantly, the object stores **only its atom and bond records**. Decoding traverses those records rather than returning a saved input string. Consequently, modifying the representation changes the decoded molecule.

## 3. What the local utilities actually do

### Parsing: SMILES → explicit graph

`parse_smiles()` maintains a current atom, a branch stack, a pending bond type, and a table of open ring labels. It creates atom records as it reads atom tokens and creates bond records as it encounters connections.

This follows the SMILES distinction between atoms, bond symbols, branches, ring closures, and disconnected components. Those constructs are specified independently of RDKit in OpenSMILES. [GitHub](https://raw.githubusercontent.com/opensmiles/OpenSMILES/master/opensmiles.asciidoc)

### Writing: explicit graph → SMILES

`write_smiles()` builds a depth-first traversal of the current graph. Tree edges become the main traversal and branches; non-tree edges become numbered ring closures. Separate connected components are joined with dots.

The writer does **not** need to recover the original traversal. Its responsibility is to serialize the current graph. In this prototype, its output is deterministic for a fixed ordered graph, but it is **not canonical**.

### Queries: explicit graph → answers

Functions such as `neighbors()`, `bond_between()`, and `shortest_path()` operate directly on the atom and bond records.

For your research goal, I would make these operations inspectable. The LLM can then reason over both the molecular data and the algorithms that access it, instead of calling an opaque chemistry object.

## 4. One important attribute: hydrogen policy

The new `Atom` fields are:

```
id
element
formal_charge
isotope
aromatic
explicit_hs
atom_map
```

Here, `explicit_hs` distinguishes an inferred hydrogen count from a specified bracket hydrogen count:

| Input atom | Stored `explicit_hs` | Interpretation                                  |
| ---------- | -------------------- | ----------------------------------------------- |
| `C`        | `None`               | Use unbracketed-style implicit-hydrogen rules   |
| `[C]`      | `0`                  | Exactly zero bracket hydrogens                  |
| `[CH3]`    | `3`                  | Exactly three bracket hydrogens                 |
| `[nH]`     | `1`                  | One hydrogen specified on the aromatic nitrogen |

This distinction is necessary: under OpenSMILES, omitting the hydrogen count **inside brackets** means zero, whereas an unbracketed organic-subset atom uses implicit-hydrogen rules. Hydrogen atoms written as `[H]` are separate graph vertices. [GitHub](https://raw.githubusercontent.com/opensmiles/OpenSMILES/master/opensmiles.asciidoc)

For example, ammonium can be constructed directly:

```
ammonium = Molecule(
    atoms=(
        Atom(
            id="n",
            element="N",
            formal_charge=1,
            explicit_hs=4,
        ),
    ),
)

assert ammonium.decode_to_smiles() == "[NH4+]"
```

The prototype deliberately requires an explicit hydrogen policy for manually created atoms that require brackets. It does not silently guess when a charge or isotope annotation forces bracket notation.

## 5. What is supported—and what still needs work

| Feature                                                   | Current standalone implementation          |
| --------------------------------------------------------- | ------------------------------------------ |
| Chains, branches, rings, disconnected components          | Implemented                                |
| Single, double, triple, and aromatic bonds                | Implemented                                |
| Elements, isotopes, charges, atom maps, bracket hydrogens | Implemented                                |
| Aromatic notation                                         | Preserved; aromaticity is not recalculated |
| Specified stereochemistry                                 | Explicitly rejected                        |
| Canonical SMILES and full chemical sanitization           | Not implemented                            |

**Stereochemistry requires a proper extension, not just an extra string attribute.** SMILES tetrahedral markers are interpreted relative to neighbor ordering. Reordering the traversal while copying the original `@` or `@@` can change the represented configuration. The extension should therefore store local stereo references and adjust the serialized marker when the traversal changes. [GitHub](https://raw.githubusercontent.com/opensmiles/OpenSMILES/master/opensmiles.asciidoc)

Until that is implemented, the codec raises `UnsupportedSmilesError` for `@`, `/`, and `\`. It never silently removes them.

**Chemical validity must also remain separate from graph consistency.** The current validator checks atom IDs, bond endpoints, duplicate edges, and field consistency. RDKit sanitization does much more, including valence checks, Kekulization, radical assignment, aromaticity perception, and hybridization assignment. [RDKit](https://www.rdkit.org/docs/RDKit_Book.html)

A structurally consistent five-bond carbon graph can therefore be serialized by this prototype. That does not make it chemically valid.

## 6. Verification performed

The standard-library-only suite passed **16 test methods**, including **60 named SMILES cases**, **200 atom/bond-order permutations**, code reconstruction, structural edits, and explicit rejection of unsupported stereochemistry. It also passed with Python’s site initialization disabled. unit_test_output.txt

Separately, I used RDKit **only as an external development comparator**. All **4,991 RDKit-accepted molecules** in its bundled NCI sample matched after direct decoding and after shuffled atom/bond ordering. The other **eight records rejected by RDKit were accepted by this graph codec**, which demonstrates the chemical-validation limitation. These results do not establish full SMILES coverage. verification_results.json

The package includes the optional comparison script in a separate directory. **Only that optional script requires RDKit; the molecular class, conversion utilities, demo, and main tests do not.**

An existing alternative worth knowing is **PySMILES**, a pure-Python reader/writer that uses NetworkX. However, its documentation states that the writer ignores stereochemical information, so it should not be substituted silently when stereo preservation matters. [GitHub](https://github.com/pckroon/pysmiles?utm_source=chatgpt.com)

**My recommendation for your project is to keep `Atom`, `Bond`, and `Molecule` fully independent, use local utilities for conversion, and treat chemical validation as a separate interface.** That lets the first experiment focus on LLM reasoning over explicit molecular code without making a complete chemistry toolkit a prerequisite.

---

如有需要，我可以：

- Explain how parse_smiles() processes SMILES strings in detail
- Demonstrate editing a molecule and re-encoding to SMILES
- Describe limitations related to stereochemistry and validation in this prototype

molecule_standalone.py程式碼開啟檔案unit_test_output.txt文件開啟檔案verification_results.json程式碼
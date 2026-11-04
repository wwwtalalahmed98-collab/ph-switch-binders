"""Build the three target structures the counter-screen scores every design against.

    g12d   HLA + b2m + VVVGADGVGK     must bind
    wt     HLA + b2m + VVVGAGGVGK     must not bind
    apo    HLA + b2m, no peptide      must not bind

THE WILD TYPE IS MODELLED, NOT CRYSTALLOGRAPHIC, AND THAT IS DELIBERATE. 7OW5 is the
real wild-type structure, but using it would mean the mutant and wild-type states differ
by a whole crystal as well as by one residue, and any margin we measured could come from
either. So the wild-type state is 7OW6 with Asp6 cut back to Gly - side-chain atoms
beyond the backbone removed and the residue renamed - leaving the single substitution as
the only difference between the two states. 7OW5 is still used to fix the shared residue
set, as in phase 0.

The cost of that choice: the wild-type state carries the mutant's backbone. A real
wild-type complex may sit slightly differently. This buys a controlled comparison at the
price of a small structural inaccuracy, and the direction of that trade is recorded here
because it is not obviously the right one - it is the one that makes the margin
interpretable.

THE APO STATE IS THE WEAKEST OF THE THREE. Deleting the peptide leaves a groove still
shaped as though a peptide were in it. Real peptide-free HLA is conformationally
unstable and partly collapsed. So our apo test is easier than the assay's, and a design
that passes it has not been shown to ignore empty HLA. That belongs in the methods text,
not in a footnote.
"""
import copy
import json
import os

from Bio import PDB

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data", "pmhc")
OUT = os.path.join(D, "states")

MUT_POS = 6
PEP_G12D = "VVVGADGVGK"
PEP_WT = "VVVGAGGVGK"
BACKBONE = ("N", "CA", "C", "O", "OXT")   # what a glycine keeps

T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
      "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
      "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
      "TYR": "Y", "VAL": "V"}


def residues(chain):
    out, seen = [], set()
    for r in chain:
        if r.id[0] == " " and r.get_resname() in T3 and r.id[1] not in seen:
            seen.add(r.id[1])
            out.append(r)
    return out


def seq(rs):
    return "".join(T3[r.get_resname()] for r in rs)


class Sel(PDB.Select):
    """Keep the listed residues of the listed chains."""

    def __init__(self, keep):
        self.keep = {c: set(n) for c, n in keep.items()}

    def accept_model(self, m):
        return m.id == 0

    def accept_chain(self, c):
        return c.id in self.keep

    def accept_residue(self, r):
        return r.id[0] == " " and r.id[1] in self.keep[r.get_parent().id]

    def accept_atom(self, a):
        return (not a.is_disordered()) or a.get_altloc() in (" ", "A")


def main():
    man = json.load(open(os.path.join(OUT, "states.json")))
    m6 = PDB.MMCIFParser(QUIET=True).get_structure("7OW6",
                                                   os.path.join(D, "7OW6.cif"))[0]
    hv_ok, b2m_ok = set(man["kept"]["A"]), set(man["kept"]["B"])
    hv = [r for r in residues(m6["A"]) if r.id[1] in hv_ok]
    b2m = [r for r in residues(m6["B"]) if r.id[1] in b2m_ok]
    pep = residues(m6["C"])
    if seq(pep) != PEP_G12D:
        raise SystemExit("7OW6 chain C is %s, expected %s" % (seq(pep), PEP_G12D))

    keep_holo = {"A": [r.id[1] for r in hv], "B": [r.id[1] for r in b2m],
                 "C": [r.id[1] for r in pep]}
    keep_apo = {k: v for k, v in keep_holo.items() if k != "C"}
    mut_resnum = pep[MUT_POS - 1].id[1]

    written = {}
    io = PDB.PDBIO()
    io.set_structure(m6)
    for tag, keep in (("target_state_g12d", keep_holo),
                      ("target_state_apo", keep_apo)):
        path = os.path.join(OUT, "%s.pdb" % tag)
        io.save(path, Sel(keep))
        written[tag] = path

    # The wild-type state needs the residue RENAMED as well as stripped - dropping the
    # side chain alone leaves it reading as Asp, which is how the first version of this
    # script silently produced two identical states. Select cannot rename, so mutate a
    # copy of the model and save that.
    mut = copy.deepcopy(m6)
    target = mut["C"][mut_resnum]
    for a in [a.get_name() for a in target if a.get_name() not in BACKBONE]:
        target.detach_child(a)
    target.resname = "GLY"
    io.set_structure(mut)
    path = os.path.join(OUT, "target_state_wt.pdb")
    io.save(path, Sel(keep_holo))
    written["target_state_wt"] = path

    # ---- verify what was written, rather than trusting the selector ----
    print("%-20s %6s %6s %6s  %s" % ("state", "heavy", "b2m", "pep", "peptide sequence"))
    report = {}
    for tag, path in written.items():
        chk = PDB.PDBParser(QUIET=True).get_structure("t", path)[0]
        ch = {c.id: residues(c) for c in chk}
        pseq = seq(ch["C"]) if "C" in ch else "-"
        print("%-20s %6d %6d %6d  %s"
              % (tag, len(ch.get("A", [])), len(ch.get("B", [])),
                 len(ch.get("C", [])), pseq))
        report[tag] = dict(heavy=len(ch.get("A", [])), b2m=len(ch.get("B", [])),
                           peptide=len(ch.get("C", [])), peptide_seq=pseq)

    ok = True
    if report["target_state_g12d"]["peptide_seq"] != PEP_G12D:
        print("FAIL: g12d state peptide is not %s" % PEP_G12D)
        ok = False
    if report["target_state_wt"]["peptide_seq"] != PEP_WT:
        print("FAIL: wt state peptide is %s, not %s"
              % (report["target_state_wt"]["peptide_seq"], PEP_WT))
        ok = False
    if report["target_state_apo"]["peptide"] != 0:
        print("FAIL: apo state still has a peptide")
        ok = False
    for a, b in (("target_state_g12d", "target_state_wt"),
                 ("target_state_g12d", "target_state_apo")):
        if report[a]["heavy"] != report[b]["heavy"] or report[a]["b2m"] != report[b]["b2m"]:
            print("FAIL: %s and %s do not share an identical background" % (a, b))
            ok = False

    # the mutated residue must have lost exactly its side chain and nothing else
    g = PDB.PDBParser(QUIET=True).get_structure("g", written["target_state_g12d"])[0]
    w = PDB.PDBParser(QUIET=True).get_structure("w", written["target_state_wt"])[0]
    g_at = {a.get_name() for a in g["C"][mut_resnum]}
    w_at = {a.get_name() for a in w["C"][mut_resnum]}
    print("\nresidue %d: g12d has %s" % (mut_resnum, sorted(g_at)))
    print("residue %d: wt   has %s" % (mut_resnum, sorted(w_at)))
    if not w_at <= set(BACKBONE):
        print("FAIL: wt residue kept non-backbone atoms %s" % sorted(w_at - set(BACKBONE)))
        ok = False
    if "CA" not in w_at or "N" not in w_at or "C" not in w_at:
        print("FAIL: wt residue lost backbone atoms")
        ok = False

    man["counter_states"] = dict(
        files={k: os.path.basename(v) for k, v in written.items()},
        wild_type="modelled from 7OW6 by truncating Asp6 to Gly, so the only difference "
                  "from the g12d state is that one substitution",
        apo_caveat="peptide deleted from a holo structure; the groove keeps its "
                   "peptide-bound shape, so this is an easier test than the assay",
        counts=report)
    json.dump(man, open(os.path.join(OUT, "states.json"), "w"), indent=1)

    print("\n%s" % ("PASS - three states share one background and differ only as intended"
                    if ok else "TEST FAILED"))
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()

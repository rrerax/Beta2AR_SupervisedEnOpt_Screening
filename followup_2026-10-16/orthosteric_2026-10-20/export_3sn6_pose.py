#!/usr/bin/env python3
"""Write the already-accepted 3SN6 orciprenaline pose as PDB files.

Chi1 is chosen by the same representative rule as place_representative.py,
before the ligand moves. The script stops at the first arrangement that
passes the written hydrogen-bond and 2.6 A rules, then saves coordinates.
"""

import math
from pathlib import Path

import place_representative as pr
from rdkit import Chem
from rdkit.Chem import AllChem

HERE = Path(__file__).resolve().parent
PAIR = ("203", "207")


def pdb_line(serial, record, name, resn, chain, resi, xyz, element):
    atom = f"{name:>4}" if len(name) < 4 else name[:4]
    return (
        f"{record:<6}{serial:5d} {atom} {resn:>3} {chain:1}{int(resi):4d}    "
        f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {element:>2}"
    )


def write_protein(path, protein):
    lines = ["HEADER    3SN6 CHAIN R WITH SER203 CHI1 AND SER207 CHI1 REPLACED"]
    for serial, (resn, resi, name, xyz) in enumerate(protein, start=1):
        element = "".join(ch for ch in name if ch.isalpha())[:1] or "C"
        lines.append(pdb_line(serial, "ATOM", name, resn, "R", resi, xyz, element))
    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_ligand(path, mol, positions):
    lines = ["HEADER    PROTONATED ORCIPRENALINE ON THE PASSING 3SN6 POSE"]
    for idx, atom in enumerate(mol.GetAtoms()):
        pos = positions[idx]
        lines.append(
            pdb_line(idx + 1, "HETATM", atom.GetSymbol(), "LIG", "L", 1, pos, atom.GetSymbol())
        )
    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def full_positions(mol, cid, rot, trans):
    conf = mol.GetConformer(int(cid))
    out = []
    for idx in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(int(idx))
        out.append(pr.sp.apply_rt(rot, trans, (pos.x, pos.y, pos.z)))
    return out


def main():
    path = pr.sp.sc.ROOT / pr.FILES["3SN6"]
    serines = pr.sp.motif_serines(path)
    asp = pr.find_asp(path)
    _, protein, backbone, side, asp_atoms = pr.load_side(path, serines, asp)
    clear = pr.passing_angles(backbone, side, PAIR)
    chosen = pr.representative(clear)
    chi = {PAIR[0]: chosen[0], PAIR[1]: chosen[1]}
    if chi != {"203": -167, "207": 49}:
        raise SystemExit(f"representative changed: {chi}")
    ogs = {crystal: pr.og_at(side, crystal, chi[crystal]) for crystal in PAIR}
    posed, cells = pr.receptor_with_ogs(protein, side, PAIR, ogs)
    write_protein(HERE / "3sn6_ser203_207_receptor.pdb", posed)

    mol = Chem.AddHs(Chem.MolFromSmiles(pr.sp.SMILES))
    params = AllChem.ETKDGv3()
    params.randomSeed = 20261020
    params.numThreads = 1
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=12, params=params))
    n_idx = next(atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "N")
    phenolic, phenol_c, phenol_h = [], {}, {}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "O":
            continue
        carbons = [nbr for nbr in atom.GetNeighbors() if nbr.GetIsAromatic()]
        if carbons:
            phenolic.append(atom.GetIdx())
            phenol_c[atom.GetIdx()] = carbons[0].GetIdx()
            phenol_h[atom.GetIdx()] = next(nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetSymbol() == "H")
    amine_h = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(n_idx).GetNeighbors() if nbr.GetSymbol() == "H"]
    heavy = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() != "H"]
    confs = []
    for cid in ids:
        embedded = mol.GetConformer(int(cid))

        def xyz(idx, embedded=embedded):
            pos = embedded.GetAtomPosition(int(idx))
            return (pos.x, pos.y, pos.z)

        needed = set(heavy + amine_h + list(phenol_h.values()) + list(phenol_c.values()))
        confs.append(
            {
                "n": n_idx,
                "phenolic": phenolic,
                "phenol_c": phenol_c,
                "phenol_h": phenol_h,
                "amine_h": amine_h,
                "heavy": heavy,
                "pts": {idx: xyz(idx) for idx in needed},
            }
        )

    ods = [asp_atoms["OD1"], asp_atoms["OD2"]]
    clouds = {radius: (pr.sphere(ogs["203"], radius), pr.sphere(ogs["207"], radius)) for radius in (2.8, 3.0, 3.2)}
    tested = 0
    for radius, (cloud_a, cloud_b) in clouds.items():
        pairs = [(pr.sp.sc._dist(p1, p2), p1, p2) for p1 in cloud_a for p2 in cloud_b]
        for cid, conf in zip(ids, confs):
            pts = conf["pts"]
            oo = pr.sp.sc._dist(pts[conf["phenolic"][0]], pts[conf["phenolic"][1]])
            radii = {
                conf["phenolic"][0]: pr.sp.sc._dist(pts[conf["n"]], pts[conf["phenolic"][0]]),
                conf["phenolic"][1]: pr.sp.sc._dist(pts[conf["n"]], pts[conf["phenolic"][1]]),
            }
            near = [row for row in pairs if abs(row[0] - oo) <= 0.3]
            for _, p1, p2 in near:
                for oa, ob, t1, t2 in (
                    (conf["phenolic"][0], conf["phenolic"][1], p1, p2),
                    (conf["phenolic"][1], conf["phenolic"][0], p1, p2),
                ):
                    for target_n in pr.nitrogen_circle(t1, radii[oa], t2, radii[ob]):
                        if pr.sp.sc._dist(target_n, asp_atoms["CG"]) < 2.4:
                            continue
                        if not 2.6 <= min(pr.sp.sc._dist(target_n, od) for od in ods) <= 3.5:
                            continue
                        tested += 1
                        rot, trans = pr.sp.kabsch((pts[conf["n"]], pts[oa], pts[ob]), (target_n, t1, t2))
                        placed = {idx: pr.sp.apply_rt(rot, trans, point) for idx, point in pts.items()}
                        amine_hits = [
                            pr.sp.hbond(placed[conf["n"]], placed[h], od) for h in conf["amine_h"] for od in ods
                        ]
                        amine = max(amine_hits, key=lambda hit: (hit[0], -hit[1], hit[2]))
                        phenols = []
                        for oxy, target in ((oa, ogs["203"]), (ob, ogs["207"])):
                            oxygen = placed[oxy]
                            hydrogen = placed[conf["phenol_h"][oxy]]
                            flipped = pr.sp.flip_h(placed[conf["phenol_c"][oxy]], oxygen, hydrogen)
                            phenols.append(
                                max(
                                    (pr.sp.hbond(oxygen, hydrogen, target), pr.sp.hbond(oxygen, flipped, target)),
                                    key=lambda hit: (hit[0], -hit[1], hit[2]),
                                )
                            )
                        if int(amine[0]) + sum(int(hit[0]) for hit in phenols) < 3:
                            continue
                        allowed = {
                            (conf["n"], asp["file_resi"], "OD1"),
                            (conf["n"], asp["file_resi"], "OD2"),
                            (oa, side["203"]["file_resi"], "OG"),
                            (ob, side["207"]["file_resi"], "OG"),
                        }
                        clearance, _atom = pr.other_clearance(cells, placed, conf["heavy"], allowed)
                        if clearance < 2.6:
                            continue
                        brute, brute_atom = pr.brute_clearance(posed, placed, conf["heavy"], allowed)
                        if brute < 2.6:
                            continue
                        positions = full_positions(mol, cid, rot, trans)
                        write_ligand(HERE / "3sn6_orciprenaline.pdb", mol, positions)
                        note = HERE / "3sn6_pose_export.txt"
                        note.write_text(
                            "\n".join(
                                [
                                    f"chi1_203 {chi['203']}",
                                    f"chi1_207 {chi['207']}",
                                    f"arrangements_tested {tested}",
                                    f"amine_A {amine[1]}",
                                    f"amine_deg {amine[2]}",
                                    f"phenol_A {' '.join(str(hit[1]) for hit in phenols)}",
                                    f"phenol_deg {' '.join(str(hit[2]) for hit in phenols)}",
                                    f"other_heavy_min {round(brute, 3)}",
                                    f"closest_other {' '.join(map(str, brute_atom))}",
                                    f"radius {radius}",
                                ]
                            )
                            + "\n",
                            encoding="utf-8",
                        )
                        print("PASS", tested, round(brute, 3), brute_atom, flush=True)
                        return
    raise SystemExit(f"no passing pose in {tested} arrangements")


if __name__ == "__main__":
    main()

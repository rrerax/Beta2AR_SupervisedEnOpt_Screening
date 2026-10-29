"""Build the 2Y03 isoprenaline membrane system.

The contact rules are already fixed in PROTOCOL.md. This script only prepares
the crystal starting point. It does not analyze a trajectory.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from openmm import LangevinMiddleIntegrator, Platform, Vec3, unit
from openmm.app import (
    ForceField,
    HBonds,
    Modeller,
    NoCutoff,
    PDBFile,
    PME,
    Simulation,
    Topology,
    element,
)
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parent
RUN = Path(r"D:\repositories\md_runs\2y03_isoprenaline")
OPM = Path(r"D:\repositories\md_runs\inputs\2y03_opm.pdb")
CRYSTAL = Path(
    r"D:\repositories\Beta2AR_SupervisedEnOpt_Screening"
    r"\followup_2026-10-16\verify_2026-10-26\crystals\2Y03.pdb"
)
RUN.mkdir(parents=True, exist_ok=True)

AMINO = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
DISULFIDES = ((114, 199), (192, 198))
CHARGE_METHOD = "openff-gnn-am1bcc-1.0.0.pt"


def as_nm(positions):
    if hasattr(positions, "value_in_unit"):
        rows = positions.value_in_unit(unit.nanometer)
        return [Vec3(float(row[0]), float(row[1]), float(row[2])) for row in rows]
    out = []
    for value in positions:
        if hasattr(value, "value_in_unit"):
            row = value.value_in_unit(unit.nanometer)
            out.append(Vec3(float(row[0]), float(row[1]), float(row[2])))
        else:
            out.append(Vec3(float(value[0]), float(value[1]), float(value[2])))
    return out


def write_pdb(topology, positions, path):
    rows = as_nm(positions)
    array = np.array([[point.x, point.y, point.z] for point in rows], dtype=float)
    with path.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(topology, array * unit.nanometer, handle, keepIds=True)


def find_atom(topology, positions, chain, resname, resi, atomname):
    for atom in topology.atoms():
        residue = atom.residue
        if (
            residue.chain.id == chain
            and residue.name == resname
            and int(residue.id) == resi
            and atom.name == atomname
        ):
            return np.array(as_nm([positions[atom.index]])[0], dtype=float)
    raise SystemExit(f"missing {chain} {resname} {resi} {atomname}")


def pair_distance(topology, positions, spec_a, spec_b):
    a = find_atom(topology, positions, *spec_a)
    b = find_atom(topology, positions, *spec_b)
    return float(np.linalg.norm(a - b) * 10.0)


def check_against_crystal(opm):
    crystal = PDBFile(str(CRYSTAL))
    pairs = [
        (("A", "5FW", 601, "OAB"), ("A", "SER", 211, "OG")),
        (("A", "5FW", 601, "OAC"), ("A", "SER", 215, "OG")),
        (("A", "5FW", 601, "NAD"), ("A", "ASP", 121, "OD1")),
        (("A", "5FW", 601, "NAD"), ("A", "ASP", 121, "OD2")),
    ]
    report = []
    for left, right in pairs:
        d_opm = pair_distance(opm.topology, opm.positions, left, right)
        d_crystal = pair_distance(crystal.topology, crystal.positions, left, right)
        report.append(
            {
                "pair": [left[3], right[1] + str(right[2]), right[3]],
                "opm_A": round(d_opm, 3),
                "crystal_A": round(d_crystal, 3),
            }
        )
        if abs(d_opm - d_crystal) > 0.02:
            raise SystemExit(f"OPM moved a crystal contact: {report[-1]}")
    return report


def write_fragments(opm):
    topology = Topology()
    chains = {
        "A": topology.addChain(id="A"),
        "C": topology.addChain(id="C"),
    }
    positions = []
    counts = {"A": 0, "C": 0}
    coords = as_nm(opm.positions)
    for residue in opm.topology.residues():
        if residue.chain.id != "A" or residue.name not in AMINO:
            continue
        resi = int(residue.id)
        if 33 <= resi <= 239:
            chain_id = "A"
        elif 278 <= resi <= 358:
            chain_id = "C"
        else:
            continue
        new_residue = topology.addResidue(residue.name, chains[chain_id], id=str(resi))
        counts[chain_id] += 1
        for atom in residue.atoms():
            if atom.element is None or atom.element.symbol == "H":
                continue
            topology.addAtom(atom.name, atom.element, new_residue)
            positions.append(coords[atom.index])
    path = RUN / "fragments_heavy.pdb"
    write_pdb(topology, positions, path)
    print("fragments", counts, flush=True)
    if counts["A"] < 200 or counts["C"] < 70:
        raise SystemExit(f"unexpected fragment sizes {counts}")
    return path


def fix_heavy_atoms(path):
    fixer = PDBFixer(filename=str(path))
    fixer.findMissingResidues()
    if fixer.missingResidues:
        raise SystemExit(f"unexpected missing residues: {fixer.missingResidues!r}")
    fixer.findMissingAtoms()
    missing = {
        f"{residue.name}{residue.id}": list(atoms)
        for residue, atoms in fixer.missingAtoms.items()
        if atoms
    }
    fixer.addMissingAtoms()
    (RUN / "missing_heavy_atoms.json").write_text(
        json.dumps(missing, indent=2), encoding="utf-8"
    )
    print(f"completed heavy atoms on {len(missing)} residues", flush=True)
    return fixer


def norm(vector):
    length = np.linalg.norm(vector)
    if length < 1e-8:
        raise SystemExit("cannot place a cap from overlapping atoms")
    return vector / length


def rotate(axis, vector, degrees):
    angle = np.deg2rad(degrees)
    axis = norm(axis)
    return (
        vector * np.cos(angle)
        + np.cross(axis, vector) * np.sin(angle)
        + axis * np.dot(axis, vector) * (1.0 - np.cos(angle))
    )


def cap_positions(residue_atoms, kind):
    xyz = {atom.name: np.array(as_nm([residue_atoms[atom]])[0]) for atom in residue_atoms}
    # residue_atoms is a dict name -> position already. Handle both.
    return xyz


def place_nme(c_pos, ca_pos, o_pos):
    axis = norm(np.cross(ca_pos - c_pos, o_pos - c_pos))
    ca_dir = norm(ca_pos - c_pos)
    candidates = (rotate(axis, ca_dir, 117.0), rotate(axis, ca_dir, -117.0))
    oxygen_dir = norm(o_pos - c_pos)
    n_dir = min(candidates, key=lambda direction: np.dot(direction, oxygen_dir))
    n_pos = c_pos + 0.133 * n_dir
    back = norm(c_pos - n_pos)
    methyl_candidates = (rotate(axis, back, 120.0), rotate(axis, back, -120.0))
    methyl_dir = min(methyl_candidates, key=lambda direction: np.dot(direction, oxygen_dir))
    return n_pos, n_pos + 0.145 * methyl_dir


def place_planar_hydrogen(center, neighbor_a, neighbor_b):
    direction = norm(-(norm(neighbor_a - center) + norm(neighbor_b - center)))
    return center + 0.101 * direction


def methyl_hydrogens(carbon, neighbor):
    bond_dir = norm(neighbor - carbon)
    helper = np.array([1.0, 0.0, 0.0]) if abs(bond_dir[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    perpendicular = norm(np.cross(bond_dir, helper))
    orthogonal = norm(np.cross(bond_dir, perpendicular))
    cosine = np.cos(np.deg2rad(109.47))
    sine = np.sin(np.deg2rad(109.47))
    placed = []
    for turn in (0.0, 120.0, 240.0):
        rotated = (
            np.cos(np.deg2rad(turn)) * perpendicular
            + np.sin(np.deg2rad(turn)) * orthogonal
        )
        placed.append(carbon + 0.109 * (cosine * bond_dir + sine * rotated))
    return placed


def place_ace(n_pos, ca_pos, c_pos):
    axis = norm(np.cross(ca_pos - n_pos, c_pos - n_pos))
    bisector = norm(norm(ca_pos - n_pos) + norm(c_pos - n_pos))
    c_ace = n_pos - 0.133 * bisector
    toward_n = norm(n_pos - c_ace)
    oxygen = c_ace + 0.123 * rotate(axis, toward_n, 120.0)
    methyl = c_ace + 0.151 * rotate(axis, toward_n, -120.0)
    return methyl, c_ace, oxygen


def atom_map(topology, positions):
    coords = as_nm(positions)
    found = {}
    for atom in topology.atoms():
        point = coords[atom.index]
        found[(atom.residue.chain.id, int(atom.residue.id), atom.name)] = (
            atom,
            np.array([point.x, point.y, point.z]),
        )
    return found


def build_capped_protein(fixer):
    old_positions = as_nm(fixer.positions)
    new_topology = Topology()
    chain_a = new_topology.addChain(id="A")
    chain_c = new_topology.addChain(id="C")
    mapping = {}
    new_positions = []

    def copy_residue(old_residue, new_chain, name):
        new_residue = new_topology.addResidue(name, new_chain, id=str(old_residue.id))
        for atom in old_residue.atoms():
            if atom.element is None or atom.element.symbol == "H":
                continue
            if atom.name == "OXT" and old_residue.chain.id == "A" and int(old_residue.id) == 239:
                continue
            mapping[atom] = new_topology.addAtom(atom.name, atom.element, new_residue)
            new_positions.append(old_positions[atom.index])
        return new_residue

    disulfide_ids = {resi for pair in DISULFIDES for resi in pair}
    for residue in fixer.topology.residues():
        resi = int(residue.id)
        if residue.chain.id == "A" and 33 <= resi <= 239:
            copy_residue(residue, chain_a, residue.name)
    located = atom_map(new_topology, new_positions)
    n_pos, methyl_pos = place_nme(
        located[("A", 239, "C")][1],
        located[("A", 239, "CA")][1],
        located[("A", 239, "O")][1],
    )
    nme = new_topology.addResidue("NME", chain_a, id="439")
    n_atom = new_topology.addAtom("N", element.nitrogen, nme)
    c_atom = new_topology.addAtom("CH3", element.carbon, nme)
    h_atom = new_topology.addAtom("H", element.hydrogen, nme)
    methyl_h_atoms = [
        new_topology.addAtom(name, element.hydrogen, nme)
        for name in ("HH31", "HH32", "HH33")
    ]
    amide_h = place_planar_hydrogen(n_pos, located[("A", 239, "C")][1], methyl_pos)
    methyl_h = methyl_hydrogens(methyl_pos, n_pos)
    new_positions.extend(
        [Vec3(*map(float, point)) for point in (n_pos, methyl_pos, amide_h, *methyl_h)]
    )
    new_topology.addBond(n_atom, c_atom)
    new_topology.addBond(n_atom, h_atom)
    for atom, _point in zip(methyl_h_atoms, methyl_h):
        new_topology.addBond(c_atom, atom)
    new_topology.addBond(located[("A", 239, "C")][0], n_atom)
    for residue in fixer.topology.residues():
        resi = int(residue.id)
        if residue.chain.id == "C" and 278 <= resi <= 358:
            copy_residue(residue, chain_c, residue.name)
    for bond in fixer.topology.bonds():
        if bond[0] in mapping and bond[1] in mapping:
            new_topology.addBond(mapping[bond[0]], mapping[bond[1]])

    located = atom_map(new_topology, new_positions)
    methyl, carbonyl, oxygen = place_ace(
        located[("C", 278, "N")][1],
        located[("C", 278, "CA")][1],
        located[("C", 278, "C")][1],
    )
    return rebuild_with_ace(
        new_topology, new_positions, methyl, carbonyl, oxygen, located
    )


def rebuild_with_ace(topology, positions, methyl, carbonyl, oxygen, located):
    rebuilt = Topology()
    chain_a = rebuilt.addChain(id="A")
    chain_c = rebuilt.addChain(id="C")
    mapping = {}
    new_positions = []
    for residue in topology.residues():
        if residue.chain.id != "A":
            continue
        copy = rebuilt.addResidue(residue.name, chain_a, id=residue.id)
        for atom in residue.atoms():
            mapping[atom] = rebuilt.addAtom(atom.name, atom.element, copy)
            new_positions.append(positions[atom.index])
    ace = rebuilt.addResidue("ACE", chain_c, id="477")
    ace_atoms = {
        "CH3": rebuilt.addAtom("CH3", element.carbon, ace),
        "C": rebuilt.addAtom("C", element.carbon, ace),
        "O": rebuilt.addAtom("O", element.oxygen, ace),
    }
    ace_h = [rebuilt.addAtom(name, element.hydrogen, ace) for name in ("HH31", "HH32", "HH33")]
    ace_methyl_h = methyl_hydrogens(methyl, carbonyl)
    new_positions.extend(
        [Vec3(*map(float, point)) for point in (methyl, carbonyl, oxygen, *ace_methyl_h)]
    )
    rebuilt.addBond(ace_atoms["CH3"], ace_atoms["C"])
    rebuilt.addBond(ace_atoms["C"], ace_atoms["O"])
    for atom in ace_h:
        rebuilt.addBond(ace_atoms["CH3"], atom)
    for residue in topology.residues():
        if residue.chain.id != "C":
            continue
        copy = rebuilt.addResidue(residue.name, chain_c, id=residue.id)
        for atom in residue.atoms():
            mapping[atom] = rebuilt.addAtom(atom.name, atom.element, copy)
            new_positions.append(positions[atom.index])
    for bond in topology.bonds():
        if bond[0] in mapping and bond[1] in mapping:
            rebuilt.addBond(mapping[bond[0]], mapping[bond[1]])
    rebuilt.addBond(ace_atoms["C"], mapping[located[("C", 278, "N")][0]])
    for first, second in DISULFIDES:
        sg1 = None
        sg2 = None
        for atom in rebuilt.atoms():
            if atom.name == "SG" and atom.residue.chain.id == "A" and int(atom.residue.id) == first:
                sg1 = atom
            if atom.name == "SG" and atom.residue.chain.id == "A" and int(atom.residue.id) == second:
                sg2 = atom
        if sg1 is None or sg2 is None:
            raise SystemExit(f"disulfide atoms missing for {first}-{second}")
        xyz1 = np.array([new_positions[sg1.index].x, new_positions[sg1.index].y, new_positions[sg1.index].z])
        xyz2 = np.array([new_positions[sg2.index].x, new_positions[sg2.index].y, new_positions[sg2.index].z])
        distance = np.linalg.norm(xyz1 - xyz2)
        if not 0.18 <= distance <= 0.23:
            raise SystemExit(f"disulfide {first}-{second} is {distance:.4f} nm, expected about 0.208")
        rebuilt.addBond(sg1, sg2)
        print(f"disulfide {first}-{second} {distance:.4f} nm", flush=True)
    return rebuilt, new_positions


def add_protein_hydrogens(topology, positions):
    forcefield = ForceField("amber14/protein.ff14SB.xml")
    modeller = Modeller(topology, positions)
    modeller.addHydrogens(forcefield, pH=7.4)
    counts = {}
    for residue in modeller.topology.residues():
        if residue.name in {"HID", "HIE", "HIP"}:
            counts[residue.name] = counts.get(residue.name, 0) + 1
        if int(residue.id) == 121 and residue.name != "ASP":
            raise SystemExit(f"Asp121 was renamed to {residue.name}")
    print("histidine protonation", counts, flush=True)
    return modeller


def ligand_molecule(opm):
    atoms = []
    for atom in opm.topology.atoms():
        residue = atom.residue
        if residue.chain.id == "A" and residue.name == "5FW" and int(residue.id) == 601:
            atoms.append(
                (
                    atom.name,
                    atom.element.symbol,
                    np.array(as_nm([opm.positions[atom.index]])[0]) * 10.0,
                )
            )
    if len(atoms) != 15:
        raise SystemExit(f"expected 15 ligand heavy atoms, found {len(atoms)}")
    editable = Chem.RWMol()
    for name, symbol, _ in atoms:
        rd_atom = Chem.Atom(symbol)
        rd_atom.SetProp("crystal_name", name)
        editable.AddAtom(rd_atom)
    radii = {"C": 0.77, "N": 0.75, "O": 0.73}
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            distance = np.linalg.norm(atoms[i][2] - atoms[j][2])
            limit = (radii[atoms[i][1]] + radii[atoms[j][1]]) * 1.35
            if distance < limit:
                editable.AddBond(i, j, Chem.BondType.SINGLE)
    crystal = editable.GetMol()
    conformer = Chem.Conformer(crystal.GetNumAtoms())
    for index, (_, _, xyz) in enumerate(atoms):
        conformer.SetAtomPosition(index, xyz.tolist())
    crystal.AddConformer(conformer, assignId=True)
    template = Chem.RemoveHs(Chem.MolFromSmiles("CC(C)[NH2+]CC(O)c1ccc(O)c(O)c1"))
    ordered = AllChem.AssignBondOrdersFromTemplate(template, crystal)
    Chem.AssignStereochemistryFrom3D(ordered, replaceExistingTags=True)
    centers = Chem.FindMolChiralCenters(ordered, includeUnassigned=True)
    if len(centers) != 1:
        raise SystemExit(f"expected one chiral center, found {centers}")
    chiral_index, cip = centers[0]
    chiral_xyz = np.array(ordered.GetConformer().GetAtomPosition(chiral_index))
    chiral_name = min(atoms, key=lambda item: np.linalg.norm(item[2] - chiral_xyz))[0]
    if chiral_name != "CAE":
        raise SystemExit(f"chiral atom mapped to {chiral_name}, expected CAE")
    charged = Chem.AddHs(ordered, addCoords=True)
    if Chem.rdMolDescriptors.CalcMolFormula(charged) != "C11H18NO3+":
        raise SystemExit(f"unexpected formula {Chem.rdMolDescriptors.CalcMolFormula(charged)}")
    heavy_names = []
    for index, atom in enumerate(charged.GetAtoms()):
        if atom.GetAtomicNum() == 1:
            continue
        xyz = np.array(charged.GetConformer().GetAtomPosition(index))
        heavy_names.append(min(atoms, key=lambda item: np.linalg.norm(item[2] - xyz))[0])
    from openff.toolkit import Molecule

    molecule = Molecule.from_rdkit(charged, allow_undefined_stereo=False, hydrogens_are_explicit=True)
    molecule.assign_partial_charges(CHARGE_METHOD)
    charge = float(np.sum(molecule.partial_charges.m_as("elementary_charge")))
    if abs(charge - 1.0) > 1e-4:
        raise SystemExit(f"ligand charge is {charge}, expected +1")
    print("ligand", molecule.to_smiles(explicit_hydrogens=False), cip, "charge", round(charge, 4), flush=True)
    return charged, heavy_names, molecule


def add_ligand(modeller, rdkit_mol, heavy_names):
    fragment = Topology()
    chain = fragment.addChain(id="L")
    residue = fragment.addResidue("5FW", chain, id="601")
    conformer = rdkit_mol.GetConformer()
    rd_atoms = []
    positions = []
    for index, atom in enumerate(rdkit_mol.GetAtoms()):
        symbol = atom.GetSymbol()
        name = heavy_names.pop(0) if symbol != "H" else f"H{index:02d}"
        rd_atoms.append(fragment.addAtom(name, element.Element.getBySymbol(symbol), residue))
        point = conformer.GetAtomPosition(index)
        positions.append(Vec3(point.x / 10.0, point.y / 10.0, point.z / 10.0))
    if heavy_names:
        raise SystemExit("ligand heavy-atom names were not consumed in order")
    for bond in rdkit_mol.GetBonds():
        fragment.addBond(rd_atoms[bond.GetBeginAtomIdx()], rd_atoms[bond.GetEndAtomIdx()])
    array = np.array([[point.x, point.y, point.z] for point in positions], dtype=float)
    modeller.add(fragment, array * unit.nanometer)


def minimize_caps_and_hydrogens(modeller, molecule):
    from openmmforcefields.generators import SMIRNOFFTemplateGenerator

    generator = SMIRNOFFTemplateGenerator(molecules=molecule, forcefield="openff-2.2.1")
    forcefield = ForceField(
        "amber14/protein.ff14SB.xml",
        "amber14/lipid17.xml",
        "amber14/tip3p.xml",
    )
    forcefield.registerTemplateGenerator(generator.generator)
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
        constraints=HBonds,
        rigidWater=False,
    )
    restraint = openmm_restraint(modeller.topology, modeller.positions, 1000.0)
    system.addForce(restraint)
    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 0.004 * unit.picoseconds
    )
    platform = Platform.getPlatformByName("OpenCL")
    platform.setPropertyDefaultValue("OpenCLPlatformIndex", "0")
    platform.setPropertyDefaultValue("DeviceIndex", "0")
    platform.setPropertyDefaultValue("Precision", "mixed")
    simulation = Simulation(
        modeller.topology,
        system,
        integrator,
        platform,
        {"OpenCLPlatformIndex": "0", "DeviceIndex": "0", "Precision": "mixed"},
    )
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy(maxIterations=300)
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    print(f"vacuum minimization {energy:.1f} kJ/mol", flush=True)
    if not np.isfinite(energy):
        raise SystemExit("vacuum minimization energy is not finite")
    before = np.array([[point.x, point.y, point.z] for point in as_nm(modeller.positions)])
    after = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    moved = np.linalg.norm(after - before, axis=1)
    worst_index = max(
        (
            atom.index
            for atom in modeller.topology.atoms()
            if atom.element is not None
            and atom.element.symbol != "H"
            and atom.residue.name not in {"ACE", "NME"}
            and (atom.name == "CA" or atom.residue.name == "5FW")
        ),
        key=lambda index: moved[index],
    )
    worst_atom = list(modeller.topology.atoms())[worst_index]
    worst = float(moved[worst_index])
    print(
        f"largest restrained move {worst:.4f} nm at {worst_atom.residue.name} "
        f"{worst_atom.residue.id} {worst_atom.name}",
        flush=True,
    )
    if worst > 0.10:
        raise SystemExit("crystal heavy atoms moved by more than 1 Å during hydrogen minimization")
    modeller.positions = state.getPositions()
    return forcefield


def openmm_restraint(topology, positions, k_value):
    from openmm import CustomExternalForce

    force = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addGlobalParameter("k", k_value)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    coords = as_nm(positions)
    for atom in topology.atoms():
        if atom.element is None or atom.element.symbol == "H":
            continue
        if atom.residue.name in {"ACE", "NME"}:
            continue
        xyz = coords[atom.index]
        force.addParticle(atom.index, [xyz.x, xyz.y, xyz.z])
    return force


def membrane(modeller, forcefield):
    platform = Platform.getPlatformByName("OpenCL")
    platform.setPropertyDefaultValue("OpenCLPlatformIndex", "0")
    platform.setPropertyDefaultValue("DeviceIndex", "0")
    platform.setPropertyDefaultValue("Precision", "mixed")
    print("adding POPC membrane", flush=True)
    modeller.addMembrane(
        forcefield,
        lipidType="POPC",
        membraneCenterZ=0 * unit.nanometer,
        minimumPadding=1.0 * unit.nanometer,
        ionicStrength=0.15 * unit.molar,
        neutralize=True,
        platform=platform,
    )
    print(
        "atoms",
        modeller.topology.getNumAtoms(),
        "residues",
        modeller.topology.getNumResidues(),
        flush=True,
    )


def system_charge(system):
    total = 0.0
    for force in system.getForces():
        if force.getName() == "NonbondedForce":
            for index in range(force.getNumParticles()):
                total += force.getParticleParameters(index)[0].value_in_unit(
                    unit.elementary_charge
                )
    return total


def assert_ligand_hydrogens(topology):
    wanted = {}
    for atom in topology.atoms():
        if atom.residue.name == "5FW" and atom.name in {"NAD", "OAB", "OAC", "OAA"}:
            wanted[atom.name] = atom
    counts = {name: 0 for name in wanted}
    for bond in topology.bonds():
        first, second = bond[0], bond[1]
        for name, atom in wanted.items():
            if atom not in (first, second):
                continue
            other = second if first is atom else first
            if other.element is not None and other.element.symbol == "H":
                counts[name] += 1
    print("ligand hydrogens", counts, flush=True)
    expected = {"NAD": 2, "OAB": 1, "OAC": 1, "OAA": 1}
    if counts != expected:
        raise SystemExit(f"ligand protonation differs from the crystal rule: {counts}")


def write_selections(topology):
    wanted = {
        "asp121_od": [],
        "ser211_og": [],
        "ser212_og": [],
        "ser215_og": [],
        "ligand_heavy": [],
        "ligand_atoms": {},
        "amine_n": [],
        "phenol_oab": [],
        "phenol_oac": [],
        "benzylic_oaa": [],
        "protein_ca": [],
        "tm_core_ca": [],
        "water_o": [],
    }
    water_names = set()
    for atom in topology.atoms():
        residue = atom.residue
        resi = int(residue.id) if str(residue.id).isdigit() else None
        if residue.name == "5FW":
            if atom.element.symbol != "H":
                wanted["ligand_heavy"].append(atom.index)
                wanted["ligand_atoms"][atom.name] = atom.index
            if atom.name == "NAD":
                wanted["amine_n"].append(atom.index)
            elif atom.name == "OAB":
                wanted["phenol_oab"].append(atom.index)
            elif atom.name == "OAC":
                wanted["phenol_oac"].append(atom.index)
            elif atom.name == "OAA":
                wanted["benzylic_oaa"].append(atom.index)
        if atom.name == "CA" and residue.name in AMINO:
            wanted["protein_ca"].append(atom.index)
            if resi is not None and (43 <= resi <= 229 or 288 <= resi <= 348):
                wanted["tm_core_ca"].append(atom.index)
        if residue.name == "ASP" and resi == 121 and atom.name in {"OD1", "OD2"}:
            wanted["asp121_od"].append(atom.index)
        if residue.name == "SER" and atom.name == "OG" and resi in {211, 212, 215}:
            wanted[{"211": "ser211_og", "212": "ser212_og", "215": "ser215_og"}[str(resi)]].append(
                atom.index
            )
        if atom.name == "O" and residue.name in {"HOH", "WAT", "TIP3"}:
            water_names.add(residue.name)
            wanted["water_o"].append(atom.index)
    required = [
        "amine_n", "phenol_oab", "phenol_oac", "benzylic_oaa", "asp121_od",
        "ser211_og", "ser215_og", "ligand_heavy", "tm_core_ca", "water_o",
    ]
    for key in required:
        if not wanted[key]:
            raise SystemExit(f"selection {key} is empty")
    wanted["water_residue_names"] = sorted(water_names)
    (RUN / "selections.json").write_text(json.dumps(wanted), encoding="utf-8")
    print("waters", len(wanted["water_o"]), sorted(water_names), flush=True)


def main():
    opm = PDBFile(str(OPM))
    distances = check_against_crystal(opm)
    fragment_path = write_fragments(opm)
    fixer = fix_heavy_atoms(fragment_path)
    topology, positions = build_capped_protein(fixer)
    modeller = add_protein_hydrogens(topology, positions)
    rdkit_mol, heavy_names, molecule = ligand_molecule(opm)
    add_ligand(modeller, rdkit_mol, heavy_names.copy())
    forcefield = minimize_caps_and_hydrogens(modeller, molecule)
    assert_ligand_hydrogens(modeller.topology)
    membrane(modeller, forcefield)
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
        hydrogenMass=1.5 * unit.amu,
        removeCMMotion=False,
    )
    charge = system_charge(system)
    print(f"net charge {charge:.4f}", flush=True)
    if abs(charge) > 0.05:
        raise SystemExit("membrane system is not neutral")
    write_pdb(modeller.topology, modeller.positions, RUN / "solvated.pdb")
    from openmm import XmlSerializer

    (RUN / "system.xml").write_text(XmlSerializer.serialize(system), encoding="utf-8")
    write_selections(modeller.topology)
    info = {
        "distances_A": distances,
        "ligand_smiles": molecule.to_smiles(explicit_hydrogens=False),
        "charge_method": CHARGE_METHOD,
        "ligand_forcefield": "openff-2.2.1",
        "protein_forcefield": "amber14/protein.ff14SB.xml",
        "lipid_forcefield": "amber14/lipid17.xml",
        "water": "amber14/tip3p.xml",
        "atoms": modeller.topology.getNumAtoms(),
        "net_charge": charge,
    }
    (RUN / "system_info.json").write_text(json.dumps(info, indent=2), encoding="utf-8")
    print("wrote", RUN, flush=True)


if __name__ == "__main__":
    main()

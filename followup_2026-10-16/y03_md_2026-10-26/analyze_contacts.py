"""Contact occupancies for the three 2Y03 production replicates.

The result figure is written only after every replicate has reached the
production length in run_plan.json. Rules are in PROTOCOL.md.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


def _determinant(values):
    # numpy.linalg aborts in this environment. mdtraj only needs 3x3 determinants.
    values = np.asarray(values, dtype=float)
    return (
        values[..., 0, 0] * (values[..., 1, 1] * values[..., 2, 2] - values[..., 1, 2] * values[..., 2, 1])
        - values[..., 0, 1] * (values[..., 1, 0] * values[..., 2, 2] - values[..., 1, 2] * values[..., 2, 0])
        + values[..., 0, 2] * (values[..., 1, 0] * values[..., 2, 1] - values[..., 1, 1] * values[..., 2, 0])
    )


np.linalg.det = _determinant

import mdtraj as md
from PIL import Image, ImageDraw, ImageFont

RUN = Path(r"D:\repositories\md_runs\2y03_isoprenaline")
OUT = Path(__file__).resolve().parent
SEEDS = (20261026, 20261027, 20261028)
HBOND_NM = 0.25
ANGLE = 120.0
WATER_NM = 0.32
OCCUPANCY = 0.20


def bonded_hydrogens(topology, index):
    index = int(index)
    hydrogens = []
    for bond in topology.bonds:
        atoms = (bond.atom1, bond.atom2)
        if index not in (atoms[0].index, atoms[1].index):
            continue
        other = atoms[1] if atoms[0].index == index else atoms[0]
        if other.element.symbol == "H":
            hydrogens.append(other.index)
    return hydrogens


def angle_degrees(donor, hydrogen, acceptor):
    first = donor - hydrogen
    second = acceptor - hydrogen
    cosine = np.sum(first * second, axis=-1) / (
        np.linalg.norm(first, axis=-1) * np.linalg.norm(second, axis=-1)
    )
    return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))


def hbond_frames(xyz, donor, hydrogens, acceptors):
    passed = np.zeros(xyz.shape[0], dtype=bool)
    for hydrogen in hydrogens:
        for acceptor in acceptors:
            distance = np.linalg.norm(xyz[:, hydrogen] - xyz[:, acceptor], axis=1)
            angle = angle_degrees(xyz[:, donor], xyz[:, hydrogen], xyz[:, acceptor])
            passed |= (distance <= HBOND_NM) & (angle >= ANGLE)
    return passed


def water_frames(xyz, waters, ligand_atom, acceptors):
    passed = np.zeros(xyz.shape[0], dtype=bool)
    water_xyz = xyz[:, waters, :]
    ligand_xyz = xyz[:, ligand_atom, :]
    ligand_distance = np.linalg.norm(water_xyz - ligand_xyz[:, None, :], axis=2)
    for acceptor in acceptors:
        acceptor_xyz = xyz[:, acceptor, :]
        acceptor_distance = np.linalg.norm(water_xyz - acceptor_xyz[:, None, :], axis=2)
        passed |= np.any(
            (ligand_distance <= WATER_NM) & (acceptor_distance <= WATER_NM),
            axis=1,
        )
    return passed


def _det3(matrix):
    return (
        matrix[0, 0] * (matrix[1, 1] * matrix[2, 2] - matrix[1, 2] * matrix[2, 1])
        - matrix[0, 1] * (matrix[1, 0] * matrix[2, 2] - matrix[1, 2] * matrix[2, 0])
        + matrix[0, 2] * (matrix[1, 0] * matrix[2, 1] - matrix[1, 1] * matrix[2, 0])
    )


def _jacobi_eigen(matrix):
    values = np.array(matrix, dtype=float)
    vectors = np.eye(3)
    for _ in range(16):
        for first, second in ((0, 1), (0, 2), (1, 2)):
            if abs(values[first, second]) < 1e-15:
                continue
            tau = (values[second, second] - values[first, first]) / (2.0 * values[first, second])
            tangent = math.copysign(1.0, tau) / (abs(tau) + math.sqrt(1.0 + tau * tau))
            cosine = 1.0 / math.sqrt(1.0 + tangent * tangent)
            sine = tangent * cosine
            app = values[first, first]
            aqq = values[second, second]
            apq = values[first, second]
            values[first, first] = cosine * cosine * app - 2 * sine * cosine * apq + sine * sine * aqq
            values[second, second] = sine * sine * app + 2 * sine * cosine * apq + cosine * cosine * aqq
            values[first, second] = values[second, first] = 0.0
            for other in range(3):
                if other == first or other == second:
                    continue
                arp = values[other, first]
                arq = values[other, second]
                values[other, first] = values[first, other] = cosine * arp - sine * arq
                values[other, second] = values[second, other] = sine * arp + cosine * arq
            for other in range(3):
                vrp = vectors[other, first]
                vrq = vectors[other, second]
                vectors[other, first] = cosine * vrp - sine * vrq
                vectors[other, second] = sine * vrp + cosine * vrq
    return np.array([values[0, 0], values[1, 1], values[2, 2]]), vectors


def _kabsch_rotation(mobile, target):
    """Rotation taking centered mobile onto centered target. Avoids np.linalg."""
    covariance = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            covariance[i, j] = float(np.sum(mobile[:, i] * target[:, j]))
    gram = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            gram[i, j] = sum(covariance[k, i] * covariance[k, j] for k in range(3))
    evals, right = _jacobi_eigen(gram)
    order = np.argsort(evals)[::-1]
    right = right[:, order]
    singular = np.sqrt(np.maximum(evals[order], 0.0))
    left = np.zeros((3, 3))
    for column in range(3):
        if singular[column] < 1e-12:
            continue
        for i in range(3):
            left[i, column] = sum(covariance[i, n] * right[n, column] for n in range(3)) / singular[column]
    for column in range(3):
        if float(np.sum(left[:, column] * left[:, column])) > 0.5:
            continue
        others = [n for n in range(3) if n != column]
        a_col = left[:, others[0]]
        b_col = left[:, others[1]]
        left[:, column] = np.array(
            [
                a_col[1] * b_col[2] - a_col[2] * b_col[1],
                a_col[2] * b_col[0] - a_col[0] * b_col[2],
                a_col[0] * b_col[1] - a_col[1] * b_col[0],
            ]
        )

    def rotation(sign):
        matrix = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                matrix[i, j] = left[i, 0] * right[j, 0] + left[i, 1] * right[j, 1] + sign * left[i, 2] * right[j, 2]
        return matrix

    rot = rotation(1.0)
    if _det3(rot) < 0.0:
        rot = rotation(-1.0)
    return rot


def _minimum_image(points, anchor, lengths):
    out = np.array(points, dtype=float, copy=True)
    for axis in range(3):
        delta = out[:, axis] - float(anchor[axis])
        out[:, axis] = out[:, axis] - float(lengths[axis]) * np.round(delta / float(lengths[axis]))
    return out


def ligand_rmsd(trajectory, reference, selections):
    """Ligand heavy-atom RMSD in Å after imaging onto the ligand and aligning the TM core.

    Stored frames wrap molecules into one periodic cell, so a raw CA fit mixes
    two images and reports a box-length displacement.
    """
    ca = np.array(selections["tm_core_ca"], dtype=int)
    heavy = np.array(selections["ligand_heavy"], dtype=int)
    ref_lengths = reference.unitcell_lengths[0]
    ref_anchor = np.mean(reference.xyz[0, heavy, :], axis=0)
    target_ca = _minimum_image(reference.xyz[0, ca, :], ref_anchor, ref_lengths)
    reference_heavy = _minimum_image(reference.xyz[0, heavy, :], ref_anchor, ref_lengths)
    target_mean = np.mean(target_ca, axis=0)
    target_centered = target_ca - target_mean
    values = np.empty(trajectory.n_frames, dtype=float)
    for frame in range(trajectory.n_frames):
        anchor = np.mean(trajectory.xyz[frame, heavy, :], axis=0)
        lengths = trajectory.unitcell_lengths[frame]
        mobile = _minimum_image(trajectory.xyz[frame, ca, :], anchor, lengths)
        ligand = _minimum_image(trajectory.xyz[frame, heavy, :], anchor, lengths)
        mobile_mean = np.mean(mobile, axis=0)
        rot = _kabsch_rotation(mobile - mobile_mean, target_centered)
        ligand = ligand - mobile_mean
        aligned = np.column_stack(
            [ligand[:, 0] * rot[0, i] + ligand[:, 1] * rot[1, i] + ligand[:, 2] * rot[2, i] for i in range(3)]
        )
        aligned = aligned + target_mean
        delta = aligned - reference_heavy
        values[frame] = np.sqrt(np.mean(np.sum(delta * delta, axis=1))) * 10.0
    return values


def summarize(trajectory, reference, selections):
    topology = reference.topology
    xyz = trajectory.xyz
    amine = selections["amine_n"][0]
    oab = selections["phenol_oab"][0]
    oac = selections["phenol_oac"][0]
    contacts = {
        "amine_asp121": hbond_frames(
            xyz, amine, bonded_hydrogens(topology, amine), selections["asp121_od"]
        ),
        "oab_ser211": hbond_frames(
            xyz, oab, bonded_hydrogens(topology, oab), selections["ser211_og"]
        ),
        "oac_ser215": hbond_frames(
            xyz, oac, bonded_hydrogens(topology, oac), selections["ser215_og"]
        ),
    }
    waters = np.array(selections["water_o"], dtype=int)
    bridges = {
        "water_amine_asp121": water_frames(xyz, waters, amine, selections["asp121_od"]),
        "water_oab_ser211": water_frames(xyz, waters, oab, selections["ser211_og"]),
        "water_oac_ser215": water_frames(xyz, waters, oac, selections["ser215_og"]),
    }
    oaa = selections["benzylic_oaa"][0]
    od1 = selections["asp121_od"][0]
    ser212 = selections["ser212_og"][0]
    lengths = trajectory.unitcell_lengths

    def median_distance_a(left, right):
        delta = xyz[:, left, :] - xyz[:, right, :]
        for axis in range(3):
            delta[:, axis] -= lengths[:, axis] * np.round(delta[:, axis] / lengths[:, axis])
        return float(np.median(np.sqrt(np.sum(delta * delta, axis=1))) * 10.0)

    descriptive = {
        "oaa_od1_median_A": median_distance_a(oaa, od1),
        "oab_ser212_median_A": median_distance_a(oab, ser212),
        "oac_ser212_median_A": median_distance_a(oac, ser212),
    }
    occupancies = {
        name: float(np.mean(frames))
        for name, frames in {**contacts, **bridges}.items()
    }
    return occupancies, descriptive, ligand_rmsd(trajectory, reference, selections)


def write_figure(rows, rmsd):
    labels = [
        "amine\nAsp121",
        "OAB\nSer211",
        "OAC\nSer215",
        "water\namine",
        "water\nOAB",
        "water\nOAC",
    ]
    keys = [
        "amine_asp121",
        "oab_ser211",
        "oac_ser215",
        "water_amine_asp121",
        "water_oab_ser211",
        "water_oac_ser215",
    ]
    colors = [(0, 114, 178), (213, 94, 0), (0, 158, 115)]
    font_small = ImageFont.truetype(r"C:\Windows\Fonts\arial.ttf", 18)
    font_label = ImageFont.truetype(r"C:\Windows\Fonts\arial.ttf", 24)
    image = Image.new("RGB", (1600, 780), "white")
    draw = ImageDraw.Draw(image)
    left, top, right, bottom = 90, 70, 760, 520
    draw.line((left, top, left, bottom), fill=(20, 20, 20), width=2)
    draw.line((left, bottom, right, bottom), fill=(20, 20, 20), width=2)
    draw.text((left, 36), "production occupancy", fill=(20, 20, 20), font=font_label)
    for tick in (0.0, 0.5, 1.0):
        yy = bottom - tick * (bottom - top)
        draw.line((left - 8, yy, left, yy), fill=(20, 20, 20), width=2)
        draw.text((left - 12, yy), f"{tick:.1f}", fill=(20, 20, 20), font=font_small, anchor="rm")
    line_y = bottom - OCCUPANCY * (bottom - top)
    draw.line((left, line_y, right, line_y), fill=(120, 120, 120), width=2)
    slot = (right - left) / len(keys)
    bar = slot * 0.22
    for offset, seed in enumerate(SEEDS):
        color = colors[offset]
        for index, key in enumerate(keys):
            value = rows[str(seed)]["occupancy"][key]
            x0 = left + index * slot + slot * 0.15 + offset * bar
            y0 = bottom - value * (bottom - top)
            draw.rectangle((x0, y0, x0 + bar - 2, bottom), fill=color)
        draw.rectangle((90 + offset * 220, 700, 130 + offset * 220, 714), fill=color)
        draw.text((138 + offset * 220, 706), str(seed), fill=(20, 20, 20), font=font_small, anchor="lm")
    for index, label in enumerate(labels):
        xx = left + index * slot + slot * 0.45
        draw.multiline_text((xx, bottom + 8), label, fill=(20, 20, 20), font=font_small, anchor="ma", align="center", spacing=2)
    left2, right2 = 980, 1520
    draw.line((left2, top, left2, bottom), fill=(20, 20, 20), width=2)
    draw.line((left2, bottom, right2, bottom), fill=(20, 20, 20), width=2)
    draw.text((left2, 36), "ligand RMSD (A)", fill=(20, 20, 20), font=font_label)
    longest = max(len(rmsd[str(seed)]) for seed in SEEDS)
    peak = max(float(np.max(rmsd[str(seed)])) for seed in SEEDS)
    peak = max(peak, 0.1)
    for tick in (0.0, peak / 2.0, peak):
        yy = bottom - tick / peak * (bottom - top)
        draw.line((left2 - 8, yy, left2, yy), fill=(20, 20, 20), width=2)
        draw.text((left2 - 12, yy), f"{tick:.1f}", fill=(20, 20, 20), font=font_small, anchor="rm")
    for seed, color in zip(SEEDS, colors):
        series = rmsd[str(seed)]
        points = []
        stride = 50
        for frame in range(0, len(series), stride):
            value = series[frame]
            xx = left2 + frame / max(longest - 1, 1) * (right2 - left2)
            yy = bottom - float(value) / peak * (bottom - top)
            points.append((xx, yy))
        if len(points) >= 2:
            draw.line(points, fill=color, width=2)
    for tick_ns in range(0, int(round(longest * 0.01)) + 1, max(1, int(round(longest * 0.01)) // 4 or 1)):
        xx = left2 + (tick_ns / 0.01) / max(longest - 1, 1) * (right2 - left2)
        if left2 <= xx <= right2:
            draw.text((xx, bottom + 10), str(tick_ns), fill=(20, 20, 20), font=font_small, anchor="mt")
    draw.text(
        ((left2 + right2) / 2, bottom + 40),
        "production time (ns)",
        fill=(20, 20, 20),
        font=font_label,
        anchor="mt",
    )
    image.save(OUT / "contacts_three_replicates.png", dpi=(200, 200))
    image.save(RUN / "contacts_three_replicates.png", dpi=(200, 200))


def main():
    plan = json.loads((RUN / "run_plan.json").read_text(encoding="utf-8"))
    selections = json.loads((RUN / "selections.json").read_text(encoding="utf-8"))
    expected = int(round(plan["production_ns"] * 1000.0 / DT_PS)) if False else int(
        round(plan["production_ns"] * 1000.0 / 0.004 / 2500)
    )
    reference = md.load(str(RUN / "solvated.pdb"))
    rows = {}
    rmsd = {}
    for seed in SEEDS:
        path = RUN / f"rep_{seed}" / "production.dcd"
        progress = RUN / f"rep_{seed}" / "progress.json"
        if not path.exists() or not progress.exists():
            raise SystemExit(f"{seed} has no production trajectory")
        if json.loads(progress.read_text(encoding="utf-8")).get("phase") != "done":
            raise SystemExit(f"{seed} production is not finished")
        trajectory = md.load(str(path), top=reference.topology)
        if trajectory.n_frames < int(expected * 0.95):
            raise SystemExit(f"{seed} has {trajectory.n_frames} frames, expected about {expected}")
        occupancies, descriptive, ligand = summarize(trajectory, reference, selections)
        rows[str(seed)] = {"occupancy": occupancies, "descriptive": descriptive}
        rmsd[str(seed)] = [float(value) for value in ligand]
        print(f"replicate {seed} frames {trajectory.n_frames} occupancy {occupancies}", flush=True)
        del trajectory
    calls = {}
    keys = rows[str(SEEDS[0])]["occupancy"]
    for key in keys:
        replicates = [rows[str(seed)]["occupancy"][key] >= OCCUPANCY for seed in SEEDS]
        calls[key] = {
            "occupancy": {str(seed): rows[str(seed)]["occupancy"][key] for seed in SEEDS},
            "persists": sum(replicates) >= 2,
        }
    result = {
        "calls": calls,
        "descriptive": {seed: rows[seed]["descriptive"] for seed in rows},
        "ligand_rmsd_A": {
            str(seed): {
                "median": float(np.median(rmsd[str(seed)])),
                "max": float(np.max(rmsd[str(seed)])),
            }
            for seed in SEEDS
        },
    }
    (RUN / "contact_result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    (OUT / "contact_result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    write_figure(rows, rmsd)
    print(json.dumps(calls, indent=2))


if __name__ == "__main__":
    main()

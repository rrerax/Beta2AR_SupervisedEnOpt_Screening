"""Equilibrate and run three 2Y03 replicates.

Production length is chosen from a speed benchmark and written to run_plan.json
before any contact is calculated. Contact rules live in PROTOCOL.md.
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
from openmm import (
    CustomExternalForce,
    LangevinMiddleIntegrator,
    MonteCarloMembraneBarostat,
    Platform,
    XmlSerializer,
    unit,
)
from openmm.app import DCDReporter, PDBFile, Simulation, StateDataReporter

RUN = Path(r"D:\repositories\md_runs\2y03_isoprenaline")
SEEDS = (20261026, 20261027, 20261028)
DT_PS = 0.004
STAGES = ((125, 500.0), (125, 500.0), (250, 100.0), (250, 25.0), (500, 0.0))
REPORT_STEPS = 2500  # 10 ps
CHECKPOINT_STEPS = 250_000  # 1 ns


def steps_for_ps(picoseconds):
    return int(round(picoseconds / DT_PS))


def load_inputs():
    pdb = PDBFile(str(RUN / "solvated.pdb"))
    system = XmlSerializer.deserialize((RUN / "system.xml").read_text(encoding="utf-8"))
    selections = json.loads((RUN / "selections.json").read_text(encoding="utf-8"))
    return pdb, system, selections


def coordinates(positions):
    if hasattr(positions, "value_in_unit"):
        return positions.value_in_unit(unit.nanometer)
    rows = []
    for value in positions:
        if hasattr(value, "value_in_unit"):
            row = value.value_in_unit(unit.nanometer)
            rows.append((float(row[0]), float(row[1]), float(row[2])))
        else:
            rows.append((float(value[0]), float(value[1]), float(value[2])))
    return rows


def add_restraint(system, topology, positions, selections):
    # OpenCL folds each molecule into one periodic cell. A raw x-x0 spring then
    # pulls by a full box length. The minimum image uses the current box size.
    force = CustomExternalForce(
        "0.5*k*(dx*dx+dy*dy+dz*dz); "
        "dx=x-x0-Lx*floor((x-x0)/Lx+0.5); "
        "dy=y-y0-Ly*floor((y-y0)/Ly+0.5); "
        "dz=z-z0-Lz*floor((z-z0)/Lz+0.5)"
    )
    force.addGlobalParameter("k", 0.0)
    for name in ("Lx", "Ly", "Lz"):
        force.addGlobalParameter(name, 10.0)
    for name in ("x0", "y0", "z0"):
        force.addPerParticleParameter(name)
    coords = coordinates(positions)
    keep = set(selections["protein_ca"]) | set(selections["ligand_heavy"])
    for index in sorted(keep):
        x_coord, y_coord, z_coord = coords[index]
        force.addParticle(int(index), [x_coord, y_coord, z_coord])
    system.addForce(force)
    return force


def set_box_lengths(simulation):
    axis_x, axis_y, axis_z = simulation.context.getState().getPeriodicBoxVectors()
    simulation.context.setParameter("Lx", float(axis_x.x))
    simulation.context.setParameter("Ly", float(axis_y.y))
    simulation.context.setParameter("Lz", float(axis_z.z))


def add_barostat(system, seed):
    barostat = MonteCarloMembraneBarostat(
        1.0 * unit.bar,
        0.0 * unit.bar * unit.nanometer,
        300.0 * unit.kelvin,
        MonteCarloMembraneBarostat.XYIsotropic,
        MonteCarloMembraneBarostat.ZFree,
    )
    barostat.setFrequency(100)
    barostat.setRandomNumberSeed(int(seed))
    system.addForce(barostat)


def make_simulation(topology, system, seed):
    integrator = LangevinMiddleIntegrator(
        300.0 * unit.kelvin, 1.0 / unit.picosecond, DT_PS * unit.picoseconds
    )
    integrator.setRandomNumberSeed(int(seed))
    platform = Platform.getPlatformByName("OpenCL")
    platform.setPropertyDefaultValue("OpenCLPlatformIndex", "0")
    platform.setPropertyDefaultValue("DeviceIndex", "0")
    platform.setPropertyDefaultValue("Precision", "mixed")
    return Simulation(
        topology,
        system,
        integrator,
        platform,
        {"OpenCLPlatformIndex": "0", "DeviceIndex": "0", "Precision": "mixed"},
    )


def benchmark(topology, system, positions):
    simulation = make_simulation(topology, system, SEEDS[0])
    simulation.context.setPositions(positions)
    set_box_lengths(simulation)
    simulation.context.setParameter("k", 0.0)
    simulation.minimizeEnergy(maxIterations=200)
    simulation.context.setVelocitiesToTemperature(300.0 * unit.kelvin, SEEDS[0])
    simulation.step(500)
    started = time.perf_counter()
    simulation.step(5000)
    elapsed = time.perf_counter() - started
    del simulation
    nanoseconds = 5000 * DT_PS / 1000.0
    per_day = nanoseconds / elapsed * 86400.0
    if per_day >= 40.0:
        length = 50.0
    elif per_day >= 15.0:
        length = 20.0
    else:
        length = 10.0
    plan = {
        "benchmark_ns_per_day": per_day,
        "benchmark_seconds_for_5000_steps": elapsed,
        "production_ns": length,
        "timestep_ps": DT_PS,
        "seeds": list(SEEDS),
    }
    (RUN / "run_plan.json").write_text(json.dumps(plan, indent=2), encoding="utf-8")
    print("run plan", plan, flush=True)
    return plan


def determinant(matrix):
    return (
        matrix[0, 0] * (matrix[1, 1] * matrix[2, 2] - matrix[1, 2] * matrix[2, 1])
        - matrix[0, 1] * (matrix[1, 0] * matrix[2, 2] - matrix[1, 2] * matrix[2, 0])
        + matrix[0, 2] * (matrix[1, 0] * matrix[2, 1] - matrix[1, 1] * matrix[2, 0])
    )


def jacobi_eigen(matrix):
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


def kabsch_rmsd(mobile, target):
    """RMSD in the same units as the coordinates. Avoids numpy.linalg: its
    matrix kernels abort on this machine."""
    mobile = np.asarray(mobile, dtype=float)
    target = np.asarray(target, dtype=float)
    mobile = mobile - np.mean(mobile, axis=0)
    target = target - np.mean(target, axis=0)
    covariance = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            covariance[i, j] = float(np.sum(mobile[:, i] * target[:, j]))
    gram = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            gram[i, j] = sum(covariance[k, i] * covariance[k, j] for k in range(3))
    evals, right = jacobi_eigen(gram)
    order = np.argsort(evals)[::-1]
    right = right[:, order]
    singular = np.sqrt(np.maximum(evals[order], 0.0))
    left = np.zeros((3, 3))
    for column in range(3):
        if singular[column] < 1e-12:
            continue
        for i in range(3):
            left[i, column] = (
                sum(covariance[i, n] * right[n, column] for n in range(3)) / singular[column]
            )
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
                matrix[i, j] = (
                    left[i, 0] * right[j, 0]
                    + left[i, 1] * right[j, 1]
                    + sign * left[i, 2] * right[j, 2]
                )
        return matrix

    rot = rotation(1.0)
    if determinant(rot) < 0.0:
        rot = rotation(-1.0)
    aligned = np.column_stack(
        [
            mobile[:, 0] * rot[0, i] + mobile[:, 1] * rot[1, i] + mobile[:, 2] * rot[2, i]
            for i in range(3)
        ]
    )
    delta = aligned - target
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


def position_array(positions):
    coords = coordinates(positions)
    if isinstance(coords, np.ndarray) and coords.dtype != object:
        return np.asarray(coords, dtype=float)
    return np.asarray([(float(value[0]), float(value[1]), float(value[2])) for value in coords])


def tm_rmsd(topology, reference, current, selections):
    indices = np.array(selections["tm_core_ca"], dtype=int)
    reference_positions = position_array(PDBFile(str(reference)).positions)
    return kabsch_rmsd(np.asarray(current, dtype=float)[indices], reference_positions[indices])


def run_chunks(simulation, nsteps, checkpoint, progress_path, progress):
    done = int(progress.get("step", 0))
    while done < nsteps:
        chunk = min(CHECKPOINT_STEPS, nsteps - done)
        taken = 0
        while taken < chunk:
            set_box_lengths(simulation)
            step = min(10_000, chunk - taken)
            simulation.step(step)
            taken += step
        done += chunk
        simulation.saveCheckpoint(str(checkpoint))
        progress["step"] = done
        progress_path.write_text(json.dumps(progress), encoding="utf-8")
        energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        value = energy.value_in_unit(unit.kilojoule_per_mole)
        print(f"{progress['phase']} {progress.get('stage')} {done}/{nsteps} {value:.0f}", flush=True)
        if not np.isfinite(value):
            raise SystemExit("energy became non-finite")


def run_replicate(topology, base_system, positions, selections, seed, plan):
    folder = RUN / f"rep_{seed}"
    folder.mkdir(exist_ok=True)
    progress_path = folder / "progress.json"
    checkpoint = folder / "checkpoint.chk"
    progress = json.loads(progress_path.read_text(encoding="utf-8")) if progress_path.exists() else {}
    system = XmlSerializer.deserialize(XmlSerializer.serialize(base_system))
    add_restraint(system, topology, positions, selections)
    if progress.get("phase") in {"equil", "prod"}:
        add_barostat(system, seed)
    simulation = make_simulation(topology, system, seed)
    if checkpoint.exists() and progress.get("phase") in {"equil", "prod"}:
        simulation.loadCheckpoint(str(checkpoint))
    else:
        simulation.context.setPositions(positions)
        set_box_lengths(simulation)
        simulation.context.setParameter("k", STAGES[0][1])
        simulation.minimizeEnergy(maxIterations=1000)
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if not np.isfinite(energy):
            raise SystemExit(f"replicate {seed} minimization energy is not finite")
        minimized = state.getPositions()
        box = state.getPeriodicBoxVectors()
        print(f"replicate {seed} minimized {energy:.0f} kJ/mol", flush=True)
        del simulation
        add_barostat(system, seed)
        simulation = make_simulation(topology, system, seed)
        simulation.context.setPeriodicBoxVectors(*box)
        simulation.context.setPositions(minimized)
        set_box_lengths(simulation)
        simulation.context.setParameter("k", STAGES[0][1])
        simulation.context.setVelocitiesToTemperature(300.0 * unit.kelvin, int(seed))
        progress = {"phase": "equil", "stage": 0, "step": 0}
        progress_path.write_text(json.dumps(progress), encoding="utf-8")

    if progress["phase"] == "equil":
        for stage_index, (picoseconds, k_value) in enumerate(STAGES):
            if stage_index < int(progress["stage"]):
                continue
            if stage_index != int(progress["stage"]):
                progress["stage"] = stage_index
                progress["step"] = 0
            simulation.context.setParameter("k", k_value)
            if not any(isinstance(reporter, DCDReporter) for reporter in simulation.reporters):
                simulation.reporters.append(
                    DCDReporter(str(folder / "equilibration.dcd"), REPORT_STEPS, append=checkpoint.exists())
                )
                simulation.reporters.append(
                    StateDataReporter(
                        str(folder / "equilibration.log"),
                        REPORT_STEPS,
                        step=True,
                        temperature=True,
                        potentialEnergy=True,
                        speed=True,
                        append=checkpoint.exists(),
                    )
                )
            run_chunks(simulation, steps_for_ps(picoseconds), checkpoint, progress_path, progress)
        current_nm = position_array(simulation.context.getState(getPositions=True).getPositions())
        value = tm_rmsd(topology, RUN / "solvated.pdb", current_nm, selections)
        (folder / "equil_tm_ca_rmsd.json").write_text(
            json.dumps({"nm": value}), encoding="utf-8"
        )
        print(f"replicate {seed} TM CA RMSD {value:.3f} nm", flush=True)
        if value > 0.30:
            raise SystemExit(f"replicate {seed} TM core moved {value:.3f} nm during equilibration")
        progress = {"phase": "prod", "stage": 0, "step": 0}
        progress_path.write_text(json.dumps(progress), encoding="utf-8")
        simulation.reporters.clear()

    simulation.context.setParameter("k", 0.0)
    if progress["phase"] == "prod":
        if not any(isinstance(reporter, DCDReporter) for reporter in simulation.reporters):
            append = int(progress.get("step", 0)) > 0
            simulation.reporters.append(
                DCDReporter(str(folder / "production.dcd"), REPORT_STEPS, append=append)
            )
            simulation.reporters.append(
                StateDataReporter(
                    str(folder / "production.log"),
                    REPORT_STEPS,
                    step=True,
                    temperature=True,
                    potentialEnergy=True,
                    speed=True,
                    remainingTime=True,
                    totalSteps=steps_for_ps(plan["production_ns"] * 1000.0),
                    append=append,
                )
            )
        run_chunks(
            simulation,
            steps_for_ps(plan["production_ns"] * 1000.0),
            checkpoint,
            progress_path,
            progress,
        )
        progress["phase"] = "done"
        progress_path.write_text(json.dumps(progress), encoding="utf-8")


def main():
    if not (RUN / "system.xml").exists():
        raise SystemExit("run prepare_system.py first")
    pdb, system, selections = load_inputs()
    if (RUN / "run_plan.json").exists():
        plan = json.loads((RUN / "run_plan.json").read_text(encoding="utf-8"))
    else:
        probe = XmlSerializer.deserialize(XmlSerializer.serialize(system))
        add_restraint(probe, pdb.topology, pdb.positions, selections)
        add_barostat(probe, SEEDS[0])
        plan = benchmark(pdb.topology, probe, pdb.positions)
    for seed in SEEDS:
        done = RUN / f"rep_{seed}" / "progress.json"
        if done.exists() and json.loads(done.read_text(encoding="utf-8")).get("phase") == "done":
            print("already done", seed, flush=True)
            continue
        run_replicate(pdb.topology, system, pdb.positions, selections, seed, plan)


if __name__ == "__main__":
    main()

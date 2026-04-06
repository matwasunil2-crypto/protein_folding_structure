"""
3D Protein Folding Simulation (FINAL VERSION)
- HP Model (3D lattice)
- Monte Carlo + Simulated Annealing
- Energy + Temperature plots
- Auto-save output (Windows compatible)
"""

import numpy as np
import matplotlib.pyplot as plt
import random
import math
import copy
import os

# ─────────────────────────────────────────────
# ENERGY FUNCTION
# ─────────────────────────────────────────────
def compute_energy(positions, sequence):
    energy = 0
    pos_map = {p: i for i, p in enumerate(positions)}

    for i, (x, y, z) in enumerate(positions):
        if sequence[i] != 'H':
            continue
        for dx, dy, dz in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]:
            nb = (x+dx, y+dy, z+dz)
            j = pos_map.get(nb, -1)
            if j != -1 and abs(i - j) > 1 and sequence[j] == 'H':
                energy -= 1

    return energy // 2


# ─────────────────────────────────────────────
# RANDOM STRUCTURE
# ─────────────────────────────────────────────
def generate_random_conformation(n):
    directions = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]

    for _ in range(5000):
        pos = [(0,0,0)]
        visited = {(0,0,0)}

        for _ in range(n-1):
            x, y, z = pos[-1]
            random.shuffle(directions)

            placed = False
            for dx, dy, dz in directions:
                nxt = (x+dx, y+dy, z+dz)
                if nxt not in visited:
                    pos.append(nxt)
                    visited.add(nxt)
                    placed = True
                    break

            if not placed:
                break

        if len(pos) == n:
            return pos

    return [(i, 0, 0) for i in range(n)]


# ─────────────────────────────────────────────
# MOVES
# ─────────────────────────────────────────────
def end_move(positions):
    new_pos = positions.copy()
    dirs = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
    visited = set(new_pos)

    idx = random.choice([0, len(positions)-1])
    anchor = new_pos[1] if idx == 0 else new_pos[-2]

    candidates = [(anchor[0]+dx, anchor[1]+dy, anchor[2]+dz)
                  for dx,dy,dz in dirs if (anchor[0]+dx, anchor[1]+dy, anchor[2]+dz) not in visited]

    if not candidates:
        return None

    new_pos[idx] = random.choice(candidates)
    return new_pos


def apply_random_move(positions):
    return end_move(positions)


# ─────────────────────────────────────────────
# SIMULATED ANNEALING
# ─────────────────────────────────────────────
def simulated_annealing(sequence, T_start=5.0, T_end=0.01, steps=20000):
    n = len(sequence)
    pos = generate_random_conformation(n)
    energy = compute_energy(pos, sequence)

    best_pos = pos.copy()
    best_energy = energy

    e_hist = [energy]
    t_hist = [T_start]

    for step in range(steps):
        T = T_start * (T_end / T_start) ** (step / steps)

        new_pos = apply_random_move(pos)
        if new_pos is None:
            continue

        new_energy = compute_energy(new_pos, sequence)
        delta = new_energy - energy

        if delta < 0 or random.random() < math.exp(-delta / max(T, 1e-10)):
            pos = new_pos
            energy = new_energy

        if energy < best_energy:
            best_energy = energy
            best_pos = pos.copy()

        e_hist.append(energy)
        t_hist.append(T)

    return best_pos, best_energy, e_hist, t_hist


# ─────────────────────────────────────────────
# PLOTTING
# ─────────────────────────────────────────────
def plot_results(sequence, positions, e_hist, t_hist):
    fig = plt.figure(figsize=(14, 8))
    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(224)

    xs = [p[0] for p in positions]
    ys = [p[1] for p in positions]
    zs = [p[2] for p in positions]

    # 3D structure
    ax1.plot(xs, ys, zs, '-o')
    for i, (x, y, z) in enumerate(positions):
        color = 'red' if sequence[i] == 'H' else 'blue'
        ax1.scatter(x, y, z, c=color)

    ax1.set_title("3D Protein Fold")

    # Energy
    ax2.plot(e_hist)
    ax2.set_title("Energy")

    # Temperature
    ax3.plot(t_hist)
    ax3.set_title("Temperature")

    # SAVE FIX (IMPORTANT)
    output_dir = "outputs"
    os.makedirs(output_dir, exist_ok=True)

    file_path = os.path.join(output_dir, "protein_folding_result.png")
    plt.savefig(file_path, dpi=150)

    print(f"✅ Saved: {file_path}")

    plt.show()


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
if __name__ == "__main__":
    sequence = "HPHPPHHPHPPHPHHPPHPH"

    print("🧬 Running Simulation...")
    
    best_pos, best_e, e_hist, t_hist = simulated_annealing(sequence)

    print(f"✅ Best Energy: {best_e}")

    plot_results(sequence, best_pos, e_hist, t_hist)
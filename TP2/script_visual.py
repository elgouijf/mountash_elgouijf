#!/usr/bin/env python3
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from collections import defaultdict

def main():
    corps_data = {}
    temps = []

    with open("../systeme_solaire.txt") as f:
        for _ in range(3):
            f.readline()  # sauter les 3 lignes elapsed time

        lignes = f.readlines()

    t = None
    nom = None
    x = None

    for ligne in lignes:
        ligne = ligne.rstrip()
        if not ligne:
            continue
        _, _, value = ligne.split()
        if ligne.startswith("t"):
            t = float(value)
            temps.append(t)
        elif ligne.startswith("nom"):
            nom = value
            if nom not in corps_data:
                corps_data[nom] = {"x" : [],"y" :[]}
        elif ligne.startswith("x"):
            x = float(value)
        elif ligne.startswith("y"):
            y = float(value)
            corps_data[nom]["x"].append(x)
            corps_data[nom]["y"].append(y)


    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_facecolor("black")
    fig.patch.set_facecolor("black")

    couleurs = {
        "soleil":  "yellow",
        "terre":   "deepskyblue",
        "jupiter": "orange",
        "haley":   "white"
    }
    tailles = {
        "soleil":  120,
        "terre":   40,
        "jupiter": 80,
        "haley":   15
    }

    # lignes de trajectoire
    trajectoires = {}
    points = {}
    for nom, couleur in couleurs.items():
        traj, = ax.plot([], [], color=couleur, linewidth=0.5, alpha=0.5)
        point, = ax.plot([], [], 'o', color=couleur,
                         markersize=tailles[nom]**0.5, label=nom)
        trajectoires[nom] = traj
        points[nom] = point

    titre = ax.set_title("t = 0.000", color="white")
    ax.legend(loc="upper right", facecolor="black", labelcolor="white")
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_edgecolor("white")

    # calcul des limites
    all_x = [x for d in corps_data.values() for x in d["x"]]
    all_y = [y for d in corps_data.values() for y in d["y"]]
    margin = 2
    ax.set_xlim(min(all_x) - margin, max(all_x) + margin)
    ax.set_ylim(min(all_y) - margin, max(all_y) + margin)

    n_frames = len(temps)

    def update(frame):
        for nom in couleurs:
            xs = corps_data[nom]["x"][:frame+1]
            ys = corps_data[nom]["y"][:frame+1]
            trajectoires[nom].set_data(xs, ys)
            if xs:
                points[nom].set_data([xs[-1]], [ys[-1]])
        titre.set_text(f"t = {temps[frame]:.3f}")
        return list(trajectoires.values()) + list(points.values()) + [titre]

    ani = animation.FuncAnimation(
        fig, update, frames=n_frames,
        interval=20, blit=True, repeat = False
    )

    plt.tight_layout()
    plt.show()

main()
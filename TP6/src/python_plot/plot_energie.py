#!/usr/bin/env python3

import sys
from pathlib import Path

import matplotlib.pyplot as plt
#### Vous aurez besoin de pandas pour lire le csv, si vous ne l'avez pas installé : pip install pandas
import pandas as pd


def main():
    # 🔹 Chemin du script
    script_dir = Path(__file__).resolve().parent

    # 🔹 Racine du projet (2 niveaux au-dessus : src/python_plot → projet)
    project_root = script_dir.parent.parent

    # 🔹 Chemin par défaut
    default_path = project_root / "energy" / "energie_tp6.csv"

    # 🔹 Argument utilisateur
    if len(sys.argv) >= 2:
        csv_path = Path(sys.argv[1])
    else:
        csv_path = default_path
        print(f"[INFO] Aucun fichier fourni.")
        print(f"[INFO] Utilisation du fichier par défaut : {csv_path}")

    # 🔹 Vérification
    if not csv_path.exists():
        print(f"[ERREUR] Fichier introuvable : {csv_path}")
        return 1

    data = pd.read_csv(csv_path)

    # 🔹 Dossier plots (à la racine du projet)
    output_dir = project_root / "plots"
    output_dir.mkdir(exist_ok=True)

    # 🔹 Plot
    plt.figure()
    plt.plot(data["time"], data["Ec"], label="Ec")
    plt.plot(data["time"], data["Ep"], label="Ep")
    plt.plot(data["time"], data["Em"], label="Em")

    plt.xlabel("Temps")
    plt.ylabel("Énergie")
    plt.title(f"Énergies - {csv_path.name}")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    output_path = output_dir / f"{csv_path.stem}.png"
    plt.savefig(output_path, dpi=200)

    plt.show()

    print(f"[INFO] Figure générée : {output_path}")


if __name__ == "__main__":
    raise SystemExit(main())
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent

def load_frames(filename=None):
    """
    Charge les données des frames de simulation depuis un fichier texte.

    Le fichier doit contenir cinq colonnes :
    - identifiant de frame
    - x
    - y
    - z
    - type de particule

    Parameters
    ----------
    filename : str
        Chemin vers le fichier contenant les frames.

    Returns
    -------
    tuple
        data : np.ndarray
            Tableau contenant toutes les données.
        frame_ids : np.ndarray
            Identifiants uniques des frames.

    Raises
    ------
    FileNotFoundError
        Si le fichier n'existe pas.
    """
    if filename is None:
        filename = PROJECT_ROOT / "frames_3d" / "frames.txt"
    else:
        filename = Path(filename)
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Fichier introuvable : {filename}")

    data = np.loadtxt(filename)
    frame_ids = np.unique(data[:, 0]).astype(int)
    return data, frame_ids


def create_animation(data, frame_ids):
    """
    Crée l'animation matplotlib de la collision.

    Parameters
    ----------
    data : np.ndarray
        Données complètes de la simulation.
    frame_ids : np.ndarray
        Identifiants uniques des frames.
    """
    fig, ax = plt.subplots(figsize=(10, 4))

    scat_carre = ax.scatter([], [], s=4, c="red", label="Carré")
    scat_rect = ax.scatter([], [], s=4, c="blue", label="Rectangle")

    ax.set_xlim(0, 250)
    ax.set_ylim(0, 120)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Collision de deux ensembles de particules")
    ax.legend()
    ax.grid()

    def init():
        """
        Initialise les nuages de points vides.
        """
        scat_carre.set_offsets(np.empty((0, 2)))
        scat_rect.set_offsets(np.empty((0, 2)))
        return scat_carre, scat_rect

    def update(frame_index):
        """
        Met à jour l'affichage pour une frame donnée.

        Parameters
        ----------
        frame_index : int
            Indice de la frame dans le tableau frame_ids.
        """
        fid = frame_ids[frame_index]
        frame_data = data[data[:, 0] == fid]

        carre = frame_data[frame_data[:, 4] == 1]
        rect = frame_data[frame_data[:, 4] == 2]

        scat_carre.set_offsets(carre[:, 1:3] if len(carre) else np.empty((0, 2)))
        scat_rect.set_offsets(rect[:, 1:3] if len(rect) else np.empty((0, 2)))

        ax.set_title(f"Collision de deux ensembles de particules - Frame {fid}")
        return scat_carre, scat_rect

    ani = FuncAnimation(
        fig,
        update,
        frames=len(frame_ids),
        init_func=init,
        interval=40,
        blit=True
    )

    plt.show()
    return ani


def main():
    """
    Fonction principale :
    - charge les données de simulation,
    - crée l'animation de la collision.
    """
    try:
        data, frame_ids = load_frames()
        create_animation(data, frame_ids)
    except Exception as e:
        print(f"Erreur : {e}")


if __name__ == "__main__":
    main()
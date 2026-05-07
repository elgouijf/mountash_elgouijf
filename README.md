# Simulation de Particules en C++

Simulateur de dynamique moléculaire construit progressivement au fil de six travaux pratiques. Le code évolue d'un prototype simple vers un simulateur capable de gérer des systèmes de plusieurs dizaines de milliers de particules avec interactions locales, visualisation scientifique, conditions aux limites, gravité et parallélisation OpenMP.

Le modèle physique repose sur le **potentiel de Lennard-Jones** et le schéma d'intégration **Störmer-Verlet**.

---

## Prérequis

- Compilateur C++17 (GCC ≥ 9 ou Clang ≥ 10)
- CMake ≥ 3.20
- OpenMP
- Python 3 + Matplotlib
- ParaView
- Doxygen
- Google Test *(téléchargé automatiquement par CMake via FetchContent)*

---

## Compilation

```bash
cd TP6/
mkdir build && cd build
cmake ..
make -j
```

Les flags de compilation activés sont `-O3`, `-march=native`, `-funroll-loops`, `-ffast-math` et `-flto`.

---

## Exécution

Les exécutables sont générés dans `build/demo/` :

| Exécutable | Description |
|---|---|
| `./demo/simu_tp2` | Simulation orbitale (Soleil, Terre, Jupiter, comète de Halley) |
| `./demo/simu_tp4` | Collision 2D entre deux blocs de particules avec grille de cellules |
| `./demo/simu_tp6` | Collision boule/pavé avec gravité, conditions aux limites et limiteur d'énergie |
| `./demo/test_3d` | Collision 3D entre deux blocs de particules |
| `./demo/test_plot_ur` | Tracé du potentiel de Lennard-Jones $U(r)$ |

Depuis les modifications réalisées au TP6, la majorité des simulations utilisent OpenMP automatiquement. Les seules exceptions sont `simu_tp2` et `test_plot_ur`, qui restent séquentielles.

### Contrôle des threads

```bash
# Fixer le nombre de threads
OMP_NUM_THREADS=6 ./demo/simu_tp6

# Forcer l'exécution séquentielle
OMP_NUM_THREADS=1 ./demo/simu_tp6
```

Les mesures réalisées montrent qu'un compromis proche de **6 threads** donne les meilleures performances sur les machines utilisées pendant le projet.

---

## Visualisation

Les sorties sont générées dans des répertoires dédiés à la racine du projet. Le chemin est affiché dans le terminal à chaque exécution.

| Répertoire | Contenu |
|---|---|
| `frames/` | Sorties texte |
| `vtk_frames/` | Sorties VTK Legacy |
| `vtu_frames/` | Sorties VTU/XML |
| `energy/` | Énergie cinétique avec un .csv au cours du temps |

### Scripts Python (`src/python_plot/`)

```bash
python3 src/python_plot/plot_collision.py      # animation 2D
python3 src/python_plot/plot_collision_3d.py   # animation 3D
python3 src/python_plot/plot_energie.py        # courbes d'énergie
python3 src/python_plot/plot_ur.py             # potentiel de Lennard-Jones
python3 src/python_plot/script_visual.py       # trajectoires orbitales (TP2)
```

Pour les simulations plus lourdes, les sorties **VTK/VTU** sont directement exploitables avec **ParaView** : `File > Open` sur le fichier `.vtu.series` généré automatiquement.

---

## Tests

```bash
cd build/
ctest --output-on-failure
```

Les tests couvrent les classes `vecteur`, `particule`, `cellule` et `univers`, avec une stratégie ascendante : validation des composants de bas niveau avant les tests d'intégration.

---

## Documentation

```bash
cd build/
make doc
```

La documentation HTML est générée dans `doc/html/index.html`.

---

## Organisation du dépôt

```
TP6/
├── include/          # Déclarations des classes principales
├── src/              # Implémentations et programme principal
│   └── python_plot/  # Scripts de visualisation
├── demo/             # Simulations de démonstration
├── test/             # Tests unitaires (Google Test)
├── doc/              # Documentation Doxygen et diagrammes ACVL
└── CMakeLists.txt
```

---

## Auteurs

Faical El Gouij, Mountassir Hamza — ENSIMAG, C++ pour les Mathématiques Appliquées

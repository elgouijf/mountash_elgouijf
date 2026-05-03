#include <iostream>
#include <vector>
#include <cmath>
// Pour enregistrer les frames de la simulation
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono> // pour la perf
#include <cstdlib>
#include <string>

#include "univers.hxx"
#include "io.hxx"
#include "output_paths.hxx"

/**
 * @brief Sauvegarde une frame dans un fichier texte.
 *
 * Chaque ligne du fichier contient :
 * frame_id x y z type
 *
 * @param file Flux de sortie texte.
 * @param u Univers contenant les particules.
 * @param frame_id Identifiant de la frame courante.
 */
void sauvegarde_frame(std::ofstream& file, const univers& u, int frame_id) {
    for (particule* p : u.getParticules()) {
        const vecteur& pos = p->getPosition();
        file << frame_id << " "
             << pos[0] << " "
             << pos[1] << " "
             << pos[2] << " "
             << p->getType() << "\n";
    }
} 



/**
 * @brief Recherche le script Python de visualisation 3D.
 *
 * Plusieurs chemins sont testés afin de permettre l'exécution
 * depuis différents répertoires.
 *
 * @return Le chemin du script si trouvé, une chaîne vide sinon.
 */
std::string trouver_script_python_3d() {
    std::vector<std::string> candidats = {
        std::string(PROJECT_SOURCE_DIR) + "/src/python_plot/plot_collision_3d.py",
        "src/python_plot/plot_collision_3d.py",
        "../src/python_plot/plot_collision_3d.py"
    };

    for (const auto& chemin : candidats) {
        if (std::filesystem::exists(chemin)) {
            return chemin;
        }
    }

    return "";
}

/**
 * @brief Programme principal de simulation 3D.
 *
 * Ce programme :
 * - initialise un cube de particules au-dessus d'un pavé,
 * - simule leur collision avec un potentiel de Lennard-Jones,
 * - sauvegarde les frames soit en texte, soit en VTK,
 * - génère en mode VTK un fichier animation.vtk.series,
 * - lance en mode texte un script Python de visualisation 3D.
 *
 * Modes disponibles :
 * - 't' : export texte + visualisation Python
 * - 'v' : export VTK + génération du fichier animation.vtk.series
 * - 'x' : export VTU + génération du fichier animation.vtu.series
 *
 * @return EXIT_SUCCESS si l'exécution se termine correctement,
 *         EXIT_FAILURE en cas d'erreur.
 */
int main(){

    std::string mode;
    std::cout << "Choisir le mode : (t = txt, v = vtk, x = vtu) : ";
    std::cin >> mode;

    if (mode != 't' && mode != 'v' && mode != 'x') {
        std::cerr << "Mode invalide. Choisir 't', 'v' ou 'x'.\n";
        return EXIT_FAILURE;
    }

    // Caractéristiques de l'univers
    double sigma = 1.0;
    double epsilon = 5.0;
    double mass = 1.0;
    //double dt = 0.00005;
    double dt = 0.001;
    double rcut = 2.5 * sigma;
    int dim = 3;
    std::vector<double> Lds = {120.0, 120.0, 120.0}; // L1, L2 et L3
    double dist_entre_particules = pow(2, 1/6.0)/sigma; // distance entre les particules, donc chaque particule occupe
                                                        // une sphère de rayon dist_entre_particules/2, pour que les particules
    std::vector<particule*> particules;

    // Simulation duration and number of frames
    //double duration = 19.5; // durée de la simulation en secondes
    double duration = 10.0;
    int num_frames = duration / dt;

    // Vitesses des coprs
    // !!!!! Attention l'objet vecteur et conçu pour être un vecteur 3D, donc ici cette fois on utilise vraiment les 3 dimensions
    // Cei étant dit notre univers est équipé de projection dimensionnelle, donc il s'occupera de ne garder que les bonnes dimensions selon le cas
    vecteur v2;
    vecteur v1(0.0, 0.0, -8.0); // chute => vecteur vitesse négatif dans la direction z

    // Nombre de particules
    int N1x = 10; // nombre de particules selon x pour le cube du haut
    int N1y = 10; // nombre de particules selon y pour le cube du haut
    int N1z = 10; // nombre de particules selon z pour le cube du haut

    int N2x = 10; // nombre de particules selon x pour le pavé du bas
    int N2y = 30; // nombre de particules selon y pour le pavé du bas
    int N2z = 10;  // nombre de particules selon z pour le pavé du bas

    int N1 = N1x * N1y * N1z; // nombre de particules du cube
    int N2 = N2x * N2y * N2z; // nombre de particules du pavé

    // Création de l'univers vide (de la matière)
    univers uni(particules, Lds, rcut, dim, epsilon, sigma);
    // Création des particules
    int id = 0;

    double largeur_cube_x = (N1x - 1) * dist_entre_particules;
    double largeur_cube_y = (N1y - 1) * dist_entre_particules;

    double largeur_pave_x = (N2x - 1) * dist_entre_particules;
    double largeur_pave_y = (N2y - 1) * dist_entre_particules;

    double x0_pave = 20.0;
    double y0_pave = 20.0;
    double z0_pave = 15.0;

    double x0_cube = x0_pave + (largeur_pave_x - largeur_cube_x) / 2.0;
    double y0_cube = y0_pave + (largeur_pave_y - largeur_cube_y) / 2.0;
    double z0_cube = z0_pave + N2z * dist_entre_particules + 5.0;

    //!!!!!!!!!!!! Important: Sauvegarder chauqe les x frames !!!!!!!!!!!!!!!!!!!!!!
    int save_every = 200;

    for (int k = 0; k < N1z; ++k) {
        for (int j = 0; j < N1y; ++j) {
            for (int i = 0; i < N1x; ++i) {
                // simplification : on place les particules sur une grille régulière au lieu d'une sphère
                double x = x0_cube + i * dist_entre_particules + dist_entre_particules / 2.0; // centre du cube juste au dessus de celui du pavé
                double y = y0_cube + j * dist_entre_particules + dist_entre_particules / 2.0;
                double z = z0_cube + k * dist_entre_particules + dist_entre_particules / 2.0;
                vecteur pos = vecteur(x, y, z);
                uni.ajoute_particule(new particule(id++, 1, mass, pos, v1));
            }
        }
    }

    for (int k = 0; k < N2z; ++k) {
        for (int j = 0; j < N2y; ++j) {
            for (int i = 0; i < N2x; ++i) {
                double x = x0_pave + i * dist_entre_particules + dist_entre_particules / 2.0;
                double y = y0_pave + j * dist_entre_particules + dist_entre_particules / 2.0;
                double z = z0_pave + k * dist_entre_particules + dist_entre_particules / 2.0;
                vecteur pos = vecteur(x, y, z);
                uni.ajoute_particule(new particule(id++, 2, mass, pos, v2));
            }
        }
    }

    std::ofstream file;

    std::filesystem::path dossier_frames;
    std::filesystem::path dossier_vtk;
    std::filesystem::path dossier_vtu;

    if (mode == 't') {
        dossier_frames = ensure_output_dir(OutputType::FramesTxt3D);

        file.open(dossier_frames / "frames.txt");

        if (!file.is_open()) {
            std::cerr << "Impossible d'ouvrir "
                    << dossier_frames / "frames.txt" << "\n";
            return EXIT_FAILURE;
        }
    }

    if (mode == 'v') {
        dossier_vtk = ensure_output_dir(OutputType::FramesVTK3D);
    }

    if (mode == 'x') {
        dossier_vtu = ensure_output_dir(OutputType::FramesVTU3D);
    }

    int frame_id = 0;
    auto start = std::chrono::high_resolution_clock::now();

    uni.calcule_forces();

    // Boucle de simulation
    for (int frame = 0; frame < num_frames; ++frame) {
        //uni.calcule_forces();
        uni.evolue_particules(dt);

        if (frame % 1000 == 0) {
            std::cout << "Frame " << frame << "/" << num_frames << "\n";
        }

        if (frame % save_every == 0) {
            if (mode == 't') {
                sauvegarde_frame(file, uni, frame_id);
            }
            else if (mode == 'v') {
                sauvegarde_frame_vtk(uni, frame_id, dossier_vtk.string());
            }
            else if (mode == 'x') {
                sauvegarde_frame_vtu(uni, frame_id, dossier_vtu.string());
            }
            frame_id++;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Temps de simulation : " << elapsed.count() << " secondes\n";

    if (mode == 't') {
        file.close();

        std::string script_python = trouver_script_python_3d();
        if (script_python.empty()) {
            std::cerr << "Impossible de trouver plot_collision_3d.py\n";
            return EXIT_FAILURE;
        }

        std::cout << "Simulation 3D terminee.\n";
        std::cout << "Pour visualiser les frames 3D texte, lancez :\n"
                  << "python3 src/python_plot/plot_collision_3d.py\n";
    }

    if (mode == 'v') {
    ecrire_fichier_series_json(frame_id, dt, save_every, dossier_vtk.string(), "vtk");
    std::cout << "Fichier de series genere : "
              << (dossier_vtk / "animation.vtk.series") << "\n";
    }

    if (mode == 'x') {
        ecrire_fichier_series_json(frame_id, dt, save_every, dossier_vtu.string(), "vtu");
        std::cout << "Fichier de series genere : "
                << (dossier_vtu / "animation.vtu.series") << "\n";
    }

    return EXIT_SUCCESS;
}
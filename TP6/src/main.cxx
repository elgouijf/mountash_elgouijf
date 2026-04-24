/**
 * @file test_simulation_collision_2d.cxx
 * @brief Simulation 2D d'une collision entre un carré et un rectangle de particules.
 *
 * Ce programme :
 * - initialise deux ensembles de particules en 2D,
 * - simule leur interaction via un potentiel de Lennard-Jones,
 * - sauvegarde les résultats soit en format texte, soit en format VTK,
 * - permet une visualisation avec Python ou ParaView.
 *
 * Modes disponibles :
 * - 't' : export texte + visualisation Python
 * - 'v' : export VTK + visualisation ParaView
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <cstdlib>


#include "univers.hxx"
#include "io.hxx"



/**
 * @brief Recherche le script Python de visualisation 2D.
 *
 * Plusieurs chemins possibles sont testés pour permettre
 * une exécution depuis différents répertoires.
 *
 * @return Chemin vers le script si trouvé, chaîne vide sinon.
 */
std::string trouver_script_python() {
    std::vector<std::string> candidats = {
        std::string(PROJECT_SOURCE_DIR) + "/src/python_plot/plot_collision.py",
        "src/python_plot/plot_collision.py",
        "../src/python_plot/plot_collision.py"
    };

    for (const auto& chemin : candidats) {
        if (std::filesystem::exists(chemin)) {
            return chemin;
        }
    }

    return "";
}

/**
 * @brief Programme principal de simulation 2D.
 *
 * Étapes :
 * - choix du mode de sortie,
 * - initialisation des particules,
 * - simulation temporelle avec Störmer-Verlet,
 * - sauvegarde des frames,
 * - visualisation avec Python ou ParaView.
 *
 * @return EXIT_SUCCESS si la simulation s'exécute correctement,
 *         EXIT_FAILURE sinon.
 */
int main(){

    char mode;
    std::cout << "Choisir le mode : (t = txt, v = vtk legacy, x = vtu xml) : ";
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
    int dim = 2;
    std::vector<double> Lds = {250.0, 150.0};
    double dist_entre_particules = pow(2, 1/6.0)/sigma;
    std::vector<particule*> particules;

    // Durée de la simulation
    //double duration = 19.5;
    double duration = 13.0;
    int num_frames = duration / dt;

    // Vitesses initiales
    vecteur v2;
    vecteur v1(0.0, -10, 0.0);

    // Nombre de particules
    int N1 = 40*40;
    int N2 = 160*40;

    // Création de l'univers
    univers uni(particules, Lds, rcut, dim, epsilon, sigma);

    int id = 0;
    double largeur_carre = 39 * dist_entre_particules;
    double largeur_rect = 159 * dist_entre_particules;

    // Sauvegarde d'une frame toutes les save_every itérations
    int save_every = 100;

    // Création du carré supérieur
    for (int i = 0; i < N1; ++i) {
        double x = 20.0 + (largeur_rect - largeur_carre) / 2.0 + (i % 40) * dist_entre_particules + dist_entre_particules / 2.0;
        double y = 73.0 + (i / 40) * dist_entre_particules + dist_entre_particules / 2.0;
        vecteur pos = vecteur(x, y, 0.0);
        uni.ajoute_particule(new particule(id++, 1, mass, pos, v1));
    }

    // Création du rectangle inférieur
    for (int i = 0; i < N2; ++i) {
        double x = 20.0 + (i % 160) * dist_entre_particules + dist_entre_particules / 2.0;
        double y = 20.5 + (i / 160) * dist_entre_particules + dist_entre_particules / 2.0;
        vecteur pos = vecteur(x, y, 0.0);
        uni.ajoute_particule(new particule(id++, 2, mass, pos, v2));
    }

    std::ofstream file;
    if (mode == 't') {
        std::filesystem::create_directories("frames");
        file.open("frames/frames.txt");
        if (!file.is_open()) {
            std::cerr << "Impossible d'ouvrir frames/frames.txt\n";
            return EXIT_FAILURE;
        }
    }
    std::string dossier_vtk = std::string(PROJECT_SOURCE_DIR) + "/vtk_frames";
    std::string dossier_vtu = std::string(PROJECT_SOURCE_DIR) + "/vtu_frames";

    if (mode == 'v') {
        std::filesystem::create_directories(dossier_vtk);
    }

    if (mode == 'x') {
        std::filesystem::create_directories(dossier_vtu);
    }

    int frame_id = 0;
    auto start = std::chrono::high_resolution_clock::now();

    // Initialisation des forces au temps initial
    uni.calcule_forces();

    // Boucle de simulation
    for (int frame = 0; frame < num_frames; ++frame) {
        uni.evolue_particules(dt);

        if (frame % 1000 == 0) {
            std::cout << "Frame " << frame << "/" << num_frames << "\n";
        }

    if (frame % save_every == 0) {
        if (mode == 't') {
            sauvegarde_frame_txt(file, uni, frame_id);
        } else if (mode == 'v') {
            sauvegarde_frame_vtk(uni, frame_id, dossier_vtk);
        } else if (mode == 'x') {
            sauvegarde_frame_vtu(uni, frame_id, dossier_vtu);
        }
        frame_id++;
    }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Temps de simulation : " << elapsed.count() << " secondes\n";

    if (mode == 't') {
        file.close();

        std::string script_python = trouver_script_python();
        if (script_python.empty()) {
            std::cerr << "Impossible de trouver src/plot_collision.py\n";
            return EXIT_FAILURE;
        }

        std::string commande_python = "python3 " + script_python;
        int code_python = system(commande_python.c_str());
        if (code_python != 0) {
            std::cerr << "Erreur lors de l'execution du script Python.\n";
            return EXIT_FAILURE;
        }
    }

    if (mode == 'v') {
        ecrire_fichier_series_json(frame_id, dt, save_every, dossier_vtk, "vtk");
        std::cout << "Fichier de series genere : "
          << std::filesystem::absolute(dossier_vtk + "/animation.vtk.series")
          << "\n";
    }

    if (mode == 'x') {
        ecrire_fichier_series_json(frame_id, dt, save_every, dossier_vtu, "vtu");
        std::cout << "Fichier de series genere : "
          << std::filesystem::absolute(dossier_vtu + "/animation.vtk.series")
          << "\n";
    }



    return EXIT_SUCCESS;
}
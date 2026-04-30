/**
 * @file simu_tp6.cxx
 * @brief Simulation 2D d'une collision entre une boule et un pavé de particules.
 *
 * Ce programme :
 * - initialise une boule et un pavé de particules en 2D,
 * - simule leur interaction via un potentiel de Lennard-Jones,
 * - ajoute un champ gravitationnel,
 * - sauvegarde les résultats soit en format texte, VTK legacy ou VTU XML.
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

using namespace std;

std::string trouver_script_python() {
    std::vector<std::string> candidats = {
        std::string(PROJECT_SOURCE_DIR) + "/src/python_plot/plot_collision.py",
        "src/python_plot/plot_collision.py",
        "../src/python_plot/plot_collision.py"
    };
    for (const auto& chemin : candidats) {
        if (std::filesystem::exists(chemin)) return chemin;
    }
    return "";
}

ConditionLimite convert_string_to_Cdl(string s) {
    if (s == "a") return ConditionLimite::Absorbante;
    if (s == "p") return ConditionLimite::Periodique;
    return ConditionLimite::Reflexive;
}

int main() {
    std::string mode;
    std::cout << "Choisir le mode : (t = txt, v = vtk legacy, x = vtu xml) : ";
    std::cin >> mode;

    if (mode != "t" && mode != "v" && mode != "x") {
        std::cerr << "Mode invalide. Choisir 't', 'v' ou 'x'.\n";
        return EXIT_FAILURE;
    }

    int dim = 2;

    double sigma   = 1.0;
    double epsilon = 1.0;
    double mass_pave    = 2.0;
    double mass_boule   = 1.0;
    double dt      = 0.0001;
    double rcut    = 2.5 * sigma;
    double G       = -12.0;

    std::vector<double> Lds = {250.0, 180.0};
    std::vector<particule*> particules;

    double dist_entre_particules = std::pow(2.0, 1.0 / 6.0) / sigma;

    // Debug : durée courte. Pour le rendu final, mettre 29.5.
    double duration = 30.0;
    int num_frames = static_cast<int>(duration / dt);

    vecteur v2(0.0, 0.0, 0.0);
    vecteur v1(0.0, -10.0, 0.0);

    cout << "Conditions aux limites (r/a/p) pour xmin xmax ymin ymax : ";
    vector<string> liste_condition;
    for (int i = 0; i < 2 * dim; i++) {
        string condition_l;
        cin >> condition_l;
        liste_condition.push_back(condition_l);
    }

    /* ++++++++++++++++++++ Potentiel mur ++++++++++++++++++ */
    bool utiliser_potentiel_mur = true;
    univers uni(particules, Lds, rcut, dim, epsilon, sigma, G, utiliser_potentiel_mur);

    uni.setConditionsLimites(
        convert_string_to_Cdl(liste_condition[0]),
        convert_string_to_Cdl(liste_condition[1]),
        convert_string_to_Cdl(liste_condition[2]),
        convert_string_to_Cdl(liste_condition[3]),
        ConditionLimite::Reflexive, // peut importe car on est en 2D
        ConditionLimite::Reflexive
    );

    // Pour debugger da'abord la réflexion géométrique :
    //uni.setUtiliserPotentielMur(false);

    int id = 0;

    // Petites tailles pour debug
    int N2x = 20;
    int N2y = 5;

    double x0_pave = 80.0;
    double y0_pave = 5.0;

    // Pavé inférieur
    int N2 = N2x * N2y;
    for (int j = 0; j < N2y; ++j) {
        for (int i = 0; i < N2x; ++i) {
            double x = x0_pave + i * dist_entre_particules;
            double y = y0_pave + j * dist_entre_particules;

            uni.ajoute_particule(
                new particule(id++, 2, mass_pave, vecteur(x, y, 0.0), v2)
            );
        }
    }

    // Boule supérieure
    double cx = x0_pave + 0.5 * (N2x - 1) * dist_entre_particules;
    double cy = 100.0;
    double rayon_boule = 0.0;

    int N1 = 0;
    for (double y = cy - rayon_boule; y <= cy + rayon_boule; y += dist_entre_particules) {
        for (double x = cx - rayon_boule; x <= cx + rayon_boule; x += dist_entre_particules) {
            double dx = x - cx;
            double dy = y - cy;

            if (dx * dx + dy * dy <= rayon_boule * rayon_boule) {
                uni.ajoute_particule(
                    new particule(id++, 1, mass_boule, vecteur(x, y, 0.0), v1)
                );
                N1++;
            }
        }
    }

    std::cout << "Particules boule : " << N1 << "\n";
    std::cout << "Particules pave  : " << N2 << "\n";

    std::ofstream file;
    if (mode == "t") {
        std::filesystem::create_directories("frames");
        file.open("frames/frames.txt");
        if (!file.is_open()) {
            std::cerr << "Impossible d'ouvrir frames/frames.txt\n";
            return EXIT_FAILURE;
        }
    }

    string dossier_vtk;
    string dossier_vtu;

    if (mode == "v") {
        dossier_vtk = "vtk_frames";
        std::filesystem::create_directories(dossier_vtk);
    }

    if (mode == "x") {
        dossier_vtu = "vtu_frames";
        std::filesystem::create_directories(dossier_vtu);
    }

    int save_every = 2000;
    int frame_id = 0;
    auto start = std::chrono::high_resolution_clock::now();


    //uni.calcule_forces();
    // Initialisation des forces
    uni.calcule_forces();

    if (utiliser_potentiel_mur) {
        uni.applique_potentiel_mur();
    }

    uni.applique_gravite();

    bool activer_limiteur = false; // msut be true for larger systems

    for (int frame = 0; frame < num_frames; ++frame) {
        uni.evolue_particules(dt);

        if (frame % 1000 == 0) {
            std::cout << "Frame " << frame << "/" << num_frames
                    << "  particules : " << uni.getNumParticules() << "\n";

            if (activer_limiteur && frame > 0) {
                uni.limite_vitesses(N1, N2);
            }
        }

        if (frame % save_every == 0) {
            if (mode == "t") {
                sauvegarde_frame_txt(file, uni, frame_id);
            } else if (mode == "v") {
                sauvegarde_frame_vtk(uni, frame_id, dossier_vtk);
            } else if (mode == "x") {
                sauvegarde_frame_vtu(uni, frame_id, dossier_vtu);
            }
            frame_id++;
        }
    }


    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Temps de simulation : " << elapsed.count() << " secondes\n";

    if (mode == "t") {
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

    if (mode == "v") {
        ecrire_fichier_series_json(frame_id, dt, save_every, dossier_vtk, "vtk");
        std::cout << "Fichier de series genere : " << dossier_vtk << "/animation.vtk.series\n";
    }

    if (mode == "x") {
        ecrire_fichier_series_json(frame_id, dt, save_every, dossier_vtu, "vtu");
        std::cout << "Fichier de series genere : " << dossier_vtu << "/animation.vtk.series\n";
    }

    return EXIT_SUCCESS;
}
/**
 * @file simu_tp6.cxx
 * @brief TP6 Question 6 : collision entre deux objets avec gravité et limiteur d'énergie.
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
#include <algorithm>

#include "univers.hxx"
#include "io.hxx"

using namespace std;

int main() {
    std::string mode;
    std::cout << "Choisir le mode : (t = txt, v = vtk legacy, x = vtu xml) : ";
    std::cin >> mode;

    if (mode != "t" && mode != "v" && mode != "x") {
        std::cerr << "Mode invalide. Choisir 't', 'v' ou 'x'.\n";
        return EXIT_FAILURE;
    }

    // ==========================
    // Paramètres TP6 Q6
    // ==========================
    int dim = 2;

    double sigma   = 1.0;
    double epsilon = 1.0;
    double masse   = 1.0;
    double dt      = 0.0005;
    double rcut    = 2.5 * sigma;
    double G       = -12.0;

    std::vector<double> Lds = {250.0, 180.0};
    std::vector<particule*> particules;

    double dist_entre_particules = std::pow(2.0, 1.0 / 6.0) * sigma;

    double duration = 29.5;
    int num_frames = static_cast<int>(duration / dt);

    vecteur v_pave(0.0, 0.0, 0.0);
    vecteur v_boule(0.0, -10.0, 0.0);

    bool utiliser_potentiel_mur = true;

    univers uni(particules, Lds, rcut, dim, epsilon, sigma, G, utiliser_potentiel_mur);

    uni.setConditionsLimites(
        ConditionLimite::Reflexive,
        ConditionLimite::Reflexive,
        ConditionLimite::Reflexive,
        ConditionLimite::Reflexive,
        ConditionLimite::Reflexive,
        ConditionLimite::Reflexive
    );

    int id = 0;

    // ==========================
    // Pavé inférieur : N2 = 17227 = 161 * 107
    // ==========================
    int N2x = 161;
    int N2y = 107;
    int N2 = N2x * N2y;

    double largeur_pave = (N2x - 1) * dist_entre_particules;
    double x0_pave = 0.5 * (Lds[0] - largeur_pave);
    double y0_pave = 1.0;

    for (int j = 0; j < N2y; ++j) {
        for (int i = 0; i < N2x; ++i) {
            double x = x0_pave + i * dist_entre_particules;
            double y = y0_pave + j * dist_entre_particules;

            uni.ajoute_particule(
                new particule(id++, 2, masse, vecteur(x, y, 0.0), v_pave)
            );
        }
    }

    // ==========================
    // Boule supérieure : N1 = 395
    // ==========================
    int N1_cible = 395;
    int N1 = 0;

    double cx = Lds[0] / 2.0;
    double cy = 155.0;
    double rayon_recherche = 20.0;

    std::vector<std::pair<double, vecteur>> candidats;

    for (double y = cy - rayon_recherche; y <= cy + rayon_recherche; y += dist_entre_particules) {
        for (double x = cx - rayon_recherche; x <= cx + rayon_recherche; x += dist_entre_particules) {
            double dx = x - cx;
            double dy = y - cy;
            double d2 = dx * dx + dy * dy;

            candidats.push_back({d2, vecteur(x, y, 0.0)});
        }
    }

    std::sort(candidats.begin(), candidats.end(),
              [](const auto& a, const auto& b) {
                  return a.first < b.first;
              });

    for (int k = 0; k < N1_cible; ++k) {
        uni.ajoute_particule(
            new particule(id++, 1, masse, candidats[k].second, v_boule)
        );
        N1++;
    }

    std::cout << "Particules boule : " << N1 << "\n";
    std::cout << "Particules pave  : " << N2 << "\n";
    std::cout << "Total            : " << uni.getNumParticules() << "\n";

    // ==========================
    // Fichiers de sortie
    // ==========================
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

    // ==========================
    // Initialisation des forces
    // ==========================
    uni.calcule_forces();

    if (utiliser_potentiel_mur) {
        uni.applique_potentiel_mur();
    }

    uni.applique_gravite();

    // ==========================
    // Boucle temporelle
    // ==========================
    int save_every = 100;
    int frame_id = 0;

    bool activer_limiteur = true;

    auto start = std::chrono::high_resolution_clock::now();

    for (int frame = 0; frame < num_frames; ++frame) {
        uni.evolue_particules(dt);

        if (frame % 1000 == 0) {
            std::cout << "Frame " << frame << "/" << num_frames
                      << "  particules : " << uni.getNumParticules()
                      << "  Ec = " << uni.energie_cinetique()
                      << "\n";

            if (activer_limiteur && frame > 0) {
                uni.limite_vitesses(N1, N2);
            }
        }

        if (frame % save_every == 0) {
            if (mode == "t") {
                sauvegarde_frame_txt(file, uni, frame_id);
            }
            else if (mode == "v") {
                sauvegarde_frame_vtk(uni, frame_id, dossier_vtk);
            }
            else if (mode == "x") {
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
    }

    if (mode == "v") {
        ecrire_fichier_series_json(frame_id, dt, save_every, dossier_vtk, "vtk");
        std::cout << "Fichier de series genere : "
                  << dossier_vtk << "/animation.vtk.series\n";
    }

    if (mode == "x") {
        ecrire_fichier_series_json(frame_id, dt, save_every, dossier_vtu, "vtu");
        std::cout << "Fichier de series genere : "
                  << dossier_vtu << "/animation.vtk.series\n";
    }

    return EXIT_SUCCESS;
}
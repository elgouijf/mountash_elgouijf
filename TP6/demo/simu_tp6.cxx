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
#include "output_paths.hxx"

using namespace std;

ConditionLimite lire_condition(char c) {
    switch (c) {
        case 'r': return ConditionLimite::Reflexive;
        case 'a': return ConditionLimite::Absorbante;
        case 'p': return ConditionLimite::Periodique;
        case 'z': return ConditionLimite::Aucune;
        default:
            std::cerr << "Condition inconnue, Reflexive par defaut.\n";
            return ConditionLimite::Reflexive;
    }
}


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
    int frame_debut_limiteur = 6000;
    int periode_limiteur = 1000;          // appliquer tous les 1000 pas après ce début
    int periode_debug_energie = 1000;     
    int dim = 2;

    double sigma   = 1.0;
    double epsilon = 1.0;
    double masse   = 2.0;
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

    char cxmin = 'r';
    char cxmax = 'r';
    
    char cymin = 'r';
    char cymax = 'r';

    std::cout << "\nConditions limites possibles :\n";
    std::cout << "  r = reflexive\n";
    std::cout << "  a = absorbante\n";
    std::cout << "  p = periodique\n";
    std::cout << "  z = aucune\n";

    std::cout << "Condition xmin [r/a/p/z] (defaut r) : ";
    std::cin >> cxmin;

    std::cout << "Condition xmax [r/a/p/z] (defaut r) : ";
    std::cin >> cxmax;

    std::cout << "Condition ymin [r/a/p/z] (defaut r) : ";
    std::cin >> cymin;

    std::cout << "Condition ymax [r/a/p/z] (defaut r) : ";
    std::cin >> cymax;

    uni.setConditionsLimites(
        lire_condition(cxmin),
        lire_condition(cxmax),
        lire_condition(cymin),
        lire_condition(cymax),
        ConditionLimite::Reflexive, // zmin
        ConditionLimite::Reflexive  // zmax
    );

    uni.setUtiliserListeVerlet(true, 0.5);

    std::cout << "\nConditions limites appliquees :\n";
    uni.afficherConditionsLimites();
    int id = 0;

    // ==========================
    // Pavé inférieur : N2 = 17227 = 161 * 107
    // ==========================
    double largeur_pave_voulue = 200.0;

    int N2_cible = 17227;
    int N2x = static_cast<int>(largeur_pave_voulue / dist_entre_particules) + 1;
    int N2y = static_cast<int>(std::round(static_cast<double>(N2_cible) / N2x));
    int N2 = N2x * N2y;

    double largeur_pave = (N2x - 1) * dist_entre_particules;
    double hauteur_pave = (N2y - 1) * dist_entre_particules;

    double x0_pave = 0.5 * (Lds[0] - largeur_pave);
    double y0_pave = 1.0;
    double y_top_pave = y0_pave + hauteur_pave;

    for (int j = 0; j < N2y; ++j) {
        for (int i = 0; i < N2x; ++i) {
            double x = x0_pave + i * dist_entre_particules;
            double y = y0_pave + j * dist_entre_particules;

            uni.ajoute_particule(
                /* les particules du pavé ont une masse légèrement supérieure pour mieux différencier
                  les deux objets dans ParaView en utilisant la coloeation par masse */
                new particule(id++, 2, masse + 0.001, vecteur(x, y, 0.0), v_pave)
            );
        }
    }

    // ==========================
    // Boule supérieure : N1 = 395
    // ==========================
    int N1_cible = 395;
    int N1 = 0;

    double cx = Lds[0] / 2.0;
    double cy = y_top_pave + 25.0;
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
    // Enérgies mécaniques
    std::filesystem::path dossier_energy = ensure_output_dir(OutputType::Energy);

    std::ofstream energy_file(dossier_energy / "energie_tp6.csv"); // fichier d'énergie 

    if (!energy_file.is_open()) {
        std::cerr << "Impossible d'ouvrir "
                << dossier_energy / "energie_tp6.csv" << "\n";
        return EXIT_FAILURE;
    }

    ecrire_entete_energie(energy_file);

    std::ofstream file;

    std::filesystem::path dossier_frames;
    std::filesystem::path dossier_vtk;
    std::filesystem::path dossier_vtu;

    if (mode == "t") {
        dossier_frames = ensure_output_dir(OutputType::FramesTxt);

        file.open(dossier_frames / "frames.txt");

        if (!file.is_open()) {
            std::cerr << "Impossible d'ouvrir "
                    << dossier_frames / "frames.txt" << "\n";
            return EXIT_FAILURE;
        }
    }

    if (mode == "v") {
        dossier_vtk = ensure_output_dir(OutputType::FramesVTK);
    }

    if (mode == "x") {
        dossier_vtu = ensure_output_dir(OutputType::FramesVTU);
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

        // Limiteur indépendant du debug/énergie
        if (activer_limiteur &&
            frame >= frame_debut_limiteur &&
            (frame - frame_debut_limiteur) % periode_limiteur == 0) {
            
            uni.limite_vitesses(N1, N2);
        }

        // Debug + énergie indépendants
        if (frame % periode_debug_energie == 0) {
/*             uni.debug_cellules();

            double Ec = uni.energie_cinetique();
            double Ep = uni.energie_potentielle();
            double Em = Ec + Ep;

            ecrire_energie(energy_file, frame, frame * dt, Ec, Ep, Em);

            std::cout << "Frame " << frame << "/" << num_frames
                    << "  Em = " << Em << "\n"; */
            std::cout << "Frame " << frame << "/" << num_frames << "\n";
        }


        if (frame % save_every == 0) {
            if (mode == "t") {
                sauvegarde_frame_txt(file, uni, frame_id);
            }
            else if (mode == "v") {
                sauvegarde_frame_vtk(uni, frame_id, dossier_vtk.string());
            }
            else if (mode == "x") {
                sauvegarde_frame_vtu(uni, frame_id, dossier_vtu.string());
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
        ecrire_fichier_series_json(frame_id, dt, save_every, dossier_vtk.string(), "vtk");
        std::cout << "Fichier de series genere : "
                  << (dossier_vtk / "animation.vtk.series") << "\n";
    }

    if (mode == "x") {
        ecrire_fichier_series_json(frame_id, dt, save_every, dossier_vtu.string(), "vtu");
        std::cout << "Fichier de series genere : "
                  << (dossier_vtu / "animation.vtu.series") << "\n";
    }

    energy_file.close();

    std::cout << "Fichier energie genere : "
            << (dossier_energy / "energie_tp6.csv") << "\n";
    return EXIT_SUCCESS;
}
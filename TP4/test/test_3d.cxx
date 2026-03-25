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

#include "../include/univers.hxx"
#include "../include/io.hxx"

#ifndef PROJECT_SOURCE_DIR
#define PROJECT_SOURCE_DIR "."
#endif


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


// Fonction VTK
void sauvegarde_frame_vtk(const univers& u, int frame_id, const std::string& dossier = "vtk_frames_3d") {
    std::ostringstream nom;
    nom << dossier << "/frame_"
        << std::setw(6) << std::setfill('0') << frame_id
        << ".vtk";

    std::ofstream file(nom.str());
    if (!file.is_open()) {
        std::cerr << "Impossible d'ouvrir " << nom.str() << "\n";
        return;
    }

    const std::vector<particule*>& particules = u.getParticules();
    int N = particules.size();

    file << "# vtk DataFile Version 3.0\n";
    file << "Frame " << frame_id << "\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    file << "POINTS " << N << " float\n";
    for (particule* p : particules) {
        const vecteur& pos = p->getPosition();
        file << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
    }

    file << "\nVERTICES " << N << " " << 2 * N << "\n";
    for (int i = 0; i < N; ++i) {
        file << "1 " << i << "\n";
    }

    file << "\nPOINT_DATA " << N << "\n";

    file << "SCALARS type int 1\n";
    file << "LOOKUP_TABLE default\n";
    for (particule* p : particules) {
        file << p->getType() << "\n";
    }
}


void ecrire_fichier_series_json(int nb_frames, double dt, int save_every, const std::string& dossier = "vtk_frames_3d") {
    std::ofstream file(dossier + "/animation.vtk.series");
    if (!file.is_open()) {
        std::cerr << "Impossible d'ouvrir " << dossier << "/animation.vtk.series\n";
        return;
    }

    file << "{\n";
    file << "  \"file-series-version\" : \"1.0\",\n";
    file << "  \"files\" : [\n";

    for (int i = 0; i < nb_frames; ++i) {
        std::ostringstream nom;
        nom << "frame_" << std::setw(6) << std::setfill('0') << i << ".vtk";

        double temps = i * save_every * dt;

        file << "    { \"name\" : \"" << nom.str() << "\", \"time\" : " << temps << " }";
        if (i != nb_frames - 1) file << ",";
        file << "\n";
    }

    file << "  ]\n";
    file << "}\n";
}


// Vérifie si ParaView existe
bool paraview_disponible() {
    return std::system("command -v paraview >/dev/null 2>&1") == 0;
}


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


int main(){

    char mode;
    std::cout << "Choisir le mode : (t = txt, v = vtk) : ";
    std::cin >> mode;

    bool lancer_paraview = false;

    if (mode == 'v') {
        if (paraview_disponible()) {
            lancer_paraview = true;
        } else {
            std::cout << "ParaView n'est pas detecte : soit vouus ne l'avez pas, ou on arrive pas à le trouver.\n";
            std::cout << "Voulez-vous quand meme generer les fichiers VTK (et lancer la simu à la main)? (y/n) : ";
            char rep;
            std::cin >> rep;

            if (rep != 'y') {
                return EXIT_FAILURE;
            }
        }
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
    double duration = 12.0;
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

    int N2x = 20; // nombre de particules selon x pour le pavé du bas
    int N2y = 20; // nombre de particules selon y pour le pavé du bas
    int N2z = 20;  // nombre de particules selon z pour le pavé du bas

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
    int save_every = 100;


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
    if (mode == 't') {
        std::filesystem::create_directories("frames_3d");
        file.open("frames_3d/frames.txt");
        if (!file.is_open()) {
            std::cerr << "Impossible d'ouvrir frames_3d/frames.txt\n";
            return EXIT_FAILURE;
        }
    }

    if (mode == 'v') {
        std::filesystem::create_directories("vtk_frames_3d");
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
            } else {
                sauvegarde_frame_vtk(uni, frame_id, "vtk_frames_3d");
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

        std::string commande_python = "python3 " + script_python;
        int code_python = system(commande_python.c_str());
        if (code_python != 0) {
            std::cerr << "Erreur lors de l'execution du script Python 3D.\n";
            return EXIT_FAILURE;
        }
    }

    if (mode == 'v') {
        ecrire_fichier_series_json(frame_id, dt, save_every, "vtk_frames_3d");
    }

    if (mode == 'v' && lancer_paraview) {
        int code_paraview = system("paraview vtk_frames_3d/animation.vtk.series &");
        if (code_paraview != 0) {
            std::cerr << "Erreur lors du lancement de ParaView.\n";
        }
    }

    return EXIT_SUCCESS;
}
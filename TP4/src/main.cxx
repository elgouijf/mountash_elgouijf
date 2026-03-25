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
void sauvegarde_frame_vtk(const univers& u, int frame_id, const std::string& dossier = "vtk_frames") {
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


void ecrire_fichier_series_json(int nb_frames, double dt, int save_every, const std::string& dossier = "vtk_frames") {
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


// Trouver le script Python dans plusieurs cas possibles
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
    int dim = 2;
    std::vector<double> Lds = {250.0, 150.0}; // L1 et L2
    double dist_entre_particules = pow(2, 1/6.0)/sigma; // distance entre les particules, donc chaque particule occupe
                                                        // une sphère(ici un cercle) de rayon dist_entre_particules/2, pour que les particules
    std::vector<particule*> particules;


    // Simulation duration and number of frames
    //double duration = 19.5; // durée de la simulation en secondes
    double duration = 13.0;
    int num_frames = duration / dt;


    // Vitesses des coprs
    // !!!!! Attention l'objet vecteur et conçu pour être un vecteur 3D, donc même pour un univers 2D ou 1D, il faut lui donner une composante z (même si elle ne sera pas utilisée) 
    // Cei étant dit notre univers est équipé de projection dimensionnelle, donc on peut lui donner des composantes dans les 3 dimensions et il s'occupera de les projeter dans les bonnes dimensions selon le cas
    vecteur v2;
    vecteur v1(0.0, -10, 0.0); // chute => vecteur vitesse négatif dans la direction y

    // Nombre de particules
    int N1 = 40*40; // nombre de particule du carré
    int N2 = 160*40; // nombre de particule du rectangle

    // Création de l'univers vide (de la matière)
    univers uni(particules, Lds, rcut, dim, epsilon, sigma);
    // Création des particules
    int id = 0;
    double largeur_carre = 39 * dist_entre_particules; //39 intervalles entre les 40 particules du carré, donc 39*dist_entre_particules
    double largeur_rect = 159 * dist_entre_particules; // idem


    //!!!!!!!!!!!! Important: Sauvegarder chauqe les x frames !!!!!!!!!!!!!!!!!!!!!!
    int save_every = 100;
    


    for (int i = 0; i < N1; ++i) {
        // simplification : on place les particules sur une grille régulière au lieu d'un cercle
        double x = 20.0 + (largeur_rect - largeur_carre) / 2.0 + (i % 40) * dist_entre_particules + dist_entre_particules / 2.0; // centre du carrée juste au dessus du celui du rectangle
        double y = 73.0 + (i / 40) * dist_entre_particules + dist_entre_particules / 2.0;
        vecteur pos = vecteur(x, y, 0.0);
        uni.ajoute_particule(new particule(id++, 1, mass, pos, v1));
    }

    
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

    if (mode == 'v') {
        std::filesystem::create_directories("vtk_frames");
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
                sauvegarde_frame_vtk(uni, frame_id);
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
        ecrire_fichier_series_json(frame_id, dt, save_every, "vtk_frames");
    }

    if (mode == 'v' && lancer_paraview) {
        int code_paraview = system("paraview vtk_frames/animation.vtk.series &");
        if (code_paraview != 0) {
            std::cerr << "Erreur lors du lancement de ParaView.\n";
        }
    }

    return EXIT_SUCCESS;
}
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <filesystem>
#include <vector>


/**
 * @brief Calcule le potentiel de Lennard-Jones pour une distance donnée.
 *
 * @param r Distance entre deux particules.
 * @return Valeur du potentiel U(r).
 */
double U_r(double r) {
    double epsilon = 1.0;
    double sigma = 1.0;
    return 4 * epsilon * (std::pow(sigma / r, 12) - std::pow(sigma / r, 6));
}

/**
 * @brief Recherche le script Python de tracé du potentiel de Lennard-Jones.
 *
 * Teste plusieurs chemins possibles pour permettre l'exécution
 * depuis différents répertoires.
 *
 * @return Chemin du script s'il est trouvé, chaîne vide sinon.
 */
std::string trouver_script_python_lj() {
    std::vector<std::string> candidats = {
        std::string(PROJECT_SOURCE_DIR) + "/src/python_plot/plot_ur.py",
        "src/python_plot/plot_ur.py",
        "../src/python_plot/plot_ur.py"
    };

    for (const auto& chemin : candidats) {
        if (std::filesystem::exists(chemin)) {
            return chemin;
        }
    }

    return "";
}

/**
 * @brief Génère des données du potentiel de Lennard-Jones et lance leur tracé.
 *
 * Ce programme :
 * - génère un fichier lj.txt contenant (r, U(r))
 * - appelle un script Python pour tracer le potentiel
 *
 * @return 0 si succès, 1 sinon.
 */
int main() {
    std::ofstream file("lj.txt");
    if (!file.is_open()) {
        std::cerr << "Impossible d'ouvrir lj.txt\n";
        return 1;
    }

    for (double r = 0.8; r <= 3.0; r += 0.01) {
        file << r << " " << U_r(r) << "\n";
    }

    file.close();

    std::string script_python = trouver_script_python_lj();
    if (script_python.empty()) {
        std::cerr << "Impossible de trouver plot_ur.py\n";
        return 1;
    }

    std::string commande = "python3 " + script_python;
    int code = std::system(commande.c_str());

    if (code != 0) {
        std::cerr << "Erreur lors de l'execution de plot_ur.py\n";
        return 1;
    }

    return 0;
}
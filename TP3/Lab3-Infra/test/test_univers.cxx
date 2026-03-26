#include "../src/univers.hxx"
#include <chrono>

int main() {
    univers u;
    int n = 1 << 5; // (2^5)^3 particules
    double masse = 1.0;
    vecteur vitesse(0.0, 0.0, 0.0); // vitesse uniforme

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                vecteur position(i / double(n-1), j / double(n-1), k / double(n-1));
                u.ajoute_particule(particule(i*n*n + j*n + k, 0, masse, position, vitesse));
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();
    std::cout << "Temps d'insertion : " << elapsed << " s" << std::endl;

    std::cout << "Univers créé avec " << u.num_particules << " particules" << std::endl;
    
}
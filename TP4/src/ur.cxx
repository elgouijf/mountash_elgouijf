#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

double U_r(double r) {
    double epsilon = 1.0;
    double sigma = 1.0;
    return 4 * epsilon * (std::pow(sigma / r, 12) - std::pow(sigma / r, 6));
}

int main() {
    std::ofstream file("lj.txt");

    for (double r = 0.8; r <= 3.0; r += 0.01) {
        file << r << " " << U_r(r) << "\n";
    }

    file.close();

    int code = std::system("python3 src/python_plot/plot_ur.py");
    if (code != 0) {
        std::cerr << "Erreur lors de l'execution de plot_ur.py\n";
        return 1;
    }

    return 0;
}
#pragma once
#include "../src/univers.hxx"
#include <iostream>
#include <chrono>
#include <atomic>
#include <thread>

void interaction_univ(int k) {
    int tps_limite = 16;
    int n = 1 << k;

    univers u;
    vecteur vitesse(0.0, 0.0, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int l = 0; l < n; ++l) {
                vecteur pos(i / double(n - 1), j / double(n - 1), l / double(n - 1));
                u.ajoute_particule(particule(i*n*n + j*n + l, 0, 1.0, pos, vitesse));
            }
        }
    }

    auto start = std::chrono::high_resolution_clock::now();
    u.calcule_forces(1.0);
    auto end = std::chrono::high_resolution_clock::now();

    double elapsed = std::chrono::duration<double>(end - start).count();

    if (elapsed > tps_limite) {
        std::cout << "Temps limite dépassé de 16s" << std::endl;
    }

    std::cout << "pour 2^" << 3*k << " -> " << n*n*n
              << " particules : " << elapsed << " s" << std::endl;
}

int run_simulation(int k){
    std::atomic<bool> finished(false);

    std::thread t([&]() {
        interaction_univ(k);
        finished = true;
    });

    // attendre max 16s
    for (int i = 0; i < 16; ++i) {
        if (finished) break;
        std::this_thread::sleep_for(std::chrono::seconds(1));
    }

    if (!finished) {
        std::cout << "Timeout global dépassé (16s)" << std::endl;
        std::terminate(); // arrêt brutal
    }

    t.join();

    return 0;
}
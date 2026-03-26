#ifndef UNIVERS_HXX
#define UNIVERS_HXX

#include "particule.hxx"
#include <vector>
#include <iostream> // pour std::ostream
/**
 * @mainpage Simulation de particules en C++
 *
 * @section intro Introduction
 *
 * Ce projet implémente une simulation de particules en interaction.
 * Chaque particule possède une position, une vitesse et une masse,
 * et évolue dans le temps sous l'effet des forces appliquées.
 *
 * @section structure Structure du projet
 *
 * Le projet est organisé autour de plusieurs classes principales :
 *
 * - vecteur : permet de manipuler des vecteurs en 3 dimensions
 * - particule : représente une particule avec ses propriétés physiques
 * - univers : contient un ensemble de particules et gère leurs interactions
 *
 * @section physique Modèle physique
 *
 * Les particules interagissent entre elles via des forces (ex : gravitation).
 * La force entre deux particules dépend de leur distance.
 *
 * L'énergie du système est composée de :
 * - énergie cinétique
 * - énergie potentielle
 *
 * @section algo Algorithme utilisé
 *
 * L'évolution du système est calculée à l'aide de la méthode de Störmer-Verlet.
 *
 * À chaque pas de temps :
 * - mise à jour des positions
 * - calcul des forces
 * - mise à jour des vitesses
 *
 * Cet algorithme est utilisé pour sa stabilité numérique.
 *
 * @section utilisation Utilisation
 *
 * Le programme permet :
 * - de créer un ensemble de particules
 * - de simuler leur évolution dans le temps
 * - d'observer leurs positions et trajectoires
 *
 * @section perf Performances
 *
 * Des tests peuvent être réalisés avec différents nombres de particules
 * afin d'étudier le temps de calcul et la complexité.
 *
 * @section auteurs Auteurs
 *
 * Projet réalisé dans le cadre d'un TP de simulation numérique en C++.
 */

/**
 * @brief Représente un ensemble de particules
 * 
 * Cette classe gère une collection de particules et permet
 * de calculer leurs interactions et leur évolution dans le temps.
 */
class univers {
public:
    std::vector<particule> particules; ///< Liste des particules
    int num_particules;                ///< Nombre de particules

    /// Constructeur par défaut
    univers();

    /**
     * @brief Constructeur à partir d'un vecteur de particules
     * 
     * @param v vecteur de pointeurs vers des particules
     */
    univers(std::vector<particule>& v);

    /**
     * @brief Ajoute une particule à l'univers
     * 
     * @param p particule à ajouter
     */
    void ajoute_particule(particule p);

    /**
     * @brief Fait évoluer toutes les particules
     * 
     * Met à jour leur position et leur vitesse.
     * 
     * @param dt pas de temps
     */
    void evolue_particules(double dt);

    /**
     * @brief Calcule les forces entre particules
     * 
     * Applique les interactions (ex : gravitation) entre toutes les particules.
     * 
     * @param G constante d'interaction
     * @param eps petite valeur pour éviter les divisions par zéro
     */
    void calcule_forces(double G, double eps=1e-9);
};

/**
 * @brief Affichage d'un univers
 * 
 * @param os flux de sortie
 * @param u univers à afficher
 * @return flux modifié
 */
std::ostream& operator<<(std::ostream& os, const univers& u);

#endif // UNIVERS_HXX
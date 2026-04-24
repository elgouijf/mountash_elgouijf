#include "cellule.hxx"

cellule::cellule(){
    // vecteur vide
    this->num_particules = 0;
    this->particules = std::vector<particule*>();
    this->voisins = std::vector<cellule*>(); // pas de copies → pointeurs
}

/**
 * @brief Construit une cellule avec des particules et des voisins donnés.
 *
 * Utilisé lors de l'initialisation de la grille spatiale.
 *
 * @param v Particules initiales de la cellule.
 * @param voisins Cellules voisines.
 */
cellule::cellule(std::vector<particule*>& v, std::vector<cellule*>& voisins){
    this->particules = v;
    this->num_particules = v.size();
    this->voisins = voisins;
}

/**
 * @brief Ajoute une particule dans la cellule.
 *
 * Cette opération est appelée lors du placement des particules
 * dans la grille spatiale.
 *
 * @param p Particule à ajouter.
 */
void cellule::ajoute_particule(particule* p) {
    this->particules.push_back(p);
    this->num_particules++;
}

/**
 * @brief Vide la cellule.
 *
 * Supprime toutes les particules qu'elle contient.
 * Utilisé à chaque pas de temps avant de replacer les particules.
 */
void cellule::vide(){
    this->particules.clear();
    this->num_particules = 0;
}

/**
 * @brief Ajoute une cellule voisine.
 *
 * Les voisins sont utilisés pour limiter les interactions
 * aux cellules adjacentes lors du calcul des forces.
 *
 * @param c Cellule voisine.
 */
void cellule::ajoute_voisin(cellule* c){
    this->voisins.push_back(c);
}

const std::vector<particule*>& cellule::getParticules() const {
    return particules;
}

std::vector<particule*>& cellule::getParticules() {
    return particules;
}

const std::vector<cellule*>& cellule::getVoisins() const {
    return this->voisins;
}

std::vector<cellule*>& cellule::getVoisins() {
    return this->voisins;
}
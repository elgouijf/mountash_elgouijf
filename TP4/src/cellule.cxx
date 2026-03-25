#include "../include/cellule.hxx"

cellule::cellule(){
    // vecteur vide
    this->num_particules = 0;
    this->particules = std::vector<particule*>();
    this->voisins = std::vector<cellule*>(); // je veux pas de copies   
}

cellule::cellule(std::vector<particule*>& v, std::vector<cellule*>& voisins){
    this->particules = v;
    this->num_particules = v.size();
    this->voisins = voisins;
}

void cellule::ajoute_particule(particule* p) {
    this->particules.push_back(p);
    this->num_particules++;
}

void cellule::vide(){
    this->particules.clear();
    this->num_particules = 0;
}

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

#include "univers.hxx"
#include <cmath>
#include <chrono>
univers::univers(){
     // vecteur vide
    this->num_particules = 0;
    this->particules = std::vector<particule>();
}

univers::univers(std::vector<particule>& v){
    this->particules = v;
    this->num_particules = v.size();
}


void univers::ajoute_particule(particule p) {
    this->particules.push_back(p);
    this->num_particules++;
}

void univers::evolue_particules(double dt){
    for (particule p : particules) {
        // mettre à jour la particule
        p.evolue(dt);
    }
        
}
void univers::calcule_forces(double G, double eps){
    // Calculer les forces gravitationnelles entre les particules

    for (size_t i = 0; i < particules.size(); ++i){
    vecteur Fi;
    for (size_t j = i+1; j < particules.size(); ++j){ // j = i+1
        /*if (elapsed > tmps_limite){
            return;
        }*/
        vecteur diff = particules[j].getPosition() - particules[i].getPosition();
        double dist = sqrt(diff.norme()) + eps;
        //double dist2 = diff.norme() * diff.norme() + eps;
        vecteur Fij = diff * (G * particules[i].getMasse() * particules[j].getMasse() / dist / dist / dist);
        Fi = Fi + Fij;
        particules[j].ajouterForce(vecteur()-Fij); // force opposée sur j
    }
    particules[i].setForce(Fi);
}
}

std::ostream& operator<<(std::ostream& os, const univers& u) {
    for (auto p : u.particules){
        os << p << std::endl;
    }
    return os;
}

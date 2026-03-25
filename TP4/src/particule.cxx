#include "../include/particule.hxx"

#include <iostream>


particule:: particule(){
    id = 0;
    type = 0;
    masse = 0.0;
    position = vecteur();
    vitesse = vecteur();
    force = vecteur();
    force_old = vecteur();
    //voisins = std::vector<particule>();
}

particule:: particule(int id,int type,double m,vecteur p,vecteur v){
    this->id = id;
    this->type = type;
    this->masse = m;
    this->position = p;
    this->vitesse = v;
    this->force = vecteur();
    this->force_old = vecteur();
    //this->voisins = voisins;
}


int particule::getId() const { return id; }
const vecteur& particule::getPosition() const { return position; }
const vecteur& particule::getVitesse() const { return vitesse; }
const vecteur& particule::getForce() const { return force; }
const vecteur& particule::getForceOld() const { return force_old; }
double particule::getMasse() const { return masse; }      
int particule::getType() const { return type; }  

void particule::setForce(const vecteur& f) { this->force = f; }
void particule::setForceOld(const vecteur& f) { this->force_old = f; }
void particule::setPosition(const vecteur& p) { this->position = p; }
void particule::setVitesse(const vecteur& v) { this->vitesse = v; }


void particule::evolue(double dt){
    if (masse == 0.0) return;
    vitesse = vitesse + force * dt / masse;
    position = position + vitesse * dt;
    force = vecteur();
}

void particule::avance_position_verlet(double dt) {
    if (masse == 0.0) return;

    force_old = force;
    position = position + vitesse * dt + force * (0.5 * dt * dt / masse);
}

void particule::avance_vitesse_verlet(double dt) {
    if (masse == 0.0) return;

    vitesse = vitesse + (force + force_old) * (0.5 * dt / masse);
}

void particule::ajouterForce(double fx, double fy, double fz) {
    force.ajoute(fx, fy, fz);
}

std::ostream& operator<<(std::ostream& os, const particule& p){
    os << "Particule(ID: " << p.getId() << ", Type: " << p.getType()
       << ", Masse: " << p.getMasse() << ", Position: " << p.getPosition() 
       << ", Vitesse: " << p.getVitesse() << ", Force: " << p.getForce() << ")";
    return os;
}
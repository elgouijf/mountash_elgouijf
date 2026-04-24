#include "particule.hxx"

#include <iostream>

/**
 * @brief Construit une particule nulle.
 *
 * Toutes les grandeurs physiques sont initialisées à zéro.
 * Ce constructeur est principalement utile pour créer une particule
 * par défaut avant initialisation complète.
 */
particule::particule(){
    id = 0;
    type = 0;
    masse = 0.0;
    position = vecteur{};
    vitesse = vecteur{};
    force = vecteur{};
    force_old = vecteur{};
}

/**
 * @brief Construit une particule avec ses caractéristiques physiques initiales.
 *
 * @param id Identifiant de la particule.
 * @param type Type de la particule.
 * @param m Masse de la particule.
 * @param p Position initiale.
 * @param v Vitesse initiale.
 */
particule::particule(int id, int type, double m, vecteur p, vecteur v){
    this->id = id;
    this->type = type;
    this->masse = m;
    this->position = p;
    this->vitesse = v;
    this->force = vecteur{};
    this->force_old = vecteur{};
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

/**
 * @brief Fait évoluer la particule sur un pas de temps avec un schéma simple.
 *
 * La vitesse est d'abord mise à jour à partir de la force actuelle,
 * puis la position est mise à jour avec la nouvelle vitesse.
 *
 * @param dt Pas de temps.
 */
void particule::evolue(double dt){
    if (masse == 0.0) return;
    vitesse = vitesse + force * dt / masse;
    position = position + vitesse * dt;
    force = vecteur{};
}

/**
 * @brief Met à jour la position selon la première étape de Störmer-Verlet.
 *
 * L'ancienne force est sauvegardée dans force_old afin d'être réutilisée
 * lors de la mise à jour de la vitesse.
 *
 * @param dt Pas de temps.
 */
void particule::avance_position_verlet(double dt) {
    if (masse == 0.0) return;

    force_old = force;
    position = position + vitesse * dt + force * (0.5 * dt * dt / masse);
}

/**
 * @brief Met à jour la vitesse selon la seconde étape de Störmer-Verlet.
 *
 * La nouvelle vitesse est calculée à partir de la moyenne
 * entre l'ancienne force et la nouvelle force.
 *
 * @param dt Pas de temps.
 */
void particule::avance_vitesse_verlet(double dt) {
    if (masse == 0.0) return;

    vitesse = vitesse + (force + force_old) * (0.5 * dt / masse);
}

/**
 * @brief Ajoute directement des composantes de force à la force courante.
 *
 * Cette version évite la création d'un vecteur temporaire et est donc
 * plus efficace lorsqu'elle est appelée très fréquemment dans le calcul
 * des interactions.
 *
 * @param fx Composante selon x.
 * @param fy Composante selon y.
 * @param fz Composante selon z.
 */
void particule::ajouterForce(double fx, double fy, double fz) {
    force.ajoute(fx, fy, fz);
}

/**
 * @brief Ajoute une force à la particule.
 *
 * Cette méthode ajoute un vecteur force à la force actuelle de la particule.
 * Elle est moins optimisée que la version composante par composante,
 * car elle implique la création d’un vecteur temporaire.
 *
 * @param f Vecteur force à ajouter.
 */
void particule::ajouterForce(const vecteur& f) {
    force = force + f;
}

/**
 * @brief Affiche les principales caractéristiques d'une particule.
 *
 * @param os Flux de sortie.
 * @param p Particule à afficher.
 * @return Le flux de sortie.
 */
std::ostream& operator<<(std::ostream& os, const particule& p){
    os << "Particule(ID: " << p.getId() << ", Type: " << p.getType()
       << ", Masse: " << p.getMasse() << ", Position: " << p.getPosition()
       << ", Vitesse: " << p.getVitesse() << ", Force: " << p.getForce() << ")";
    return os;
}
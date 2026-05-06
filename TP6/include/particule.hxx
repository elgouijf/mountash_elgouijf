#ifndef PARTICULE_HXX
#define PARTICULE_HXX

#include "vecteur.hxx"
#include <iostream>

class particule {
private:
    int id;
    int type;
    double masse;

    vecteur position;
    vecteur vitesse;
    vecteur force;
    vecteur force_old;

    int index_univers = -1;

public:
    /**
     * @brief Construit une particule nulle.
     */
    particule();

    /**
     * @brief Construit une particule avec ses caractéristiques physiques initiales.
     *
     * @param id Identifiant logique de la particule.
     * @param type Type de la particule.
     * @param m Masse de la particule.
     * @param p Position initiale.
     * @param v Vitesse initiale.
     */
    particule(int id, int type, double m, vecteur p, vecteur v);

    /**
     * @brief Retourne l'identifiant logique de la particule.
     *
     * @return Identifiant de la particule.
     */
    inline int getId() const noexcept {
        return id;
    }

    /**
     * @brief Retourne la position de la particule.
     *
     * @return Référence constante vers la position.
     */
    inline const vecteur& getPosition() const noexcept {
        return position;
    }

    /**
     * @brief Retourne la vitesse de la particule.
     *
     * @return Référence constante vers la vitesse.
     */
    inline const vecteur& getVitesse() const noexcept {
        return vitesse;
    }

    /**
     * @brief Retourne la force actuelle appliquée à la particule.
     *
     * @return Référence constante vers la force.
     */
    inline const vecteur& getForce() const noexcept {
        return force;
    }

    /**
     * @brief Retourne la force précédente.
     *
     * @return Référence constante vers l'ancienne force.
     */
    inline const vecteur& getForceOld() const noexcept {
        return force_old;
    }

    /**
     * @brief Retourne la masse de la particule.
     *
     * @return Masse de la particule.
     */
    inline double getMasse() const noexcept {
        return masse;
    }

    /**
     * @brief Retourne le type de la particule.
     *
     * @return Type de la particule.
     */
    inline int getType() const noexcept {
        return type;
    }

    /**
     * @brief Modifie la force actuelle.
     *
     * @param f Nouvelle force.
     */
    inline void setForce(const vecteur& f) noexcept {
        force = f;
    }

    /**
     * @brief Modifie l'ancienne force.
     *
     * @param f Nouvelle ancienne force.
     */
    inline void setForceOld(const vecteur& f) noexcept {
        force_old = f;
    }

    /**
     * @brief Modifie la position.
     *
     * @param p Nouvelle position.
     */
    inline void setPosition(const vecteur& p) noexcept {
        position = p;
    }

    /**
     * @brief Modifie la vitesse.
     *
     * @param v Nouvelle vitesse.
     */
    inline void setVitesse(const vecteur& v) noexcept {
        vitesse = v;
    }

    /**
     * @brief Retourne l'index actuel de la particule dans l'univers.
     *
     * Cet index sert aux buffers OpenMP. Contrairement à l'identifiant logique,
     * il est toujours contigu après reconstruction du tableau des particules.
     *
     * @return Index courant dans le vecteur de particules.
     */
    inline int getIndexUnivers() const noexcept {
        return index_univers;
    }

    /**
     * @brief Définit l'index actuel de la particule dans l'univers.
     *
     * @param i Nouvel index dans le vecteur de particules.
     */
    inline void setIndexUnivers(int i) noexcept {
        index_univers = i;
    }

    /**
     * @brief Fait évoluer la particule avec un schéma simple.
     *
     * @param dt Pas de temps.
     */
    void evolue(double dt);

    /**
     * @brief Met à jour la position avec la première étape de Störmer-Verlet.
     *
     * @param dt Pas de temps.
     */
    void avance_position_verlet(double dt);

    /**
     * @brief Met à jour la vitesse avec la seconde étape de Störmer-Verlet.
     *
     * @param dt Pas de temps.
     */
    void avance_vitesse_verlet(double dt);

    /**
     * @brief Ajoute directement des composantes de force.
     *
     * @param fx Composante selon x.
     * @param fy Composante selon y.
     * @param fz Composante selon z.
     */
    void ajouterForce(double fx, double fy, double fz);

    /**
     * @brief Ajoute une force vectorielle.
     *
     * @param f Force à ajouter.
     */
    void ajouterForce(const vecteur& f);
};

/**
 * @brief Affiche les principales caractéristiques d'une particule.
 *
 * @param os Flux de sortie.
 * @param p Particule à afficher.
 * @return Flux de sortie modifié.
 */
std::ostream& operator<<(std::ostream& os, const particule& p);

#endif
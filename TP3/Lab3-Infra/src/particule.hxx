#ifndef PARTICULE_HXX
#define PARTICULE_HXX
#include "vecteur.hxx"

/**
 * @brief Représente une particule dans l'espace
 * 
 * Une particule possède une position, une vitesse, une force,
 * ainsi qu'une masse, un identifiant et un type.
 */
class particule {

public:
    vecteur position; ///< Position de la particule
    vecteur vitesse;  ///< Vitesse de la particule
    vecteur force;    ///< Force appliquée sur la particule

    double masse;     ///< Masse de la particule
    int id;           ///< Identifiant unique
    int type;         ///< Type de la particule

    /// Constructeur par défaut
    particule();

    /**
     * @brief Constructeur avec paramètres
     * 
     * @param id identifiant de la particule
     * @param type type de la particule
     * @param m masse
     * @param p position initiale
     * @param v vitesse initiale
     */
    particule(int id,int type,double m,vecteur p,vecteur v);

    // getters

    /// Retourne la position
    const vecteur& getPosition() const;

    /// Retourne la vitesse
    const vecteur& getVitesse() const;

    /// Retourne la force
    const vecteur& getForce() const;

    /// Retourne la masse
    double getMasse() const;

    // setters

    /// Modifie la force
    void setForce(const vecteur& f);

    /// Modifie la position
    void setPosition(const vecteur& p);

    /// Modifie la vitesse
    void setVitesse(const vecteur& v);

    /**
     * @brief Fait évoluer la particule dans le temps
     * 
     * Met à jour la position et la vitesse en fonction de la force.
     * 
     * @param dt pas de temps
     */
    void evolue(double dt);

    /**
     * @brief Ajoute une force à la particule
     * 
     * @param f force à ajouter
     */
    void ajouterForce(const vecteur& f) { this->force = this->force + f; }
};

/**
 * @brief Affichage d'une particule
 * 
 * @param os flux de sortie
 * @param p particule à afficher
 * @return flux modifié
 */
std::ostream& operator<<(std::ostream& os, const particule& p);

#endif // PARTICULE_HXX
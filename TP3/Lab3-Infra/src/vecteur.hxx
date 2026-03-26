#ifndef VECTEUR_HXX
#define VECTEUR_HXX
#include <ostream>

/**
 * @brief Représente un vecteur en 3 dimensions
 * 
 * Permet de manipuler des vecteurs (addition, soustraction,
 * multiplication par un scalaire, norme, etc.).
 */
class vecteur {
private:
    double x; ///< Coordonnée x
    double y; ///< Coordonnée y
    double z; ///< Coordonnée z

public:
    /// Constructeur par défaut (vecteur nul)
    vecteur();

    /**
     * @brief Constructeur avec coordonnées
     * 
     * @param x_val coordonnée x
     * @param y_val coordonnée y
     * @param z_val coordonnée z
     */
    vecteur(double x_val, double y_val, double z_val);

    /**
     * @brief Calcule la norme du vecteur
     * 
     * @return longueur du vecteur
     */
    double norme() const;

    /**
     * @brief Addition de deux vecteurs
     */
    vecteur operator+(const vecteur& other) const;

    /**
     * @brief Soustraction de deux vecteurs
     */
    vecteur operator-(const vecteur& other) const;

    /**
     * @brief Multiplication par un scalaire
     */
    vecteur operator*(double scalar) const;

    /**
     * @brief Division par un scalaire
     */
    vecteur operator/(double scalar) const;

    /// Retourne la coordonnée x
    double getX() const;

    /// Retourne la coordonnée y
    double getY() const;

    /// Retourne la coordonnée z
    double getZ() const;
};

/**
 * @brief Affichage d'un vecteur
 * 
 * @param os flux de sortie
 * @param v vecteur à afficher
 * @return flux modifié
 */
std::ostream& operator<<(std::ostream& os, const vecteur& v);

#endif // VECTEUR_HXX
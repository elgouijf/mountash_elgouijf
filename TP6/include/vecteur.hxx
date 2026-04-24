#ifndef VECTEUR_HXX
#define VECTEUR_HXX

#include <ostream>

/**
 * @class vecteur
 * @brief Représente un vecteur de dimension 3.
 *
 * Cette classe est utilisée pour représenter les positions, vitesses et forces.
 * Même dans une simulation 1D ou 2D, les vecteurs restent stockés avec
 * trois composantes.
 */
class vecteur {
private:
    /**
     * @brief Composante selon x.
     */
    double x;

    /**
     * @brief Composante selon y.
     */
    double y;

    /**
     * @brief Composante selon z.
     */
    double z;

public:
    /**
     * @brief Constructeur avec initialisation des composantes.
     * @param x_val Valeur selon x.
     * @param y_val Valeur selon y.
     * @param z_val Valeur selon z.
     */
    vecteur(double x_val = 0.0, double y_val = 0.0, double z_val = 0.0);

    /**
     * @brief Calcule la norme euclidienne du vecteur.
     * @return La norme du vecteur.
     */
    double norme() const;

    /**
     * @brief Addition de deux vecteurs.
     * @param other Vecteur à additionner.
     * @return Le vecteur somme.
     */
    vecteur operator+(const vecteur& other) const;

    /**
     * @brief Soustraction de deux vecteurs.
     * @param other Vecteur à soustraire.
     * @return Le vecteur différence.
     */
    vecteur operator-(const vecteur& other) const;

    /**
     * @brief Multiplication par un scalaire.
     * @param scalar Scalaire multiplicatif.
     * @return Le vecteur multiplié.
     */
    vecteur operator*(double scalar) const;

    /**
     * @brief Division par un scalaire.
     * @param scalar Scalaire diviseur.
     * @return Le vecteur divisé.
     */
    vecteur operator/(double scalar) const;

    /**
     * @brief Retourne la composante x.
     * @return La composante x.
     */
    inline double getX() const { return x; }

    /**
     * @brief Retourne la composante y.
     * @return La composante y.
     */
    inline double getY() const { return y; }

    /**
     * @brief Retourne la composante z.
     * @return La composante z.
     */
    inline double getZ() const { return z; }

    /**
     * @brief Modifie la composante x.
     * @param x_val Nouvelle valeur.
     */
    inline void setX(double x_val) { x = x_val; }

    /**
     * @brief Modifie la composante y.
     * @param y_val Nouvelle valeur.
     */
    inline void setY(double y_val) { y = y_val; }

    /**
     * @brief Modifie la composante z.
     * @param z_val Nouvelle valeur.
     */
    inline void setZ(double z_val) { z = z_val; }

    /**
     * @brief Ajoute des composantes au vecteur courant.
     * @param dx Incrément selon x.
     * @param dy Incrément selon y.
     * @param dz Incrément selon z.
     */
    inline void ajoute(double dx, double dy, double dz) {
        x += dx;
        y += dy;
        z += dz;
    }

    /**
     * @brief Calcule le carré de la norme.
     * @return La norme au carré.
     */
    inline double norme2() const {
        return x*x + y*y + z*z;
    }

    /**
     * @brief Accès à une composante par indice.
     * @param index Indice de la composante.
     * @return Référence vers la composante demandée.
     */
    double& operator[](int index);

    /**
     * @brief Accès constant à une composante par indice.
     * @param index Indice de la composante.
     * @return Référence constante vers la composante demandée.
     */
    const double& operator[](int index) const;

    /**
 * @brief Constructeur de copie par défaut.
 */
    vecteur(const vecteur&) = default;

    /**
    * @brief Opérateur d'affectation par défaut.
    */
    vecteur& operator=(const vecteur&) = default;

    /**
    * @brief Destructeur par défaut.
    */
    ~vecteur() = default;
};

/**
 * @brief Affichage textuel d'un vecteur.
 * @param os Flux de sortie.
 * @param v Vecteur à afficher.
 * @return Le flux de sortie.
 */
std::ostream& operator<<(std::ostream& os, const vecteur& v);

#endif // VECTEUR_HXX
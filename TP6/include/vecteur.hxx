#ifndef VECTEUR_HXX
#define VECTEUR_HXX

#include <cmath>
#include <cassert>
#include <stdexcept>
#include <ostream>

class vecteur {

private:
    double x;
    double y;
    double z;

public:

    /**
     * @brief Construit un vecteur à partir de ses composantes.
     *
     * @param x_val Composante selon x.
     * @param y_val Composante selon y.
     * @param z_val Composante selon z.
     */
    inline vecteur(double x_val = 0.0,
                   double y_val = 0.0,
                   double z_val = 0.0) noexcept
        : x(x_val), y(y_val), z(z_val) {}

    /**
     * @brief Retourne la composante x.
     *
     * @return Valeur de x.
     */
    inline double getX() const noexcept {
        return x;
    }

    /**
     * @brief Retourne la composante y.
     *
     * @return Valeur de y.
     */
    inline double getY() const noexcept {
        return y;
    }

    /**
     * @brief Retourne la composante z.
     *
     * @return Valeur de z.
     */
    inline double getZ() const noexcept {
        return z;
    }

    /**
     * @brief Modifie la composante x.
     *
     * @param value Nouvelle valeur.
     */
    inline void setX(double value) noexcept {
        x = value;
    }

    /**
     * @brief Modifie la composante y.
     *
     * @param value Nouvelle valeur.
     */
    inline void setY(double value) noexcept {
        y = value;
    }

    /**
     * @brief Modifie la composante z.
     *
     * @param value Nouvelle valeur.
     */
    inline void setZ(double value) noexcept {
        z = value;
    }

    /**
     * @brief Calcule la norme euclidienne au carré.
     *
     * @return Valeur de x²+y²+z².
     */
    inline double norme2() const noexcept {
        return x * x + y * y + z * z;
    }

    /**
     * @brief Calcule la norme euclidienne du vecteur.
     *
     * @return Norme du vecteur.
     */
    inline double norme() const noexcept {
        return std::sqrt(norme2());
    }

    /**
     * @brief Addition de deux vecteurs.
     *
     * @param other Vecteur à ajouter.
     * @return Nouveau vecteur somme.
     */
    inline vecteur operator+(const vecteur& other) const noexcept {
        return vecteur(x + other.x,
                       y + other.y,
                       z + other.z);
    }

    /**
     * @brief Soustraction de deux vecteurs.
     *
     * @param other Vecteur à soustraire.
     * @return Nouveau vecteur différence.
     */
    inline vecteur operator-(const vecteur& other) const noexcept {
        return vecteur(x - other.x,
                       y - other.y,
                       z - other.z);
    }

    /**
     * @brief Multiplication par un scalaire.
     *
     * @param scalar Facteur multiplicatif.
     * @return Nouveau vecteur résultant.
     */
    inline vecteur operator*(double scalar) const noexcept {
        return vecteur(x * scalar,
                       y * scalar,
                       z * scalar);
    }

    /**
     * @brief Division par un scalaire.
     *
     * @param scalar Diviseur.
     * @return Nouveau vecteur résultant.
     */
    inline vecteur operator/(double scalar) const noexcept {
        return vecteur(x / scalar,
                       y / scalar,
                       z / scalar);
    }

    /**
     * @brief Accès en lecture à une composante.
     *
     * @param index Indice de composante.
     * @return Référence constante vers la composante.
     */
    inline const double& operator[](int index) const {
        if (index == 0) return x;
        if (index == 1) return y;
        if (index == 2) return z;

        throw std::out_of_range("Index out of range for vecteur");
    }

    /**
     * @brief Accès en écriture à une composante.
     *
     * @param index Indice de composante.
     * @return Référence vers la composante.
     */
    inline double& operator[](int index) noexcept {
        assert(index >= 0 && index < 3);

        if (index == 0) return x;
        if (index == 1) return y;

        return z;
    }

    /**
    * @brief Ajoute directement des composantes au vecteur courant.
    *
    * Cette méthode modifie le vecteur en place et évite la création
    * d'un vecteur temporaire.
    *
    * @param dx Incrément selon x.
    * @param dy Incrément selon y.
    * @param dz Incrément selon z.
    */
    inline void ajoute(double dx, double dy, double dz) noexcept {
        x += dx;
        y += dy;
        z += dz;
    }
};

/**
 * @brief Affichage d'un vecteur.
 *
 * @param os Flux de sortie.
 * @param v Vecteur à afficher.
 * @return Flux modifié.
 */
std::ostream& operator<<(std::ostream& os, const vecteur& v);

#endif
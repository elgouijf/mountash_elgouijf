#include "vecteur.hxx"
#include <cmath>
#include <iostream>
#include <cassert>


/**
 * @brief Construit un vecteur à partir de ses composantes, avec (0,0,0) par défaut.
 *
 * @param x_val Composante selon x.
 * @param y_val Composante selon y.
 * @param z_val Composante selon z.
 */
vecteur::vecteur(double x_val, double y_val, double z_val){
    x = x_val;
    y = y_val;
    z = z_val;
}

/**
 * @brief Calcule la norme euclidienne du vecteur.
 *
 * @return La norme du vecteur.
 */
double vecteur::norme() const{
    return std::sqrt(x*x + y*y + z*z);
}

/**
 * @brief Addition de deux vecteurs.
 *
 * @param other Vecteur à ajouter.
 * @return Nouveau vecteur résultant de l'addition.
 */
vecteur vecteur::operator+(const vecteur& other) const{
    return vecteur(x + other.x, y + other.y, z + other.z);
}

/**
 * @brief Soustraction de deux vecteurs.
 *
 * @param other Vecteur à soustraire.
 * @return Nouveau vecteur résultant de la soustraction.
 */
vecteur vecteur::operator-(const vecteur& other) const {
    return vecteur(x - other.x, y - other.y, z - other.z);
}

/**
 * @brief Multiplication par un scalaire.
 *
 * @param scalar Facteur multiplicatif.
 * @return Nouveau vecteur résultant de la multiplication.
 */
vecteur vecteur::operator*(double scalar) const{
    return vecteur(x * scalar, y * scalar, z * scalar);
}

/**
 * @brief Division par un scalaire.
 *
 * @param scalar Diviseur.
 * @return Nouveau vecteur résultant de la division.
 */
vecteur vecteur::operator/(double scalar) const{
    return vecteur(x / scalar, y / scalar, z / scalar);
}

/**
 * @brief Accès en lecture à une composante du vecteur.
 *
 * @param index Indice (0 = x, 1 = y, 2 = z).
 * @return Référence constante vers la composante.
 * @throws std::out_of_range si l'indice est invalide.
 */
const double& vecteur::operator[](int index) const {
    if (index == 0) return x;
    else if (index == 1) return y;
    else if (index == 2) return z;
    else throw std::out_of_range("Index out of range for vecteur");
}

/**
 * @brief Accès en écriture à une composante du vecteur.
 *
 * Permet de modifier directement une composante.
 *
 * @param index Indice (0 = x, 1 = y, 2 = z).
 * @return Référence vers la composante.
 * @throws std::out_of_range si l'indice est invalide.
 */
double& vecteur::operator[](int index) {
    assert(index >= 0 && index < 3);
    if (index == 0) return x;
    if (index == 1) return y;
    return z;
}

/**
 * @brief Affiche un vecteur sous forme lisible.
 *
 * @param os Flux de sortie.
 * @param v Vecteur à afficher.
 * @return Le flux de sortie.
 */
std::ostream& operator<<(std::ostream& os, const vecteur& v){
    os << "Vecteur(" << v.getX() << ", " << v.getY() << ", " << v.getZ() << ")";
    return os;
}
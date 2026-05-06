#include "vecteur.hxx"
#include <iostream>

/**
 * @brief Affiche un vecteur sous forme lisible.
 *
 * @param os Flux de sortie.
 * @param v Vecteur à afficher.
 * @return Le flux de sortie.
 */
std::ostream& operator<<(std::ostream& os, const vecteur& v) {
    os << "Vecteur("
       << v.getX() << ", "
       << v.getY() << ", "
       << v.getZ() << ")";

    return os;
}
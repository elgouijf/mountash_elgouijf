#ifndef CELLULE_HXX
#define CELLULE_HXX

#include "particule.hxx"
#include <vector>
#include <iostream>

/**
 * @class cellule
 * @brief Représente une cellule de la grille spatiale utilisée pour accélérer le calcul des interactions.
 *
 * Une cellule contient un ensemble de particules ainsi qu'une liste de cellules voisines.
 * Cette structure permet de limiter le calcul des interactions aux particules
 * situées dans la cellule courante et dans les cellules adjacentes.
 */
class cellule {
    private:
        /**
         * @brief Particules actuellement contenues dans la cellule.
         */
        std::vector<particule*> particules;

        /**
         * @brief Nombre de particules dans la cellule.
         */
        int num_particules;

        /**
         * @brief Liste des cellules voisines.
         */
        std::vector<cellule*> voisins;

    public:
        /**
         * @brief Constructeur par défaut.
         */
        cellule();

        /**
         * @brief Constructeur à partir d'une liste de particules et de voisins.
         * @param v Vecteur des particules contenues dans la cellule.
         * @param voisins Vecteur des cellules voisines.
         */
        cellule(std::vector<particule*>& v, std::vector<cellule*>& voisins);

        /**
         * @brief Ajoute une particule dans la cellule.
         * @param p Pointeur vers la particule à ajouter.
         */
        void ajoute_particule(particule* p);

        /**
         * @brief Vide la cellule de toutes ses particules.
         */
        void vide();

        /**
         * @brief Ajoute une cellule voisine à la liste des voisins.
         * @param voisin Pointeur vers la cellule voisine à ajouter.
         */
        void ajoute_voisin(cellule* voisin);

        /**
         * @brief Retourne les particules contenues dans la cellule.
         * @return Référence constante vers le vecteur de particules.
         */
        const std::vector<particule*>& getParticules() const;


        /**
         * @brief Retourne la liste des cellules voisines.
         * @return Référence constante vers le vecteur des voisins.
         */
        const std::vector<cellule*>& getVoisins() const;
};

#endif
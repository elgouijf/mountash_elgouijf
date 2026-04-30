#ifndef UNIVERS_HXX
#define UNIVERS_HXX

#include "particule.hxx"
#include <vector>
#include "cellule.hxx"
#include <iostream>
#include "conditionLimite.hxx"



/**
 * @class univers
 * @brief Représente l'ensemble du système physique simulé.
 *
 * L'univers contient toutes les particules ainsi qu'une grille de cellules
 * permettant d'accélérer le calcul des interactions de Lennard-Jones.
 *
 * Il stocke également les paramètres physiques et géométriques du problème :
 * - distance de coupure,
 * - dimension,
 * - tailles du domaine,
 * - nombre de cellules par direction,
 * - paramètres du potentiel de Lennard-Jones.
 */
class univers {
    private:
        /**
         * @brief Ensemble des particules du système.
         */
        std::vector<particule*> particules;

        /**
         * @brief Nombre total de particules.
         */
        int num_particules;

        /**
         * @brief Distance de coupure des interactions.
         */
        double r_cut;

        double G;

        /**
         * @brief Dimension de l'univers (1, 2 ou 3).
         */
        int dim;

        /**
         * @brief Tailles du domaine selon chaque direction.
         */
        std::vector<double> Lds;

        /**
         * @brief Nombre de cellules dans chaque direction.
         */
        std::vector<int> ncd;

        /**
         * @brief Grille des cellules de l'univers.
         */
        std::vector<cellule> cellules;

        /**
         * @brief Paramètre epsilon du potentiel de Lennard-Jones.
         */
        double eps;

        /**
         * @brief Paramètre sigma du potentiel de Lennard-Jones.
         */
        double sigma;

        ConditionLimite condl_xmin;

        ConditionLimite condl_xmax;

        ConditionLimite condl_ymin;

        ConditionLimite condl_ymax;

        ConditionLimite condl_zmin;

        ConditionLimite condl_zmax;

        bool utiliser_potentiel_mur;

    public:
        /**
         * @brief Constructeur par défaut.
         */
        univers();

        univers(std::vector<particule*>& v,
                std::vector<double> Lds,
                double r_cut,
                int dim,
                double eps,
                double sigma);

        univers(std::vector<particule*>& v,
                std::vector<double> Lds,
                double r_cut,
                int dim,
                double eps,
                double sigma,
                double G);

        univers(std::vector<particule*>& v,
                std::vector<double> Lds,
                double r_cut,
                int dim,
                double eps,
                double sigma,
                double G,
                bool utiliser_potentiel_mur);

        /**
         * @brief Destructeur.
         *
         * Libère la mémoire des particules allouées dynamiquement.
         */
        ~univers();

        /**
         * @brief Ajoute une particule à l'univers.
         * @param p Pointeur vers la particule à ajouter.
         */
        void ajoute_particule(particule* p);

        /**
         * @brief Fait évoluer toutes les particules sur un pas de temps.
         * @param dt Pas de temps.
         */
        void evolue_particules(double dt);

        /**
         * @brief Calcule les forces exercées sur toutes les particules.
         *
         * Les interactions sont limitées aux particules situées dans la cellule
         * courante et dans ses voisines, avec une coupure à la distance r_cut.
         */
        void calcule_forces();

        /**
         * @brief Vide toutes les cellules de l'univers.
         */
        void vide_cellules();

        /**
         * @brief Place une particule dans la cellule correspondant à sa position.
         * @param p Pointeur vers la particule à placer.
         */
        void place_particule_dans_cellule(particule* p);

        /**
         * @brief Replace toutes les particules dans la grille de cellules.
         */
        void place_particules_dans_cellules();

        /**
         * @brief Initialise la grille de cellules et leurs relations de voisinage.
         */
        void initialise_cellules();

        /**
         * @brief Calcule l'indice d'une cellule dans le tableau linéarisé.
         * @param ix Indice selon x.
         * @param iy Indice selon y.
         * @param iz Indice selon z.
         * @return L'indice dans le vecteur des cellules.
         */
        int indice_cellule(int ix, int iy = 0, int iz = 0) const;

        /**
         * @brief Retourne l'ensemble des particules.
         * @return Référence constante vers le vecteur de particules.
         */
        const std::vector<particule*>& getParticules() const;

        /**
         * @brief Retourne le nombre de particules.
         * @return Nombre de particules.
         */
        int getNumParticules() const;

        /**
         * @brief Retourne la distance de coupure.
         * @return Valeur de r_cut.
         */
        double getRCut() const;

        /**
         * @brief Retourne la dimension de l'univers.
         * @return Dimension.
         */
        int getDim() const;

        /**
         * @brief Retourne les tailles du domaine.
         * @return Référence constante vers le vecteur des tailles.
         */
        const std::vector<double>& getLds() const;

        /**
         * @brief Retourne le nombre de cellules par direction.
         * @return Référence constante vers le vecteur ncd.
         */
        const std::vector<int>& getNcd() const;

        /**
         * @brief Retourne la grille des cellules.
         * @return Référence constante vers le vecteur des cellules.
         */
        const std::vector<cellule>& getCellules() const;

        /**
         * @brief Retourne la grille des cellules.
         * @return Référence modifiable vers le vecteur des cellules.
         */
        std::vector<cellule>& getCellules();

        /**
         * @brief Modifie le paramètre epsilon.
         * @param eps Nouvelle valeur de epsilon.
         */
        void setEps(double eps);

        /**
         * @brief Modifie le paramètre sigma.
         * @param sigma Nouvelle valeur de sigma.
         */
        void setSigma(double sigma);

        double energie_cinetique() const;

        void setConditionsLimites(ConditionLimite xmin, ConditionLimite xmax,
                              ConditionLimite ymin, ConditionLimite ymax,
                              ConditionLimite zmin, ConditionLimite zmax);

        void applique_conditions_limites();

        bool applique_conditions_limites_particule(particule* p);

        double calcule_force_mur(double r) const;

        void limite_vitesses(int N1,int N2);

        void setUtiliserPotentielMur(bool actif);

        void applique_gravite();

        void applique_potentiel_mur();
};

/**
 * @brief Affichage textuel de l'univers.
 * @param os Flux de sortie.
 * @param u Univers à afficher.
 * @return Le flux de sortie.
 */
std::ostream& operator<<(std::ostream& os, const univers& u);

#endif // UNIVERS_HXX
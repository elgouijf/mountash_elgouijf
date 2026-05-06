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
         * @brief Vecteur des cellules occupées.
         */
        std::vector<cellule*> cellules_occupees;

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

        /** @brief Vide toutes les particules de l'univers. */
        void clear_particules();

        /** @brief Copie les paramètres d'un autre univers. */
        void copier_parametres_depuis(const univers& other);

        /** @brief Indique si aucune condition limite n'est appliquée. */
        bool aucune_cond_limite = false;

        std::vector<double> omp_fx;
        std::vector<double> omp_fy;
        std::vector<double> omp_fz;

        int omp_threads_alloc = 0;
        int omp_particles_alloc = 0;

        void prepare_omp_force_buffers();

        

    public:
        /**
         * @brief Constructeur par défaut.
         */
        univers();

        /**
         * @brief Construit un univers paramétré.
         *
         * Initialise les paramètres physiques et géométriques de l'univers,
         * construit la grille de cellules, puis y ajoute les particules fournies.
         *
         * @param v Vecteur initial de particules.
         * @param Lds Tailles du domaine selon chaque direction.
         * @param r_cut Distance de coupure des interactions.
         * @param dim Dimension de l'univers.
         * @param eps Paramètre epsilon du potentiel de Lennard-Jones.
         * @param sigma Paramètre sigma du potentiel de Lennard-Jones.
         */
        univers(std::vector<particule*>& v,
                std::vector<double> Lds,
                double r_cut,
                int dim,
                double eps,
                double sigma);
        
        
        /** @brief Construit un univers paramétré avec un champ gravitationnel. */
        univers(std::vector<particule*>& v,
                std::vector<double> Lds,
                double r_cut,
                int dim,
                double eps,
                double sigma,
                double G);
        
        /** @brief Construit un univers paramétré avec un champ gravitationnel et potentiel de mur. */
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
        * @brief Constructeur de copie.
        * @param other L'objet à copier.
        */
        univers(const univers& other);

        /**
         * @brief Opérateur d'affectation par copie.
         * @param other L'objet à copier.
         * @return Référence à l'objet copié.
         */
        univers& operator=(const univers& other);

        /** 
         * @brief Constructeur de déplacement.
         * @param other L'objet à déplacer.
         */
        univers(univers&& other) noexcept;

        /** 
         * @brief Opérateur d'affectation par déplacement.
         * @param other L'objet à déplacer.
         * @return Référence à l'objet déplacé.
         */
        univers& operator=(univers&& other) noexcept;

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

        void calcule_forces_sequentiel();

        /**
         * @brief Vide toutes les cellules de l'univers.
         */
        void vide_cellules();

        /**
         * @brief Place une particule dans la cellule correspondant à sa position.
         * @param p Pointeur vers la particule à placer.
         * @return Pointeur vers la cellule dans laquelle la particule a été placée.
         */
        cellule* place_particule_dans_cellule(particule* p);

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
         * @brief Modifie le paramètre epsilon.
         * @param eps Nouvelle valeur de epsilon.
         */
        void setEps(double eps);

        /**
         * @brief Modifie le paramètre sigma.
         * @param sigma Nouvelle valeur de sigma.
         */
        void setSigma(double sigma);

        /**
         * @brief Calcule l'énergie cinétique totale de l'univers.
         * @return Énergie cinétique.
         */
        double energie_cinetique() const;
        
        /**
         * @brief Calcule l'énergie potentielle totale de l'univers.
         * @return Énergie potentielle.
         */
        double energie_potentielle() const;

        /**
         * @brief Calcule l'énergie mécanique totale de l'univers.
         * @return Énergie mécanique.
         */
        double energie_mecanique() const;

        /**
        * @brief Modifie les conditions aux limites selon chaque direction.
        * @param xmin Condition aux limites pour x_min.
        * @param xmax Condition aux limites pour x_max.
        * @param ymin Condition aux limites pour y_min.
        * @param ymax Condition aux limites pour y_max.
        * @param zmin Condition aux limites pour z_min.
        * @param zmax Condition aux limites pour z_max.
        */
        void setConditionsLimites(ConditionLimite xmin, ConditionLimite xmax,
                              ConditionLimite ymin, ConditionLimite ymax,
                              ConditionLimite zmin, ConditionLimite zmax);


        /**
         * @brief Applique les conditions aux limites à toutes les particules.
         * @return true si au moins une particule a été supprimée, false sinon.
         */
        void applique_conditions_limites();
        
        /** @brief Applique les conditions aux limites à une particule donnée.
         * @param p Pointeur vers la particule.
         * @return true si la particule a été supprimée, false sinon.
        */
        bool applique_conditions_limites_particule(particule* p);

        /** @brief Calcule la force exercée par un mur sur une particule.
         * @param r Distance entre la particule et le mur.
         * @return Force exercée.
         */
        double calcule_force_mur(double r) const;
        

        /** @brief Limite les vitesses des particules.
         * @param N1 Index de la première particule.
         * @param N2 Index de la deuxième particule.
         */
        void limite_vitesses(int N1,int N2);

        /** @brief Active ou désactive l'utilisation du potentiel de mur.
         * @param actif true pour activer, false pour désactiver.
         */
        void setUtiliserPotentielMur(bool actif);

        /** @brief Applique la gravité à toutes les particules.
         */
        void applique_gravite();

        /** @brief Applique le potentiel de mur à toutes les particules.
         */
        void applique_potentiel_mur();

        /** @brief Affiche les conditions aux limites.
         */
        void afficherConditionsLimites() const;

        void debug_cellules() const;

        void reserveParticules(size_t n);
};

/**
 * @brief Affichage textuel de l'univers.
 * @param os Flux de sortie.
 * @param u Univers à afficher.
 * @return Le flux de sortie.
 */
std::ostream& operator<<(std::ostream& os, const univers& u);

#endif // UNIVERS_HXX
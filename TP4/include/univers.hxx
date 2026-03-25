#ifndef UNIVERS_HXX
#define UNIVERS_HXX

#include "particule.hxx"
#include <vector>
#include "../include/cellule.hxx"
#include <iostream> // pour std::ostream

class univers {
    private:
        std::vector<particule*> particules;
        int num_particules;
        double r_cut;
        int dim;
        std::vector<double> Lds; // tailles de domaine pour chaque dimension
        std::vector<int> ncd; // nombre de cellules par dimension
        std::vector<cellule> cellules; // grille de cellules, je met pas de pointeurs car il seront réer par univers de toute façon et ça me sera plus facile pour les gérer
        // autrement, pour moi les cellule est une partie de l'univers
        double eps;
        double sigma;
        
    public:
        univers();
        univers(std::vector<particule*>& v, std::vector<double> Lds, double r_cut, int dim, double eps, double sigma);
        ~univers();

        
        void ajoute_particule(particule* p);
        void evolue_particules(double dt);
        void calcule_forces();
        void vide_cellules();
        void place_particule_dans_cellule(particule* p);
        void place_particules_dans_cellules();
        void initialise_cellules();
        int indice_cellule(int ix, int iy = 0, int iz = 0) const;
        // getters
        const std::vector<particule*>& getParticules() const;
        int getNumParticules() const;
        double getRCut() const;
        int getDim() const;
        const std::vector<double>& getLds() const;
        const std::vector<int>& getNcd() const;
        const std::vector<cellule>& getCellules() const;
        std::vector<cellule>& getCellules(); // version non const pour pouvoir modifié les cellules dans l'
        
        // setters
        void setEps(double eps);
        void setSigma(double sigma);


        
};

std::ostream& operator<<(std::ostream& os, const univers& u);
#endif // PARTICULE_HXX
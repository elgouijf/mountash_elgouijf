#ifndef PARTICULE_HXX
#define PARTICULE_HXX
#include "vecteur.hxx"
#include <iostream>
#include <vector>


class particule {
    private:
        vecteur position;
        vecteur vitesse;
        vecteur force;
        vecteur force_old;
        //std::vector<vecteur> voisins;


        double masse;
        int id;
        int type;

        
    public:

        particule();
        particule(int id,int type,double m,vecteur p,vecteur v);

        // getters
        int getId() const;
        int getType() const;
        const vecteur& getPosition() const; // 
        const vecteur& getVitesse() const;
        const vecteur& getForce() const;
        const vecteur& getForceOld() const;
        double getMasse() const;


        // setters
        void setForce(const vecteur& f);
        void setPosition(const vecteur& p);
        void setVitesse(const vecteur& v);
        void setForceOld(const vecteur& f);


        void evolue(double dt);
        void ajouterForce(const vecteur& f);
        void ajouterForce(double fx, double fy, double fz); // Version allégee : on évite de créer un vecteur temporaire pour la force, et on ajoute directement les composantes à la force actuelle
        void avance_position_verlet(double dt);
        void avance_vitesse_verlet(double dt);

    };

std::ostream& operator<<(std::ostream& os, const particule& p);
#endif // PARTICULE_HXX

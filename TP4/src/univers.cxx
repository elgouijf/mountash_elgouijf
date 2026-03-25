#include "../include/univers.hxx"
#include <cmath>
univers::univers(){
     // vecteur vide
    this->num_particules = 0;
    //this->particules = std::vector<particule*>(); c'est déjà vide
    this->dim = 3; // exemple de dimension
    this->Lds = std::vector<double>(this->dim, 10.0);
    this->r_cut = 2.5; // exemple de distance de coupure
    this->ncd = std::vector<int>(this->dim); // exemple de nombre de cellules par dimension
    this->eps = 1.0;
    this->sigma = 1.0;
    for (int i = 0; i < this->dim; ++i) {
        this->ncd[i] = static_cast<int>(this->Lds[i] / this->r_cut);
        if (this->ncd[i] < 1) this->ncd[i] = 1; // au moins une cellule par dimension
    }

    initialise_cellules();

}

univers::univers(std::vector<particule*>& v,
    std::vector<double> Lds,
    double r_cut,
    int dim,
    double eps,
    double sigma
){
    this->num_particules = 0;

    this->r_cut = r_cut;
    this->eps = eps;
    this->sigma = sigma;

    if (dim < 1 || dim > 3) {
        std::cerr << "Dimension must be between 1 and 3. Setting to 3." << std::endl;
        this->dim = 3;
    } else {
        this->dim = dim;
    }

    if (Lds.size() != this->dim) {
        std::cerr << "Size of Lds must match the dimension. Setting to default values." << std::endl;
        this->Lds = std::vector<double>(this->dim, 10.0);
    } else {
        this->Lds = Lds;
    }

    this->ncd = std::vector<int>(this->dim);
    for (int i = 0; i < this->dim; ++i) {
        this->ncd[i] = static_cast<int>(this->Lds[i] / r_cut);
        if (this->ncd[i] < 1) this->ncd[i] = 1;
    }

    initialise_cellules();

    this->particules.reserve(v.size());
    for (particule* p : v) {
        ajoute_particule(p);
    }
}

univers::~univers() {
    for (particule* p : particules) {
        delete p;
    }
    particules.clear();
}


const std::vector<particule*>& univers::getParticules() const {
    return particules;
}

int univers::getNumParticules() const {
    return num_particules;
}

double univers::getRCut() const {
    return r_cut;
}

int univers::getDim() const {
    return dim;
}

const std::vector<double>& univers::getLds() const {
    return Lds;
}

const std::vector<int>& univers::getNcd() const {
    return ncd;
}

const std::vector<cellule>& univers::getCellules() const {
    return cellules;
}

std::vector<cellule>& univers::getCellules() {
    return cellules;
}


void univers::ajoute_particule(particule* p) {
    vecteur pos = p->getPosition();
    vecteur vit = p->getVitesse();
    vecteur frc = p->getForce();

    if (dim == 1) {
        pos.setY(0.0); pos.setZ(0.0);
        vit.setY(0.0); vit.setZ(0.0);
        frc.setY(0.0); frc.setZ(0.0);

    } 
    else if (dim == 2) {
        pos.setZ(0.0); vit.setZ(0.0); frc.setZ(0.0);
    }

    p->setPosition(pos);
    p->setVitesse(vit);
    p->setForce(frc);

    this->particules.push_back(p);
    this->num_particules++;

    place_particule_dans_cellule(p);
}

void univers::evolue_particules(double dt) {
    for (particule* p : particules) {
        p->avance_position_verlet(dt);
    }
    place_particules_dans_cellules();
    calcule_forces();

    for (particule* p : particules) {
        p->avance_vitesse_verlet(dt);
    }
}

void univers::vide_cellules() {
    for (cellule& c : cellules) {
        c.vide();
    }
}

void univers::place_particule_dans_cellule(particule* p) {
    const vecteur& pos = p->getPosition();

    int ix = static_cast<int>(pos.getX() / r_cut);
    if (ix < 0) ix = 0;
    if (ix >= ncd[0]) ix = ncd[0] - 1;

    int iy = 0;
    int iz = 0;

    if (dim >= 2) {
        iy = static_cast<int>(pos.getY() / r_cut);
        if (iy < 0) iy = 0;
        if (iy >= ncd[1]) iy = ncd[1] - 1;
    }

    if (dim == 3) {
        iz = static_cast<int>(pos.getZ() / r_cut);
        if (iz < 0) iz = 0;
        if (iz >= ncd[2]) iz = ncd[2] - 1;
    }

    int idx = indice_cellule(ix, iy, iz);
    cellules[idx].ajoute_particule(p);
}

void univers::place_particules_dans_cellules() {
    vide_cellules();

    for (particule* p : particules) {
        place_particule_dans_cellule(p);
    }
}



int univers::indice_cellule(int ix, int iy, int iz) const {
    if (dim == 1) {
        return ix;
    }
    if (dim == 2) {
        return ix + ncd[0] * iy;
    }
    return ix + ncd[0] * (iy + ncd[1] * iz);
}

void univers::initialise_cellules() {
    int nb_total = 1;
    for (int d = 0; d < dim; ++d) {
        nb_total *= this->ncd[d];
    }

    int nb_voisins_max = 1;
    for (int d = 0; d < dim; ++d) {
        nb_voisins_max *= 3;
    }

    cellules.clear();
    cellules.resize(nb_total); // on créer les cellules vides

    for (cellule& c : cellules) {
        c.getVoisins().reserve(nb_voisins_max);
    }

    if (dim == 1) {
        for (int ix = 0; ix < ncd[0]; ++ix) {
            int idx = indice_cellule(ix);

            for (int dx = -1; dx <= 1; ++dx) {
                int nix = ix + dx;
                if (nix >= 0 && nix < ncd[0]) {
                    int idx_voisin = indice_cellule(nix);
                    cellules[idx].ajoute_voisin(&cellules[idx_voisin]);
                }
            }
        }
    }

    else if (dim == 2) {
        for (int ix = 0; ix < ncd[0]; ++ix) {
            for (int iy = 0; iy < ncd[1]; ++iy) {
                int idx = indice_cellule(ix, iy);

                for (int dx = -1; dx <= 1; ++dx) {
                    for (int dy = -1; dy <= 1; ++dy) {
                        int nix = ix + dx;
                        int niy = iy + dy;

                        if (nix >= 0 && nix < ncd[0] &&
                            niy >= 0 && niy < ncd[1]) {
                            int idx_voisin = indice_cellule(nix, niy);
                            cellules[idx].ajoute_voisin(&cellules[idx_voisin]);
                        }
                    }
                }
            }
        }
    }

    else if (dim == 3) {
        for (int ix = 0; ix < ncd[0]; ++ix) {
            for (int iy = 0; iy < ncd[1]; ++iy) {
                for (int iz = 0; iz < ncd[2]; ++iz) {
                    int idx = indice_cellule(ix, iy, iz);

                    for (int dx = -1; dx <= 1; ++dx) {
                        for (int dy = -1; dy <= 1; ++dy) {
                            for (int dz = -1; dz <= 1; ++dz) {
                                int nix = ix + dx;
                                int niy = iy + dy;
                                int niz = iz + dz;

                                if (nix >= 0 && nix < ncd[0] &&
                                    niy >= 0 && niy < ncd[1] &&
                                    niz >= 0 && niz < ncd[2]) {
                                    int idx_voisin = indice_cellule(nix, niy, niz);
                                    cellules[idx].ajoute_voisin(&cellules[idx_voisin]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}




void univers::calcule_forces(){
    /* Version optimisée : 
        - on ne calcule que des carrés
        - pas de pow ni sqrt
        - pas de vecteurs inutiles, mais usage de dx, dy et dz
        - pas de projections */
    //  Remise à zéro de toutes les forces
    for (particule* p : particules) {
        p->setForce(vecteur());
    }

    // Placer les particules dans les cellules
    //place_particules_dans_cellules(); pas la peine car evolue met chcune des particules dans leurs bonnes endroits

    const double r_cut2 = r_cut * r_cut;
    const double sigma2 = sigma * sigma;
    const double coeff = 24.0 * this->eps;

    //  Calculer les forces en utilisant les cellules
    for (const cellule& c : this->cellules) {
        const auto& parts = c.getParticules();

        for (const cellule* v : c.getVoisins()){
            // il n' ya que les voisins qui peuvent etre soumises à l'effet des particules de cette cellule
            if (v < &c){
                // cellule déjà rencontré (car on parcour par ordre croissant de ponteurs)
                continue;
            }

            if (v == &c){
                // meme cellule, il faut calculer les interactions à l'intérieur
                // on a besoin de vérifier si dist < r_cut, car la taille de la cellule est rcut dans chaque direction, donc la diag par exemple peut etre trop grande
                
                for (size_t i = 0; i < parts.size(); ++i) {
                    for (size_t j = i + 1; j < parts.size(); ++j) {

                        particule* pi = parts[i];
                        particule* pj = parts[j];

                        const vecteur& posi = pi->getPosition();
                        const vecteur& posj = pj->getPosition();

                        double dx = posj.getX() - posi.getX();
                        double dy = posj.getY() - posi.getY();
                        double dz = posj.getZ() - posi.getZ();

                        double dist2 = dx * dx + dy * dy + dz * dz + 1e-12;

                        if (dist2 > r_cut2) continue;

                        double sr2 = sigma2 / dist2;
                        double sr6 = sr2 * sr2 * sr2;
                        double sr12 = sr6 * sr6;

                        //double coeff_force = coeff * (2.0 * sr12 - sr6) / dist2;
                        double coeff_force = coeff * sr6 * (1.0 - 2.0 * sr6) / dist2;

                        pi->ajouterForce(dx * coeff_force, dy * coeff_force, dz * coeff_force);
                        pj->ajouterForce(-dx * coeff_force, -dy * coeff_force, -dz * coeff_force);
                    }
                }
            }
            else {
                // Ici pas besoin de faire un parcours par idx, car les cellules sont différents don con aura jamais pi = pj
                const auto& vois = v->getParticules();
                
                for (particule* pi : parts) {
                    const vecteur& posi = pi->getPosition();

                    for (particule* pj : vois) {
                        const vecteur& posj = pj->getPosition();
                    
                        double dx = posj.getX() - posi.getX();
                        double dy = posj.getY() - posi.getY();
                        double dz = posj.getZ() - posi.getZ();

                        double dist2 = dx * dx + dy * dy + dz * dz + 1e-12;

                        if (dist2 > r_cut2) continue;

                        double sr2 = sigma2 / dist2;
                        double sr6 = sr2 * sr2 * sr2;
                        double sr12 = sr6 * sr6;

                        //double coeff_force = coeff * (2.0 * sr12 - sr6) / dist2;
                        double coeff_force = coeff * sr6 * (1.0 - 2.0 * sr6) / dist2;

                        pi->ajouterForce(dx * coeff_force, dy * coeff_force, dz * coeff_force);
                        pj->ajouterForce(-dx * coeff_force, -dy * coeff_force, -dz * coeff_force);
                    }
                }
            }
        }
    }
}

std::ostream& operator<<(std::ostream& os, const univers& u) {
    for (auto* p : u.getParticules()){
        os << *p << std::endl;
    }
    return os;
}

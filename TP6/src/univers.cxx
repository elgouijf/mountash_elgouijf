#include "../include/univers.hxx"
#include <cmath>

/**
 * @brief Construit un univers par défaut.
 *
 * Initialise un univers tridimensionnel avec :
 * - un domaine cubique de taille 10 dans chaque direction,
 * - une distance de coupure égale à 2.5,
 * - des paramètres de Lennard-Jones par défaut.
 *
 * La grille de cellules est ensuite initialisée.
 */
univers::univers(){
    this->num_particules = 0;
    this->dim = 3;
    this->Lds = std::vector<double>(this->dim, 10.0);
    this->r_cut = 2.5;
    this->ncd = std::vector<int>(this->dim);
    this->eps = 1.0;
    this->sigma = 1.0;

    for (int i = 0; i < this->dim; ++i) {
        this->ncd[i] = static_cast<int>(this->Lds[i] / this->r_cut);
        if (this->ncd[i] < 1) this->ncd[i] = 1;
    }

    // Par défaut, toutes les conditions aux limites sont réflexives
    this->condl_xmin = ConditionLimite::Reflexive;
    this->condl_xmax = ConditionLimite::Reflexive;
    this->condl_ymin = ConditionLimite::Reflexive;
    this->condl_ymax = ConditionLimite::Reflexive;
    this->condl_zmin = ConditionLimite::Reflexive;
    this->condl_zmax = ConditionLimite::Reflexive; 

    initialise_cellules();
}

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

/**
 * @brief Détruit l'univers et libère la mémoire des particules.
 */
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

/**
 * @brief Ajoute une particule à l'univers en respectant sa dimension effective.
 *
 * Pour un univers 1D ou 2D, les composantes inutiles de position, vitesse et force
 * sont annulées avant insertion. La particule est ensuite ajoutée à la bonne cellule.
 *
 * @param p Particule à ajouter.
 */
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

/**
 * @brief Fait évoluer l'ensemble des particules sur un pas de temps avec Störmer-Verlet.
 *
 * L'évolution se fait en quatre étapes :
 * 1. mise à jour des positions avec les forces actuelles,
 * 2. remise à jour de la grille de cellules,
 * 3. recalcul des forces,
 * 4. mise à jour des vitesses.
 *
 * @param dt Pas de temps.
 */
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


double univers::energie_cinetique() const {
    double Ec = 0.0;

    for (const particule* p : particules) {
        const vecteur& v = p->getVitesse();
        double v2 = v.getX()*v.getX() + v.getY()*v.getY() + v.getZ()*v.getZ();
        Ec += 0.5 * p->getMasse() * v2;
    }

    return Ec;
}

/**
 * @brief Vide toutes les cellules de la grille.
 *
 * Utilisé avant de replacer les particules dans les cellules correspondant
 * à leurs nouvelles positions.
 */
void univers::vide_cellules() {
    for (cellule& c : cellules) {
        c.vide();
    }
}

/**
 * @brief Place une particule dans la cellule correspondant à sa position.
 *
 * Les indices sont bornés aux limites de la grille afin d'éviter tout
 * débordement dans le tableau des cellules.
 *
 * @param p Particule à placer.
 */
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

/**
 * @brief Replace toutes les particules dans la grille de cellules.
 *
 * Toutes les cellules sont d'abord vidées, puis chaque particule
 * est replacée à partir de sa position actuelle.
 */
void univers::place_particules_dans_cellules() {
    vide_cellules();

    for (particule* p : particules) {
        place_particule_dans_cellule(p);
    }
}

/**
 * @brief Convertit des indices de cellule en indice linéaire dans le tableau.
 *
 * Cette fonction permet de stocker une grille 1D, 2D ou 3D
 * dans un vecteur linéaire unique.
 *
 * @param ix Indice selon x.
 * @param iy Indice selon y.
 * @param iz Indice selon z.
 * @return Indice de la cellule dans le vecteur cellules.
 */
int univers::indice_cellule(int ix, int iy, int iz) const {
    if (dim == 1) {
        return ix;
    }
    if (dim == 2) {
        return ix + ncd[0] * iy;
    }
    return ix + ncd[0] * (iy + ncd[1] * iz);
}

/**
 * @brief Initialise la grille de cellules et les relations de voisinage.
 *
 * Chaque cellule reçoit la liste de ses voisines immédiates.
 * Cela permet ensuite de limiter le calcul des interactions
 * aux particules proches dans l'espace.
 */
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
    cellules.resize(nb_total);

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

/**
 * @brief Définit les conditions aux limites sur chaque face du domaine.
 *
 * @param xmin Condition sur la face x = 0.
 * @param xmax Condition sur la face x = Lx.
 * @param ymin Condition sur la face y = 0.
 * @param ymax Condition sur la face y = Ly.
 * @param zmin Condition sur la face z = 0.
 * @param zmax Condition sur la face z = Lz.
 */
void univers::setConditionsLimites(ConditionLimite cond_xmin, ConditionLimite cond_xmax,
                                ConditionLimite cond_ymin, ConditionLimite cond_ymax,
                                ConditionLimite cond_zmin, ConditionLimite cond_zmax) {
    this->condl_xmin = cond_xmin;
    this->condl_xmax = cond_xmax;
    this->condl_ymin = cond_ymin;
    this->condl_ymax = cond_ymax;
    this->condl_zmin = cond_zmin;
    this->condl_zmax = cond_zmax;
}

/**
 * @brief Applique les conditions aux limites à une particule.
 *
 * @param p Particule à traiter.
 * @return false si la particule doit être absorbée, true sinon.
 */

bool univers::applique_conditions_limites_particule(particule* p){
    vecteur pos = p->getPosition();
    vecteur vit = p->getVitesse();

    // Conditions aux limites, bornes négatives

    // x
    if (pos.getX() < 0.0) {
        if (condl_xmin == ConditionLimite::Reflexive) {
            pos.setX(-pos.getX());
            vit.setX(-vit.getX());
        } else if (condl_xmin == ConditionLimite::Absorbante) {
            return false;
        } else if (condl_xmin == ConditionLimite::Periodique) {
            pos.setX(pos.getX() + Lds[0]);
        }

    }

    // y
    if (pos.getY() < 0.0) {
        if (condl_ymin == ConditionLimite::Reflexive) {
            pos.setY(-pos.getY());
            vit.setY(-vit.getY());
        } else if (condl_ymin == ConditionLimite::Absorbante) {
            return false;
        } else if (condl_ymin == ConditionLimite::Periodique) {
            pos.setY(pos.getY() + Lds[1]);
        }
    }

    // z
    if (pos.getZ() < 0.0) {
        if (condl_zmin == ConditionLimite::Reflexive) {
            pos.setZ(-pos.getZ());
            vit.setZ(-vit.getZ());
        } else if (condl_zmin == ConditionLimite::Absorbante) {
            return false;
        } else if (condl_zmin == ConditionLimite::Periodique) {
            pos.setZ(pos.getZ() + Lds[2]);
        }
    }

    // Conditions aux limites, bornes positives
    // x
    if (pos.getX() > Lds[0]) {
        if (condl_xmax == ConditionLimite::Reflexive) {
            pos.setX(2.0 * Lds[0] - pos.getX());
            vit.setX(-vit.getX());
        } else if (condl_xmax == ConditionLimite::Absorbante) {
            return false;
        } else if (condl_xmax == ConditionLimite::Periodique) {
            pos.setX(pos.getX() - Lds[0]);
        }
    }

    // y
    if (pos.getY() > Lds[1]) {
        if (condl_ymax == ConditionLimite::Reflexive) {
            pos.setY(2.0 * Lds[1] - pos.getY());
            vit.setY(-vit.getY());
        } else if (condl_ymax == ConditionLimite::Absorbante) {
            return false;
        } else if (condl_ymax == ConditionLimite::Periodique) {
            pos.setY(pos.getY() - Lds[1]);
        }
    }

    // z
    if (pos.getZ() > Lds[2]) {
        if (condl_zmax == ConditionLimite::Reflexive) {
            pos.setZ(2.0 * Lds[2] - pos.getZ());
            vit.setZ(-vit.getZ());
        } else if (condl_zmax == ConditionLimite::Absorbante) {
            return false;
        } else if (condl_zmax == ConditionLimite::Periodique) {
            pos.setZ(pos.getZ() - Lds[2]);
        }
    } 

    // Mettre à jour la position et la vitesse de la particule
    p->setPosition(pos);
    p->setVitesse(vit);

    return true;

}

/**
 * @brief Calcule les forces de Lennard-Jones sur toutes les particules.
 *
 * Cette implémentation exploite la décomposition spatiale en cellules afin de
 * limiter le calcul aux paires de particules potentiellement proches.
 *
 * Optimisations utilisées :
 * - calcul uniquement avec des distances au carré,
 * - absence de sqrt et de pow,
 * - calcul direct sur les composantes dx, dy, dz,
 * - application du principe action-réaction.
 *
 * Les interactions sont tronquées à la distance de coupure r_cut.
 *
 * La force utilisée correspond à une interaction de type Lennard-Jones.
 */
void univers::calcule_forces(){
    // Remise à zéro des forces
    for (particule* p : particules) {
        p->setForce(vecteur());
    }

    const double r_cut2 = r_cut * r_cut;
    const double sigma2 = sigma * sigma;
    const double coeff = 24.0 * this->eps;

    // Parcours des cellules
    for (const cellule& c : this->cellules) {
        const auto& parts = c.getParticules();

        for (const cellule* v : c.getVoisins()){
            // Évite de traiter deux fois la même paire de cellules
            if (v < &c){
                continue;
            }

            if (v == &c){
                // Interactions internes à une même cellule
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

                        double coeff_force = coeff * sr6 * (1.0 - 2.0 * sr6) / dist2;

                        pi->ajouterForce(dx * coeff_force, dy * coeff_force, dz * coeff_force);
                        pj->ajouterForce(-dx * coeff_force, -dy * coeff_force, -dz * coeff_force);
                    }
                }
            }
            else {
                // Interactions entre deux cellules voisines distinctes
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

                        double coeff_force = coeff * sr6 * (1.0 - 2.0 * sr6) / dist2;

                        pi->ajouterForce(dx * coeff_force, dy * coeff_force, dz * coeff_force);
                        pj->ajouterForce(-dx * coeff_force, -dy * coeff_force, -dz * coeff_force);
                    }
                }
            }
        }
    }
}

/**
 * @brief Affiche toutes les particules de l'univers.
 *
 * @param os Flux de sortie.
 * @param u Univers à afficher.
 * @return Le flux de sortie.
 */
std::ostream& operator<<(std::ostream& os, const univers& u) {
    for (auto* p : u.getParticules()){
        os << *p << std::endl;
    }
    return os;
}
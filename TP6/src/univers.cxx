#include "univers.hxx"
#include <algorithm>
#include <cmath>
#include <omp.h>

/**
 * @brief Construit un univers par défaut.
 *
 * Initialise un univers tridimensionnel avec :
 * - un domaine cubique de taille 10 dans chaque direction,
 * - une distance de coupure égale à 2.5,
 * - des paramètres de Lennard-Jones par défaut.
 *
 * La grille de cellules est ensuite initialisée.
 * verdict: univers est propriétaire des particules stockées dans particules.
 * Les cellules ne possèdent pas les particules ,
 * elles ne stockent que des pointeurs vers celles-ci afin d’éviter les copies
 * et de permettre leur reclassification spatiale à chaque pas de temps.
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
    this->G = 0.0;
    this->utiliser_potentiel_mur = false;

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
                 double sigma,
                 double G,
                 bool utiliser_potentiel_mur)
{
    this->num_particules = 0;
    this->r_cut = r_cut;
    this->eps = eps;
    this->sigma = sigma;
    this->G = G;
    this->utiliser_potentiel_mur = utiliser_potentiel_mur;

    this->condl_xmin = ConditionLimite::Reflexive;
    this->condl_xmax = ConditionLimite::Reflexive;
    this->condl_ymin = ConditionLimite::Reflexive;
    this->condl_ymax = ConditionLimite::Reflexive;
    this->condl_zmin = ConditionLimite::Reflexive;
    this->condl_zmax = ConditionLimite::Reflexive;

    if (dim < 1 || dim > 3) {
        std::cerr << "Dimension must be between 1 and 3. Setting to 3." << std::endl;
        this->dim = 3;
    } else {
        this->dim = dim;
    }

    if (Lds.size() != static_cast<size_t>(this->dim)) {
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

univers::univers(std::vector<particule*>& v,
                 std::vector<double> Lds,
                 double r_cut,
                 int dim,
                 double eps,
                 double sigma)
    : univers(v, Lds, r_cut, dim, eps, sigma, 0.0, false)
{
}

univers::univers(std::vector<particule*>& v,
                 std::vector<double> Lds,
                 double r_cut,
                 int dim,
                 double eps,
                 double sigma,
                 double G)
    : univers(v, Lds, r_cut, dim, eps, sigma, G, false)
{
}

void univers::clear_particules() {
    for (particule* p : particules) {
        delete p;
    }

    particules.clear();
    num_particules = 0;
}


void univers::copier_parametres_depuis(const univers& other) {
    num_particules = other.num_particules;
    r_cut = other.r_cut;
    G = other.G;
    dim = other.dim;
    Lds = other.Lds;
    ncd = other.ncd;
    eps = other.eps;
    sigma = other.sigma;

    condl_xmin = other.condl_xmin;
    condl_xmax = other.condl_xmax;
    condl_ymin = other.condl_ymin;
    condl_ymax = other.condl_ymax;
    condl_zmin = other.condl_zmin;
    condl_zmax = other.condl_zmax;

    utiliser_potentiel_mur = other.utiliser_potentiel_mur;
}
/**
 * @brief Détruit l'univers et libère la mémoire des particules.
 */
univers::~univers() {
    clear_particules();
}
/**
 * @brief Constructeur de copie.
 * @param other L'objet à copier.
 */
univers::univers(const univers& other) {
    copier_parametres_depuis(other);

    particules.clear();
    particules.reserve(other.particules.size());

    for (const particule* p : other.particules) {
        particules.push_back(new particule(*p));
    }

    num_particules = static_cast<int>(particules.size());

    initialise_cellules();
    place_particules_dans_cellules();
}

/** @brief Active ou désactive l'utilisation de la liste de Verlet.
 * @param actif true pour activer, false pour désactiver.
 * @param skin Valeur du skin pour la liste de Verlet.
 */
void univers::setUtiliserListeVerlet(bool actif, double skin) {
    utiliser_liste_verlet = actif;
    verlet_skin = skin;
    verlet_valide = false;
}

/**
 * @brief Opérateur d'affectation.
 * @param other L'objet à affecter.
 * @return Référence à l'objet affecté.
 */
univers& univers::operator=(const univers& other) {
    if (this == &other) {
        return *this;
    }

    clear_particules();
    cellules.clear();

    copier_parametres_depuis(other);

    particules.reserve(other.particules.size());

    for (const particule* p : other.particules) {
        particules.push_back(new particule(*p));
    }

    num_particules = static_cast<int>(particules.size());

    initialise_cellules();
    place_particules_dans_cellules();

    return *this;
}

/**
 * @brief Constructeur de déplacement.
 * @param other L'objet à déplacer.
 */
univers::univers(univers&& other) noexcept
    : particules(std::move(other.particules)),
      num_particules(other.num_particules),
      r_cut(other.r_cut),
      G(other.G),
      dim(other.dim),
      Lds(std::move(other.Lds)),
      ncd(std::move(other.ncd)),
      cellules(std::move(other.cellules)),
      eps(other.eps),
      sigma(other.sigma),
      condl_xmin(other.condl_xmin),
      condl_xmax(other.condl_xmax),
      condl_ymin(other.condl_ymin),
      condl_ymax(other.condl_ymax),
      condl_zmin(other.condl_zmin),
      condl_zmax(other.condl_zmax),
      utiliser_potentiel_mur(other.utiliser_potentiel_mur)
{
    other.num_particules = 0;
}


/** @brief Opérateur d'affectation par déplacement.
 * @param other L'objet à déplacer.
 * @return Référence à l'objet déplacé.
 */
univers& univers::operator=(univers&& other) noexcept {
    if (this == &other) {
        return *this;
    }

    clear_particules();
    cellules.clear();

    particules = std::move(other.particules);
    num_particules = other.num_particules;

    r_cut = other.r_cut;
    G = other.G;
    dim = other.dim;

    Lds = std::move(other.Lds);
    ncd = std::move(other.ncd);
    cellules = std::move(other.cellules);

    eps = other.eps;
    sigma = other.sigma;

    condl_xmin = other.condl_xmin;
    condl_xmax = other.condl_xmax;
    condl_ymin = other.condl_ymin;
    condl_ymax = other.condl_ymax;
    condl_zmin = other.condl_zmin;
    condl_zmax = other.condl_zmax;

    utiliser_potentiel_mur = other.utiliser_potentiel_mur;

    other.num_particules = 0;

    return *this;
}

void univers::reserveParticules(size_t n) {
    particules.reserve(n);
}


/** @brief Retourne l'ensemble des particules.
 * @return Référence constante vers le vecteur de particules.
 */
const std::vector<particule*>& univers::getParticules() const {
    return particules;
}
/** @brief Active ou désactive l'utilisation du potentiel de mur.
 * @param actif true pour activer, false pour désactiver.
 */
void univers::setUtiliserPotentielMur(bool actif) {
    utiliser_potentiel_mur = actif;
}


/** @brief Retourne le nombre de particules.
 * @return Nombre de particules.
 */
int univers::getNumParticules() const {
    return num_particules;
}


/** @brief Retourne la distance de coupure.
 * @return Valeur de r_cut.
 */
double univers::getRCut() const {
    return r_cut;
}


/** @brief Retourne la dimension de l'univers.
 * @return Dimension.
 */
int univers::getDim() const {
    return dim;
}


/** @brief Retourne les tailles du domaine.
 * @return Référence constante vers le vecteur des tailles.
 */
const std::vector<double>& univers::getLds() const {
    return Lds;
}

/** @brief Retourne le nombre de cellules par direction.
 * @return Référence constante vers le vecteur ncd.
 */
const std::vector<int>& univers::getNcd() const {
    return ncd;
}


/** @brief Retourne la grille des cellules.
 * @return Référence constante vers le vecteur des cellules.
 */
const std::vector<cellule>& univers::getCellules() const {
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
    p->setIndexUnivers(static_cast<int>(this->particules.size()) - 1);

    this->num_particules++;

    cellule* c = place_particule_dans_cellule(p);
    if (c->getParticules().size() == 1)
        cellules_occupees.push_back(c);

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
    //place_particules_dans_cellules();

    //calcule_forces();
    if (!utiliser_liste_verlet) {
        place_particules_dans_cellules();
    }

    calcule_forces();

    if (utiliser_potentiel_mur) {
        applique_potentiel_mur();
    }

    applique_gravite();

    for (particule* p : particules) {
        p->avance_vitesse_verlet(dt);
    }

    if (!aucune_cond_limite) {
        const bool limites_ont_modifie = applique_conditions_limites();

        if (limites_ont_modifie) {
            verlet_valide = false;
        }
    }
}


bool univers::doit_reconstruire_liste_verlet() const {
    if (!verlet_valide) return true;

    const double seuil2 = 0.25 * verlet_skin * verlet_skin;

    for (int i = 0; i < static_cast<int>(particules.size()); ++i) {
        const vecteur& p = particules[i]->getPosition();

        const double dx = p.getX() - verlet_x0[i];
        const double dy = p.getY() - verlet_y0[i];
        const double dz = p.getZ() - verlet_z0[i];

        const double d2 = dx * dx + dy * dy + dz * dz;

        if (d2 > seuil2) {
            return true;
        }
    }

    return false;
}

void univers::construire_liste_verlet() {
    paires_verlet.clear();

    const int N = static_cast<int>(particules.size());

    verlet_x0.resize(N);
    verlet_y0.resize(N);
    verlet_z0.resize(N);

    const double r_list = r_cut + verlet_skin;
    const double r_list2 = r_list * r_list;

    const int n_layers = static_cast<int>(std::ceil(r_list / r_cut));

    //paires_verlet.reserve(static_cast<size_t>(N) * 40);
    const size_t reserve_visee =
        static_cast<size_t>(N) * ((dim == 2) ? 40 : 100);

    if (paires_verlet.capacity() < reserve_visee) {
        paires_verlet.reserve(reserve_visee);
    }
    
    if (dim == 2) {
        const int nx = ncd[0];
        const int ny = ncd[1];

        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                const int idx = indice_cellule(ix, iy, 0);
                cellule& c = cellules[idx];

                const auto& parts = c.getParticules();
                if (parts.empty()) continue;

                for (int dx = -n_layers; dx <= n_layers; ++dx) {
                    for (int dy = -n_layers; dy <= n_layers; ++dy) {
                        const int jx = ix + dx;
                        const int jy = iy + dy;

                        if (jx < 0 || jx >= nx || jy < 0 || jy >= ny) continue;

                        const int idx_voisin = indice_cellule(jx, jy, 0);
                        if (idx_voisin < idx) continue;

                        cellule& cv = cellules[idx_voisin];
                        const auto& vois = cv.getParticules();
                        if (vois.empty()) continue;

                        if (idx_voisin == idx) {
                            for (size_t a = 0; a < parts.size(); ++a) {
                                particule* pi = parts[a];
                                const int idi = pi->getIndexUnivers();

                                const vecteur& posi = pi->getPosition();
                                const double xi = posi.getX();
                                const double yi = posi.getY();

                                for (size_t b = a + 1; b < parts.size(); ++b) {
                                    particule* pj = parts[b];
                                    const int idj = pj->getIndexUnivers();

                                    const vecteur& posj = pj->getPosition();

                                    const double rx = posj.getX() - xi;
                                    const double ry = posj.getY() - yi;
                                    const double dist2 = rx * rx + ry * ry;

                                    if (dist2 <= r_list2) {
                                        paires_verlet.emplace_back(idi, idj);
                                    }
                                }
                            }
                        } else {
                            for (particule* pi : parts) {
                                const int idi = pi->getIndexUnivers();

                                const vecteur& posi = pi->getPosition();
                                const double xi = posi.getX();
                                const double yi = posi.getY();

                                for (particule* pj : vois) {
                                    const int idj = pj->getIndexUnivers();

                                    const vecteur& posj = pj->getPosition();

                                    const double rx = posj.getX() - xi;
                                    const double ry = posj.getY() - yi;
                                    const double dist2 = rx * rx + ry * ry;

                                    if (dist2 <= r_list2) {
                                        paires_verlet.emplace_back(idi, idj);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    else if (dim == 3) {
        const int nx = ncd[0];
        const int ny = ncd[1];
        const int nz = ncd[2];

        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                for (int iz = 0; iz < nz; ++iz) {
                    const int idx = indice_cellule(ix, iy, iz);
                    cellule& c = cellules[idx];

                    const auto& parts = c.getParticules();
                    if (parts.empty()) continue;

                    for (int dx = -n_layers; dx <= n_layers; ++dx) {
                        for (int dy = -n_layers; dy <= n_layers; ++dy) {
                            for (int dz = -n_layers; dz <= n_layers; ++dz) {
                                const int jx = ix + dx;
                                const int jy = iy + dy;
                                const int jz = iz + dz;

                                if (jx < 0 || jx >= nx ||
                                    jy < 0 || jy >= ny ||
                                    jz < 0 || jz >= nz) {
                                    continue;
                                }

                                const int idx_voisin = indice_cellule(jx, jy, jz);
                                if (idx_voisin < idx) continue;

                                cellule& cv = cellules[idx_voisin];
                                const auto& vois = cv.getParticules();
                                if (vois.empty()) continue;

                                if (idx_voisin == idx) {
                                    for (size_t a = 0; a < parts.size(); ++a) {
                                        particule* pi = parts[a];
                                        const int idi = pi->getIndexUnivers();

                                        const vecteur& posi = pi->getPosition();
                                        const double xi = posi.getX();
                                        const double yi = posi.getY();
                                        const double zi = posi.getZ();

                                        for (size_t b = a + 1; b < parts.size(); ++b) {
                                            particule* pj = parts[b];
                                            const int idj = pj->getIndexUnivers();

                                            const vecteur& posj = pj->getPosition();

                                            const double rx = posj.getX() - xi;
                                            const double ry = posj.getY() - yi;
                                            const double rz = posj.getZ() - zi;
                                            const double dist2 =
                                                rx * rx + ry * ry + rz * rz;

                                            if (dist2 <= r_list2) {
                                                paires_verlet.emplace_back(idi, idj);
                                            }
                                        }
                                    }
                                } else {
                                    for (particule* pi : parts) {
                                        const int idi = pi->getIndexUnivers();

                                        const vecteur& posi = pi->getPosition();
                                        const double xi = posi.getX();
                                        const double yi = posi.getY();
                                        const double zi = posi.getZ();

                                        for (particule* pj : vois) {
                                            const int idj = pj->getIndexUnivers();

                                            const vecteur& posj = pj->getPosition();

                                            const double rx = posj.getX() - xi;
                                            const double ry = posj.getY() - yi;
                                            const double rz = posj.getZ() - zi;
                                            const double dist2 =
                                                rx * rx + ry * ry + rz * rz;

                                            if (dist2 <= r_list2) {
                                                paires_verlet.emplace_back(idi, idj);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i < N; ++i) {
        const vecteur& p = particules[i]->getPosition();

        verlet_x0[i] = p.getX();
        verlet_y0[i] = p.getY();
        verlet_z0[i] = p.getZ();
    }

    verlet_valide = true;
}

/** @brief Calcule l'énergie cinétique de l'univers.
 * @return Valeur de l'énergie cinétique.
 */
double univers::energie_cinetique() const {
    double Ec = 0.0;

    for (const particule* p : particules) {
        const vecteur& v = p->getVitesse();
        double v2 = v.getX()*v.getX() + v.getY()*v.getY() + v.getZ()*v.getZ();
        Ec += 0.5 * p->getMasse() * v2;
    }

    return Ec;
}


/** @brief Calcule l'énergie potentielle de l'univers.
 * @return Valeur de l'énergie potentielle.
 */
double univers::energie_potentielle() const {
    double Ep = 0.0;

    // Énergie potentielle de gravité
    for (const particule* p : particules) {

        double y = p->getPosition().getY();
        double m = p->getMasse();
        Ep += m * G * y;
    }

    // Énergie potentielle Lennard-Jones avec la grille de cellules
    const double r_cut2 = r_cut * r_cut;
    const double sigma2 = sigma * sigma;

    for (const cellule* c : cellules_occupees) {
        const auto& parts = c->getParticules();
        for (const cellule* v : c->getVoisins()) {
            if (v == c) {
                // Paires internes à la même cellule
                for (size_t i = 0; i < parts.size(); ++i) {
                    for (size_t j = i + 1; j < parts.size(); ++j) {
                        const particule* pi = parts[i];
                        const particule* pj = parts[j];


                        /* vecteur rij = pj->getPosition() - pi->getPosition();
                        double dist2 = rij.norme2() + 1e-12; */

                        // optimisation: passage par composantes
                        double rx = pj->getPosition().getX() - pi->getPosition().getX();
                        double ry = pj->getPosition().getY() - pi->getPosition().getY();
                        double rz = pj->getPosition().getZ() - pi->getPosition().getZ();
                        double dist2 = rx*rx + ry*ry + rz*rz + 1e-12;

                        if (dist2 > r_cut2) {
                            continue;
                        }

                        double sr2 = sigma2 / dist2;
                        double sr6 = sr2 * sr2 * sr2;

                        Ep += 4.0 * eps * (sr6 * sr6 - sr6);
                    }
                }
            } else {
                // Paires entre deux cellules voisines distinctes
                const auto& vois = v->getParticules();

                for (const particule* pi : parts) {
                    for (const particule* pj : vois) {

                        double rx = pj->getPosition().getX() - pi->getPosition().getX();
                        double ry = pj->getPosition().getY() - pi->getPosition().getY();
                        double rz = pj->getPosition().getZ() - pi->getPosition().getZ();
                        double dist2 = rx*rx + ry*ry + rz*rz + 1e-12;

                        if (dist2 > r_cut2) {
                            continue;
                        }

                        double sr2 = sigma2 / dist2;
                        double sr6 = sr2 * sr2 * sr2;

                        Ep += 4.0 * eps * (sr6 * sr6 - sr6);
                    }
                }
            }
        }
    }

    return Ep;
}


/** @brief Calcule l'énergie mécanique totale de l'univers.
 * @return Valeur de l'énergie mécanique.
 */
double univers::energie_mecanique() const {
    return energie_cinetique() + energie_potentielle();
}   


/** @brief Calcule la force exercée par le mur sur une particule.
 * @param r Distance à laquelle la force est calculée.
 * @return Valeur de la force.
 */
double univers::calcule_force_mur(double r) const {
    r = std::max(r, 1e-12);

    const double a = sigma / (2.0 * r);
    const double a6 = a * a * a * a * a * a;

    return -24.0 * eps * (1.0 / (2.0 * r)) * a6 * (1.0 - 2.0 * a6);
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
cellule* univers::place_particule_dans_cellule(particule* p) {
    const vecteur& pos = p->getPosition();

    int ix = static_cast<int>(pos.getX() / r_cut);
    if (ix < 0) ix = 0;
    else if (ix >= ncd[0]) ix = ncd[0] - 1;

    int iy = 0, iz = 0;

    if (dim >= 2) {
        iy = static_cast<int>(pos.getY() / r_cut);
        if (iy < 0) iy = 0;
        else if (iy >= ncd[1]) iy = ncd[1] - 1;
    }

    if (dim == 3) {
        iz = static_cast<int>(pos.getZ() / r_cut);
        if (iz < 0) iz = 0;
        else if (iz >= ncd[2]) iz = ncd[2] - 1;
    }

    cellule* c = &cellules[indice_cellule(ix, iy, iz)];
    c->ajoute_particule(p);
    return c;
}

/**
 * @brief Replace toutes les particules dans la grille de cellules.
 *
 * Toutes les cellules sont d'abord vidées, puis chaque particule
 * est replacée à partir de sa position actuelle.
 */
void univers::place_particules_dans_cellules() {
    for (cellule* c : cellules_occupees) {
        c->vide();
    }

    cellules_occupees.clear();

    for (int i = 0; i < static_cast<int>(particules.size()); ++i) {
        particule* p = particules[i];

        // Sécurité OpenMP : index toujours contigu, même après absorption.
        p->setIndexUnivers(i);

        cellule* c = place_particule_dans_cellule(p);

        if (c->getParticules().size() == 1) {
            cellules_occupees.push_back(c);
        }
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
 * Chaque cellule identifie ses voisines dans un voisinage local (cube de taille 3^dim).
 * Afin d'éviter de traiter deux fois les mêmes interactions entre particules,
 * on ne conserve que les cellules voisines dont l'indice linéaire est supérieur
 * ou égal à celui de la cellule courante.
 *
 * Cette convention impose un parcours orienté du voisinage et garantit que chaque
 * paire de cellules est traitée une seule fois. En exploitant la troisième loi
 * de Newton, la force opposée est appliquée simultanément aux deux particules,
 * éliminant ainsi les calculs redondants.
 *
 * Cette optimisation réduit significativement le coût du calcul des interactions.
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
        c.reserveVoisins(nb_voisins_max);
    }

    if (dim == 1) {
        for (int ix = 0; ix < ncd[0]; ++ix) {
            int idx = indice_cellule(ix);

            for (int dx = -1; dx <= 1; ++dx) {
                int nix = ix + dx;

                if (nix >= 0 && nix < ncd[0]) {
                    int idx_voisin = indice_cellule(nix);

                    if (idx_voisin >= idx) {
                        cellules[idx].ajoute_voisin(&cellules[idx_voisin]);
                    }
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

                            if (idx_voisin >= idx) {
                                cellules[idx].ajoute_voisin(&cellules[idx_voisin]);
                            }
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

                                    if (idx_voisin >= idx) {
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

    this->aucune_cond_limite = (cond_xmin == ConditionLimite::Aucune &&
                 cond_xmax == ConditionLimite::Aucune &&
                 cond_ymin == ConditionLimite::Aucune &&
                 cond_ymax == ConditionLimite::Aucune &&
                 cond_zmin == ConditionLimite::Aucune &&
                 cond_zmax == ConditionLimite::Aucune);
}


/** @brief Applique les conditions aux limites à toutes les particules. */
bool univers::applique_conditions_limites() {
    bool modification = false;

    const bool has_absorbante =
        condl_xmin == ConditionLimite::Absorbante ||
        condl_xmax == ConditionLimite::Absorbante ||
        condl_ymin == ConditionLimite::Absorbante ||
        condl_ymax == ConditionLimite::Absorbante ||
        condl_zmin == ConditionLimite::Absorbante ||
        condl_zmax == ConditionLimite::Absorbante;

    if (!has_absorbante) {
        for (particule* p : particules) {
            const vecteur& avant = p->getPosition();

            const double x0 = avant.getX();
            const double y0 = avant.getY();
            const double z0 = avant.getZ();

            applique_conditions_limites_particule(p);

            const vecteur& apres = p->getPosition();

            if (apres.getX() != x0 || apres.getY() != y0 || apres.getZ() != z0) {
                modification = true;
            }
        }

        return modification;
    }

    std::vector<particule*> survivantes;
    survivantes.reserve(particules.size());

    for (particule* p : particules) {
        const vecteur& avant = p->getPosition();

        const double x0 = avant.getX();
        const double y0 = avant.getY();
        const double z0 = avant.getZ();

        if (applique_conditions_limites_particule(p)) {
            const vecteur& apres = p->getPosition();

            if (apres.getX() != x0 || apres.getY() != y0 || apres.getZ() != z0) {
                modification = true;
            }

            survivantes.push_back(p);
        } else {
            delete p;
            modification = true;
        }
    }

    particules = std::move(survivantes);
    num_particules = static_cast<int>(particules.size());

    return modification;
}

/**
 * @brief Applique les conditions aux limites à une particule.
 *
 * @param p Particule à traiter.
 * @return false si la particule doit être absorbée, true sinon.
 */

bool univers::applique_conditions_limites_particule(particule* p) {
    vecteur pos = p->getPosition();
    vecteur vit = p->getVitesse();

    bool garder = true;

    double Lx = Lds[0];
    double Ly = (dim >= 2) ? Lds[1] : 0.0;
    double Lz = (dim == 3) ? Lds[2] : 0.0;

    // =====================
    // Bord xmin : x = 0
    // =====================
    if (pos.getX() < 0.0) {
        if (condl_xmin == ConditionLimite::Reflexive) {
            if (!utiliser_potentiel_mur) {
                pos.setX(-pos.getX());
                vit.setX(-vit.getX());
            }
        }
        else if (condl_xmin == ConditionLimite::Periodique) {
            pos.setX(pos.getX() + Lx);
        }
        else if (condl_xmin == ConditionLimite::Absorbante) {
            garder = false;
        }
    }

    // =====================
    // Bord xmax : x = Lx
    // =====================
    if (pos.getX() > Lx) {
        if (condl_xmax == ConditionLimite::Reflexive) {
            if (!utiliser_potentiel_mur) {
                pos.setX(2.0 * Lx - pos.getX());
                vit.setX(-vit.getX());
            }
        }
        else if (condl_xmax == ConditionLimite::Periodique) {
            pos.setX(pos.getX() - Lx);
        }
        else if (condl_xmax == ConditionLimite::Absorbante) {
            garder = false;
        }
    }

    // =====================
    // Bord ymin : y = 0
    // =====================
    if (pos.getY() < 0.0) {
        if (condl_ymin == ConditionLimite::Reflexive) {
            if (!utiliser_potentiel_mur) {
                pos.setY(-pos.getY());
                vit.setY(-vit.getY());
            }
        }
        else if (condl_ymin == ConditionLimite::Periodique) {
            pos.setY(pos.getY() + Ly);
        }
        else if (condl_ymin == ConditionLimite::Absorbante) {
            garder = false;
        }
    }

    // =====================
    // Bord ymax : y = Ly
    // =====================
    if (pos.getY() > Ly) {
        if (condl_ymax == ConditionLimite::Reflexive) {
            if (!utiliser_potentiel_mur) {
                pos.setY(2.0 * Ly - pos.getY());
                vit.setY(-vit.getY());
            }
        }
        else if (condl_ymax == ConditionLimite::Periodique) {
            pos.setY(pos.getY() - Ly);
        }
        else if (condl_ymax == ConditionLimite::Absorbante) {
            garder = false;
        }
    }

    // =====================
    // Bord zmin : z = 0
    // seulement si dim == 3
    // =====================
    if (dim == 3 && pos.getZ() < 0.0) {
        if (condl_zmin == ConditionLimite::Reflexive) {
            if (!utiliser_potentiel_mur) {
                pos.setZ(-pos.getZ());
                vit.setZ(-vit.getZ());
            }
        }
        else if (condl_zmin == ConditionLimite::Periodique) {
            pos.setZ(pos.getZ() + Lz);
        }
        else if (condl_zmin == ConditionLimite::Absorbante) {
            garder = false;
        }
    }

    // =====================
    // Bord zmax : z = Lz
    // seulement si dim == 3
    // =====================
    if (dim == 3 && pos.getZ() > Lz) {
        if (condl_zmax == ConditionLimite::Reflexive) {
            if (!utiliser_potentiel_mur) {
                pos.setZ(2.0 * Lz - pos.getZ());
                vit.setZ(-vit.getZ());
            }
        }
        else if (condl_zmax == ConditionLimite::Periodique) {
            pos.setZ(pos.getZ() - Lz);
        }
        else if (condl_zmax == ConditionLimite::Absorbante) {
            garder = false;
        }
    }

    p->setPosition(pos);
    p->setVitesse(vit);

    return garder;
}

void univers::limite_vitesses(int N1, int N2) {
    double Ec_D = 0.005 * (N1 + N2);
    double Ec = energie_cinetique();

    if (Ec <= 1e-12) return;

    double beta = std::sqrt(Ec_D / Ec);

    for (particule* p : particules) {
        vecteur v = p->getVitesse();
        p->setVitesse(v * beta);
    }
}

void univers::calcule_forces() {
    if (utiliser_liste_verlet && (dim == 2 || dim == 3)) {
        if (!verlet_valide || doit_reconstruire_liste_verlet()) {
            place_particules_dans_cellules();
            construire_liste_verlet();
        }

        calcule_forces_verlet_omp();
        return;
    }

    if (dim == 2) {
        calcule_forces_omp_2d();
    } else {
        calcule_forces_omp_3d();
    }
}

void univers::calcule_forces_sequentiel() {
    if (dim == 2) {
        calcule_forces_sequentiel_2d();
    } else {
        calcule_forces_sequentiel_3d();
    }
}

void univers::calcule_forces_verlet_omp() {
    if (dim == 2) {
        calcule_forces_verlet_omp_2d();
    } else {
        calcule_forces_verlet_omp_3d();
    }
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
void univers::calcule_forces_omp_3d() {
    const int N = static_cast<int>(particules.size());
    if (N == 0) return; // pas de particule, pas de force à calculer

    const int nb_threads = omp_get_max_threads(); // nombre de threads disponibles

    // Pour les petits systems ou peu de cellules occupées, le parallélisme peut ne être rentable (reste à vérifier)
    if (nb_threads <= 1 || N < 6000 || cellules_occupees.size() < 4 * static_cast<size_t>(nb_threads)) {
        calcule_forces_sequentiel_3d();
        return;
    }

    prepare_omp_force_buffers(); // alloue et initialise les buffers de forces pour OpenMP

    const double r_cut2 = r_cut * r_cut;
    const double sigma2 = sigma * sigma;
    const double coeff = 24.0 * eps;

    const int nb_cellules_occ = static_cast<int>(cellules_occupees.size()); // pas besoin de consulter toutes les cellules, seulement celles qui sont occupées

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num(); // identifiant du thread courant
        const size_t base = static_cast<size_t>(tid) * N; // chaque thread écrit dans sa propre section du buffer de forces pour éviter 
                                                         // les conflits d'écriture et avec un base de N*tid, on est sur qu'il y a pas 
                                                         // des chauvauchements entre les threads car chacun a N (nobre de particules) cases pour lui dans le buffer


        /* Ici, c'est beaucoup plus intéressant de distribuer les tâches dynamiquement, en fait des cellules peuvent être
         trop saturées alors que d'autres sont presque vides. Ainsi, mettre static divisera les cellules entre les threads
         une fois pour toutes sont tenir en compte leur charge  non uniforme => des thread vont finir avant des aures et 
         vont s'arreter, alors que avec (dynamic) chaque thread redemande un autre paquet une fois il a fini de traité celui d'avant.*/
        /* 32 cellules par paquets est généralement un bon compromis, et ça vaut pas le coup de déterminer la meilleur taille pour chaque setting donné*/
        #pragma omp for schedule(dynamic, 32)
        for (int ci = 0; ci < nb_cellules_occ; ++ci) {
            cellule* c = cellules_occupees[ci];
            const auto& parts = c->getParticules();

            for (const cellule* v : c->getVoisins()) {
                const auto& vois = v->getParticules();
                if (vois.empty()) continue;

                if (v == c) {
                    for (size_t i = 0; i < parts.size(); ++i) {
                        particule* pi = parts[i];
                        //const int idi = pi->getId(); // Il faut que les particules soient contigues 
                        const int idi = pi->getIndexUnivers(); // On utilise l'index univers qui est mis à jour à chaque évolution, ainsi même après absorption, les indices sont contigus et cohérents avec le buffer de forces
                        const vecteur& posi = pi->getPosition();
                        const double xi = posi.getX();
                        const double yi = posi.getY();
                        const double zi = posi.getZ();

                        for (size_t j = i + 1; j < parts.size(); ++j) {
                            particule* pj = parts[j];
                            //const int idj = pj->getId();
                            const int idj = pj->getIndexUnivers();

                            const vecteur& posj = pj->getPosition();

                            const double rx = posj.getX() - xi;
                            const double ry = posj.getY() - yi;
                            const double rz = posj.getZ() - zi;

                            const double dist2 = rx*rx + ry*ry + rz*rz + 1e-12;
                            if (dist2 > r_cut2) continue;

                            const double sr2 = sigma2 / dist2;
                            const double sr6 = sr2 * sr2 * sr2;
                            const double coeff_force = coeff * sr6 * (1.0 - 2.0 * sr6) / dist2;

                            const double fxi = rx * coeff_force;
                            const double fyi = ry * coeff_force;
                            const double fzi = rz * coeff_force;

                            /* chaque buffer est divié en nb_threads sections de N cases (N = nombre de particules), ainsi chaque
                             thread écrit dans sa propre section pour éviter les conflits d'écriture*/
                            omp_fx[base + idi] += fxi;
                            omp_fy[base + idi] += fyi;
                            omp_fz[base + idi] += fzi;

                            omp_fx[base + idj] -= fxi;
                            omp_fy[base + idj] -= fyi;
                            omp_fz[base + idj] -= fzi;
                        }
                    }
                } else {
                    for (particule* pi : parts) {
                        //const int idi = pi->getId();
                        const int idi = pi->getIndexUnivers();
                        const vecteur& posi = pi->getPosition();
                        const double xi = posi.getX();
                        const double yi = posi.getY();
                        const double zi = posi.getZ();

                        for (particule* pj : vois) {
                            //const int idj = pj->getId();
                            const int idj = pj->getIndexUnivers();

                            const vecteur& posj = pj->getPosition();

                            const double rx = posj.getX() - xi;
                            const double ry = posj.getY() - yi;
                            const double rz = posj.getZ() - zi;

                            const double dist2 = rx*rx + ry*ry + rz*rz + 1e-12;
                            if (dist2 > r_cut2) continue;

                            const double sr2 = sigma2 / dist2;
                            const double sr6 = sr2 * sr2 * sr2;
                            const double coeff_force = coeff * sr6 * (1.0 - 2.0 * sr6) / dist2;

                            const double fxi = rx * coeff_force;
                            const double fyi = ry * coeff_force;
                            const double fzi = rz * coeff_force;
                            /* chaque buffer est divié en nb_threads sections de N cases (N = nombre de particules), ainsi chaque
                             thread écrit dans sa propre section pour éviter les conflits d'écriture*/
                            omp_fx[base + idi] += fxi;
                            omp_fy[base + idi] += fyi;
                            omp_fz[base + idi] += fzi;

                            omp_fx[base + idj] -= fxi;
                            omp_fy[base + idj] -= fyi;
                            omp_fz[base + idj] -= fzi;
                        }
                    }
                }
            }
        }
        /*Ici on parcoure le buffer de forces pour sommer les forces de chaque particule en déterminant les threads correspondants */
        #pragma omp for schedule(static)
        for (int i = 0; i < N; ++i) {
            double sx = 0.0;
            double sy = 0.0;
            double sz = 0.0;

            for (int t = 0; t < nb_threads; ++t) {
                const size_t idx = static_cast<size_t>(t) * N + i;
                sx += omp_fx[idx];
                sy += omp_fy[idx];
                sz += omp_fz[idx];
            }

            particules[i]->ajouterForce(sx, sy, sz);
        }
    }
}


void univers::calcule_forces_sequentiel_3d(){
    const double r_cut2 = r_cut * r_cut;
    const double sigma2 = sigma * sigma;
    const double coeff = 24.0 * this->eps;

    for (cellule* c : cellules_occupees) {
        const auto& parts = c->getParticules();

        for (const cellule* v : c->getVoisins()) {
            if (v == c) {
                for (size_t i = 0; i < parts.size(); ++i) {
                    for (size_t j = i + 1; j < parts.size(); ++j) {

                        particule* pi = parts[i];
                        particule* pj = parts[j];

                        const vecteur& posi = pi->getPosition();
                        const vecteur& posj = pj->getPosition();

                        double rx = posj.getX() - posi.getX();
                        double ry = posj.getY() - posi.getY();
                        double rz = posj.getZ() - posi.getZ();
                        double dist2 = rx*rx + ry*ry + rz*rz + 1e-12;

                        if (dist2 > r_cut2) continue;

                        double sr2 = sigma2 / dist2;
                        double sr6 = sr2 * sr2 * sr2;
                        double coeff_force = coeff * sr6 * (1.0 - 2.0 * sr6) / dist2;

                        pi->ajouterForce( rx * coeff_force,  ry * coeff_force,  rz * coeff_force);
                        pj->ajouterForce(-rx * coeff_force, -ry * coeff_force, -rz * coeff_force);
                    }
                }
            }
            else {
                const auto& vois = v->getParticules();
                if (vois.empty()) continue;

                for (particule* pi : parts) {
                    const vecteur& posi = pi->getPosition();

                    for (particule* pj : vois) {
                        const vecteur& posj = pj->getPosition();

                        double rx = posj.getX() - posi.getX();
                        double ry = posj.getY() - posi.getY();
                        double rz = posj.getZ() - posi.getZ();
                        double dist2 = rx*rx + ry*ry + rz*rz + 1e-12;

                        if (dist2 > r_cut2) continue;

                        double sr2 = sigma2 / dist2;
                        double sr6 = sr2 * sr2 * sr2;
                        double coeff_force = coeff * sr6 * (1.0 - 2.0 * sr6) / dist2;

                        pi->ajouterForce( rx * coeff_force,  ry * coeff_force,  rz * coeff_force);
                        pj->ajouterForce(-rx * coeff_force, -ry * coeff_force, -rz * coeff_force);
                    }
                }
            }
        }
    }
}

void univers::calcule_forces_sequentiel_2d() {
    const double r_cut2 = r_cut * r_cut;
    const double sigma2 = sigma * sigma;
    const double coeff = 24.0 * eps;

    for (cellule* c : cellules_occupees) {
        const auto& parts = c->getParticules();

        for (const cellule* v : c->getVoisins()) {
            const auto& vois = v->getParticules();
            if (vois.empty()) continue;

            if (v == c) {
                for (size_t i = 0; i < parts.size(); ++i) {
                    particule* pi = parts[i];

                    const vecteur& posi = pi->getPosition();
                    const double xi = posi.getX();
                    const double yi = posi.getY();

                    for (size_t j = i + 1; j < parts.size(); ++j) {
                        particule* pj = parts[j];

                        const vecteur& posj = pj->getPosition();

                        const double rx = posj.getX() - xi;
                        const double ry = posj.getY() - yi;

                        const double dist2 = rx * rx + ry * ry + 1e-12;
                        if (dist2 > r_cut2) continue;

                        const double sr2 = sigma2 / dist2;
                        const double sr6 = sr2 * sr2 * sr2;
                        const double coeff_force =
                            coeff * sr6 * (1.0 - 2.0 * sr6) / dist2;

                        const double fx = rx * coeff_force;
                        const double fy = ry * coeff_force;

                        pi->ajouterForce( fx,  fy, 0.0);
                        pj->ajouterForce(-fx, -fy, 0.0);
                    }
                }
            } else {
                for (particule* pi : parts) {
                    const vecteur& posi = pi->getPosition();
                    const double xi = posi.getX();
                    const double yi = posi.getY();

                    for (particule* pj : vois) {
                        const vecteur& posj = pj->getPosition();

                        const double rx = posj.getX() - xi;
                        const double ry = posj.getY() - yi;

                        const double dist2 = rx * rx + ry * ry + 1e-12;
                        if (dist2 > r_cut2) continue;

                        const double sr2 = sigma2 / dist2;
                        const double sr6 = sr2 * sr2 * sr2;
                        const double coeff_force =
                            coeff * sr6 * (1.0 - 2.0 * sr6) / dist2;

                        const double fx = rx * coeff_force;
                        const double fy = ry * coeff_force;

                        pi->ajouterForce( fx,  fy, 0.0);
                        pj->ajouterForce(-fx, -fy, 0.0);
                    }
                }
            }
        }
    }
}

void univers::calcule_forces_omp_2d() {
    const int N = static_cast<int>(particules.size());
    if (N == 0) return;

    const int nb_threads = omp_get_max_threads();

    if (nb_threads <= 1 ||
        N < 6000 ||
        cellules_occupees.size() < 4 * static_cast<size_t>(nb_threads)) {
        calcule_forces_sequentiel_2d();
        return;
    }

    prepare_omp_force_buffers();

    const double r_cut2 = r_cut * r_cut;
    const double sigma2 = sigma * sigma;
    const double coeff = 24.0 * eps;

    const int nb_cellules_occ = static_cast<int>(cellules_occupees.size());

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        const size_t base = static_cast<size_t>(tid) * N;

        #pragma omp for schedule(dynamic, 32)
        for (int ci = 0; ci < nb_cellules_occ; ++ci) {
            cellule* c = cellules_occupees[ci];
            const auto& parts = c->getParticules();

            for (const cellule* v : c->getVoisins()) {
                const auto& vois = v->getParticules();
                if (vois.empty()) continue;

                if (v == c) {
                    for (size_t i = 0; i < parts.size(); ++i) {
                        particule* pi = parts[i];
                        const int idi = pi->getIndexUnivers();

                        const vecteur& posi = pi->getPosition();
                        const double xi = posi.getX();
                        const double yi = posi.getY();

                        for (size_t j = i + 1; j < parts.size(); ++j) {
                            particule* pj = parts[j];
                            const int idj = pj->getIndexUnivers();

                            const vecteur& posj = pj->getPosition();

                            const double rx = posj.getX() - xi;
                            const double ry = posj.getY() - yi;

                            const double dist2 = rx * rx + ry * ry + 1e-12;
                            if (dist2 > r_cut2) continue;

                            const double sr2 = sigma2 / dist2;
                            const double sr6 = sr2 * sr2 * sr2;
                            const double coeff_force =
                                coeff * sr6 * (1.0 - 2.0 * sr6) / dist2;

                            const double fx = rx * coeff_force;
                            const double fy = ry * coeff_force;

                            omp_fx[base + idi] += fx;
                            omp_fy[base + idi] += fy;

                            omp_fx[base + idj] -= fx;
                            omp_fy[base + idj] -= fy;
                        }
                    }
                } else {
                    for (particule* pi : parts) {
                        const int idi = pi->getIndexUnivers();

                        const vecteur& posi = pi->getPosition();
                        const double xi = posi.getX();
                        const double yi = posi.getY();

                        for (particule* pj : vois) {
                            const int idj = pj->getIndexUnivers();

                            const vecteur& posj = pj->getPosition();

                            const double rx = posj.getX() - xi;
                            const double ry = posj.getY() - yi;

                            const double dist2 = rx * rx + ry * ry + 1e-12;
                            if (dist2 > r_cut2) continue;

                            const double sr2 = sigma2 / dist2;
                            const double sr6 = sr2 * sr2 * sr2;
                            const double coeff_force =
                                coeff * sr6 * (1.0 - 2.0 * sr6) / dist2;

                            const double fx = rx * coeff_force;
                            const double fy = ry * coeff_force;

                            omp_fx[base + idi] += fx;
                            omp_fy[base + idi] += fy;

                            omp_fx[base + idj] -= fx;
                            omp_fy[base + idj] -= fy;
                        }
                    }
                }
            }
        }

        #pragma omp for schedule(static)
        for (int i = 0; i < N; ++i) {
            double sx = 0.0;
            double sy = 0.0;

            for (int t = 0; t < nb_threads; ++t) {
                const size_t idx = static_cast<size_t>(t) * N + i;
                sx += omp_fx[idx];
                sy += omp_fy[idx];
            }

            particules[i]->ajouterForce(sx, sy, 0.0);
        }
    }
}

void univers::calcule_forces_verlet_omp_2d() {
    const int N = static_cast<int>(particules.size());
    if (N == 0) return;

    const int nb_threads = omp_get_max_threads();

    if (nb_threads <= 1 || paires_verlet.size() < 10000) {
        const double r_cut2 = r_cut * r_cut;
        const double sigma2 = sigma * sigma;
        const double coeff = 24.0 * eps;

        for (const auto& paire : paires_verlet) {
            particule* pi = particules[paire.first];
            particule* pj = particules[paire.second];

            const vecteur& posi = pi->getPosition();
            const vecteur& posj = pj->getPosition();

            const double rx = posj.getX() - posi.getX();
            const double ry = posj.getY() - posi.getY();

            const double dist2 = rx * rx + ry * ry + 1e-12;
            if (dist2 > r_cut2) continue;

            const double sr2 = sigma2 / dist2;
            const double sr6 = sr2 * sr2 * sr2;
            const double coeff_force =
                coeff * sr6 * (1.0 - 2.0 * sr6) / dist2;

            const double fx = rx * coeff_force;
            const double fy = ry * coeff_force;

            pi->ajouterForce( fx,  fy, 0.0);
            pj->ajouterForce(-fx, -fy, 0.0);
        }

        return;
    }

    prepare_omp_force_buffers();

    const double r_cut2 = r_cut * r_cut;
    const double sigma2 = sigma * sigma;
    const double coeff = 24.0 * eps;

    const int M = static_cast<int>(paires_verlet.size());

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        const size_t base = static_cast<size_t>(tid) * N;

        #pragma omp for schedule(static)
        for (int k = 0; k < M; ++k) {
            const int idi = paires_verlet[k].first;
            const int idj = paires_verlet[k].second;

            particule* pi = particules[idi];
            particule* pj = particules[idj];

            const vecteur& posi = pi->getPosition();
            const vecteur& posj = pj->getPosition();

            const double rx = posj.getX() - posi.getX();
            const double ry = posj.getY() - posi.getY();

            const double dist2 = rx * rx + ry * ry + 1e-12;
            if (dist2 > r_cut2) continue;

            const double sr2 = sigma2 / dist2;
            const double sr6 = sr2 * sr2 * sr2;
            const double coeff_force =
                coeff * sr6 * (1.0 - 2.0 * sr6) / dist2;

            const double fx = rx * coeff_force;
            const double fy = ry * coeff_force;

            omp_fx[base + idi] += fx;
            omp_fy[base + idi] += fy;

            omp_fx[base + idj] -= fx;
            omp_fy[base + idj] -= fy;
        }

        #pragma omp for schedule(static)
        for (int i = 0; i < N; ++i) {
            double sx = 0.0;
            double sy = 0.0;

            for (int t = 0; t < nb_threads; ++t) {
                const size_t idx = static_cast<size_t>(t) * N + i;
                sx += omp_fx[idx];
                sy += omp_fy[idx];
            }

            particules[i]->ajouterForce(sx, sy, 0.0);
        }
    }
}

void univers::calcule_forces_verlet_omp_3d() {
    const int N = static_cast<int>(particules.size());
    if (N == 0) return;

    const int nb_threads = omp_get_max_threads();

    if (nb_threads <= 1 || paires_verlet.size() < 10000) {
        const double r_cut2 = r_cut * r_cut;
        const double sigma2 = sigma * sigma;
        const double coeff = 24.0 * eps;

        for (const auto& paire : paires_verlet) {
            particule* pi = particules[paire.first];
            particule* pj = particules[paire.second];

            const vecteur& posi = pi->getPosition();
            const vecteur& posj = pj->getPosition();

            const double rx = posj.getX() - posi.getX();
            const double ry = posj.getY() - posi.getY();
            const double rz = posj.getZ() - posi.getZ();

            const double dist2 = rx * rx + ry * ry + rz * rz + 1e-12;
            if (dist2 > r_cut2) continue;

            const double sr2 = sigma2 / dist2;
            const double sr6 = sr2 * sr2 * sr2;
            const double coeff_force =
                coeff * sr6 * (1.0 - 2.0 * sr6) / dist2;

            const double fx = rx * coeff_force;
            const double fy = ry * coeff_force;
            const double fz = rz * coeff_force;

            pi->ajouterForce( fx,  fy,  fz);
            pj->ajouterForce(-fx, -fy, -fz);
        }

        return;
    }

    prepare_omp_force_buffers();

    const double r_cut2 = r_cut * r_cut;
    const double sigma2 = sigma * sigma;
    const double coeff = 24.0 * eps;

    const int M = static_cast<int>(paires_verlet.size());

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        const size_t base = static_cast<size_t>(tid) * N;

        #pragma omp for schedule(static)
        for (int k = 0; k < M; ++k) {
            const int idi = paires_verlet[k].first;
            const int idj = paires_verlet[k].second;

            particule* pi = particules[idi];
            particule* pj = particules[idj];

            const vecteur& posi = pi->getPosition();
            const vecteur& posj = pj->getPosition();

            const double rx = posj.getX() - posi.getX();
            const double ry = posj.getY() - posi.getY();
            const double rz = posj.getZ() - posi.getZ();

            const double dist2 = rx * rx + ry * ry + rz * rz + 1e-12;
            if (dist2 > r_cut2) continue;

            const double sr2 = sigma2 / dist2;
            const double sr6 = sr2 * sr2 * sr2;
            const double coeff_force =
                coeff * sr6 * (1.0 - 2.0 * sr6) / dist2;

            const double fx = rx * coeff_force;
            const double fy = ry * coeff_force;
            const double fz = rz * coeff_force;

            omp_fx[base + idi] += fx;
            omp_fy[base + idi] += fy;
            omp_fz[base + idi] += fz;

            omp_fx[base + idj] -= fx;
            omp_fy[base + idj] -= fy;
            omp_fz[base + idj] -= fz;
        }

        #pragma omp for schedule(static)
        for (int i = 0; i < N; ++i) {
            double sx = 0.0;
            double sy = 0.0;
            double sz = 0.0;

            for (int t = 0; t < nb_threads; ++t) {
                const size_t idx = static_cast<size_t>(t) * N + i;
                sx += omp_fx[idx];
                sy += omp_fy[idx];
                sz += omp_fz[idx];
            }

            particules[i]->ajouterForce(sx, sy, sz);
        }
    }
}

void univers::prepare_omp_force_buffers() {
    const int nb_threads = omp_get_max_threads();
    const int N = static_cast<int>(particules.size());

    const bool need_resize =
        nb_threads != omp_threads_alloc ||
        N != omp_particles_alloc;

    if (need_resize) {
        omp_threads_alloc = nb_threads;
        omp_particles_alloc = N;

        const size_t total = static_cast<size_t>(nb_threads) * N;

        omp_fx.resize(total);
        omp_fy.resize(total);

        if (dim == 3) {
            omp_fz.resize(total);
        }
    }

    const size_t total = static_cast<size_t>(omp_threads_alloc) * omp_particles_alloc;

    if (dim == 2) {
        #pragma omp parallel for schedule(static)
        for (long long i = 0; i < static_cast<long long>(total); ++i) {
            omp_fx[i] = 0.0;
            omp_fy[i] = 0.0;
        }
    } else {
        #pragma omp parallel for schedule(static)
        for (long long i = 0; i < static_cast<long long>(total); ++i) {
            omp_fx[i] = 0.0;
            omp_fy[i] = 0.0;
            omp_fz[i] = 0.0;
        }
    }
}

/** @brief Applique le potentiel de mur à toutes les particules.
 * @param u Univers auquel appliquer le potentiel de mur.
 */
void univers::applique_potentiel_mur() {
    const double r_cut_mur = 0.5 * std::pow(2.0, 1.0 / 6.0) * sigma;
    const double eps_pos   = 1e-6;

    // Déterminer quelles faces ont un potentiel de mur actif dès le début pour éviter de faire des tests à l'intérieur de la boucle.
    const bool wall_ymin = (dim >= 2 && condl_ymin == ConditionLimite::Reflexive);
    const bool wall_ymax = (dim >= 2 && condl_ymax == ConditionLimite::Reflexive);
    const bool wall_xmin = (condl_xmin == ConditionLimite::Reflexive);
    const bool wall_xmax = (condl_xmax == ConditionLimite::Reflexive);
    const bool wall_zmin = (dim == 3 && condl_zmin == ConditionLimite::Reflexive);
    const bool wall_zmax = (dim == 3 && condl_zmax == ConditionLimite::Reflexive);

    for (particule* p : particules) {
        const vecteur& pos = p->getPosition();

        if (wall_ymin) {
            const double r = std::max(pos.getY(), eps_pos);
            if (r < r_cut_mur)
                p->ajouterForce(0.0, +calcule_force_mur(r), 0.0);
        }

        if (wall_ymax) {
            const double r = std::max(Lds[1] - pos.getY(), eps_pos);
            if (r < r_cut_mur)
                p->ajouterForce(0.0, -calcule_force_mur(r), 0.0);
        }

        if (wall_xmin) {
            const double r = std::max(pos.getX(), eps_pos);
            if (r < r_cut_mur)
                p->ajouterForce(+calcule_force_mur(r), 0.0, 0.0);
        }

        if (wall_xmax) {
            const double r = std::max(Lds[0] - pos.getX(), eps_pos);
            if (r < r_cut_mur)
                p->ajouterForce(-calcule_force_mur(r), 0.0, 0.0);
        }

        if (wall_zmin) {
            const double r = std::max(pos.getZ(), eps_pos);
            if (r < r_cut_mur)
                p->ajouterForce(0.0, 0.0, +calcule_force_mur(r));
        }

        if (wall_zmax) {
            const double r = std::max(Lds[2] - pos.getZ(), eps_pos);
            if (r < r_cut_mur)
                p->ajouterForce(0.0, 0.0, -calcule_force_mur(r));
        }
    }
}


/** @brief Applique la gravité à toutes les particules.
 * @param u Univers auquel appliquer la gravité.
 */
void univers::applique_gravite() {
    // Pas de direction verticale en 1D dans ce modèle.
    if (G == 0.0 || dim == 1) return;

    for (particule* p : particules) {
        if (dim == 2) {
            p->ajouterForce(0.0, p->getMasse() * G, 0.0);
        }

        else {
            p->ajouterForce(0.0, 0.0, p->getMasse() * G);
        }
    }
}

/** @brief Convertit une condition limite en chaîne de caractères.
 * @param c La condition limite à convertir.
 * @return La chaîne de caractères correspondante.
 */
static const char* condition_to_string(ConditionLimite c) {
    switch (c) {
        case ConditionLimite::Aucune:
            return "Aucune";
        case ConditionLimite::Reflexive:
            return "Reflexive";
        case ConditionLimite::Absorbante:
            return "Absorbante";
        case ConditionLimite::Periodique:
            return "Periodique";
        default:
            return "Inconnue";
    }
}


/** @brief Affiche les conditions aux limites.
 * @param u Univers dont on veut afficher les conditions limites.
 */
void univers::afficherConditionsLimites() const {
    std::cout << "\n=== Conditions limites actuelles ===\n";
    std::cout << "xmin : " << condition_to_string(condl_xmin) << "\n";
    std::cout << "xmax : " << condition_to_string(condl_xmax) << "\n";
    std::cout << "ymin : " << condition_to_string(condl_ymin) << "\n";
    std::cout << "ymax : " << condition_to_string(condl_ymax) << "\n";
    std::cout << "zmin : " << condition_to_string(condl_zmin) << "\n";
    std::cout << "zmax : " << condition_to_string(condl_zmax) << "\n";
    std::cout << "====================================\n\n";
}


void univers::debug_cellules() const {
    size_t max_occ = 0;
    int idx_max = -1;

    for (size_t i = 0; i < cellules.size(); ++i) {
        size_t occ = cellules[i].getParticules().size();
        if (occ > max_occ) {
            max_occ = occ;
            idx_max = static_cast<int>(i);
        }
    }

    int hors_x = 0;
    int hors_y = 0;

    for (const particule* p : particules) {
        const vecteur& pos = p->getPosition();

        if (pos.getX() < 0.0 || pos.getX() > Lds[0]) {
            hors_x++;
        }

        if (dim >= 2 && (pos.getY() < 0.0 || pos.getY() > Lds[1])) {
            hors_y++;
        }
    }

    std::cout << "[DEBUG CELLULES] max_occ = " << max_occ
              << " dans cellule " << idx_max
              << " | hors_x = " << hors_x
              << " | hors_y = " << hors_y
              << "\n";
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
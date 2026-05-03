#include "univers.hxx"
#include <algorithm>
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

    // 1. update positions
    for (particule* p : particules) {
        p->avance_position_verlet(dt);
    }

    // 2. rebuild cells
    place_particules_dans_cellules();

    // 3. recompute forces
    calcule_forces();

    if (utiliser_potentiel_mur) {
        applique_potentiel_mur();
    }

    applique_gravite();

    // 4. update velocities
    for (particule* p : particules) {
        p->avance_vitesse_verlet(dt);
    }

    // CONDITIONS LIMITES EN DERNIER
    applique_conditions_limites();
}


/** @brief Calcule l'énergie cinétique de l'univers.
 * @return Valeur de l'énergie cinétique.
 */
double univers::energie_cinetique() const {
    double Ec = 0.0;

    for (const particule* p : particules) {
        if (p->getMasse() < 2.0) {
            continue;
        }
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

    for (const particule* p : particules) {
        if (p->getMasse() < 2.0){
            continue;
        }
        double y = p->getPosition().getY();
        double m = p->getMasse();
        Ep += m * G * y;
    }

    double r_cut2 = r_cut * r_cut;

    for (size_t i = 0; i < particules.size(); ++i) {
        for (size_t j = i + 1; j < particules.size(); ++j) {
            if (particules[i]->getMasse() < 2.0 || particules[j]->getMasse() < 2.0){
                continue;
            }
            const vecteur& posi = particules[i]->getPosition();
            const vecteur& posj = particules[j]->getPosition();
            vecteur rij = posj - posi;
            double dist2 = rij.norme2();

            if (dist2 > r_cut2){
                continue;
            }
            double r2 = dist2 + 1e-12;
            double sr2 = (sigma * sigma) / r2;
            double sr6 = sr2 * sr2 * sr2;

            Ep += 4.0 * eps * (sr6 * sr6 - sr6);
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

    const double s  = sigma / (2.0 * r);
    const double s2 = s * s;
    const double s6 = s2 * s2 * s2;

    // Force associated with:
    // just taking of the 1/2 so i can check something
    return -24.0 * eps * (1.0 / 2*r) * s6 * (1.0 - 2.0 * s6);
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
        c.getVoisins().reserve(nb_voisins_max);
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
}


/** @brief Applique les conditions aux limites à toutes les particules. */
void univers::applique_conditions_limites() {
    std::vector<particule*> survivantes;
    survivantes.reserve(particules.size());

    for (particule* p : particules) {
        if (applique_conditions_limites_particule(p)) {
            survivantes.push_back(p);
        } else {
            delete p;
        }
    }

    particules = survivantes;
    num_particules = static_cast<int>(particules.size());
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
        p->setForce(vecteur{});
    }

    const double r_cut2 = r_cut * r_cut;
    const double sigma2 = sigma * sigma;
    const double coeff = 24.0 * this->eps;

    // Parcours des cellules
    for (const cellule& c : this->cellules) {
        const auto& parts = c.getParticules();

        for (const cellule* v : c.getVoisins()){
            if (v == &c){
                // Interactions internes à une même cellule
                for (size_t i = 0; i < parts.size(); ++i) {
                    for (size_t j = i + 1; j < parts.size(); ++j) {

                        particule* pi = parts[i];
                        particule* pj = parts[j];

                        const vecteur& posi = pi->getPosition();
                        const vecteur& posj = pj->getPosition();

                        vecteur rij = posj - posi;
                        double dist2 = rij.norme2() + 1e-12;

                        if (dist2 > r_cut2) continue;

                        double sr2 = sigma2 / dist2;
                        double sr6 = sr2 * sr2 * sr2;

                        double coeff_force = coeff * sr6 * (1.0 - 2.0 * sr6) / dist2;

                        pi->ajouterForce(rij.getX() * coeff_force,
                                        rij.getY() * coeff_force,
                                        rij.getZ() * coeff_force);

                        pj->ajouterForce(-rij.getX() * coeff_force,
                                        -rij.getY() * coeff_force,
                                        -rij.getZ() * coeff_force);
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

                        vecteur rij = posj - posi;
                        double dist2 = rij.norme2() + 1e-12;

                        if (dist2 > r_cut2) continue;

                        double sr2 = sigma2 / dist2;
                        double sr6 = sr2 * sr2 * sr2;

                        double coeff_force = coeff * sr6 * (1.0 - 2.0 * sr6) / dist2;

                        pi->ajouterForce(rij.getX() * coeff_force,
                                        rij.getY() * coeff_force,
                                        rij.getZ() * coeff_force);
                        pj->ajouterForce(-rij.getX() * coeff_force,
                                        -rij.getY() * coeff_force,
                                        -rij.getZ() * coeff_force);
                    }
                }
            }
        }
    }
}


/** @brief Applique le potentiel de mur à toutes les particules.
 * @param u Univers auquel appliquer le potentiel de mur.
 */
void univers::applique_potentiel_mur() {
    const double r_cut_mur = 0.5 * std::pow(2.0, 1.0 / 6.0) * sigma;
    const double eps_pos   = 1e-6; // numerical guard for force evaluation only

    for (particule* p : particules) {
        const vecteur& pos = p->getPosition();

        // Bottom wall: y = 0
        if (dim >= 2 && condl_ymin == ConditionLimite::Reflexive) {
            const double r = std::max(pos.getY(), eps_pos);
            if (r < r_cut_mur) {
                p->ajouterForce(0.0, +calcule_force_mur(r), 0.0);
            }
        }

        // Top wall: y = Ly
        if (dim >= 2 && condl_ymax == ConditionLimite::Reflexive) {
            const double r = std::max(Lds[1] - pos.getY(), eps_pos);
            if (r < r_cut_mur) {
                p->ajouterForce(0.0, -calcule_force_mur(r), 0.0);
            }
        }

        // Left wall: x = 0
        if (condl_xmin == ConditionLimite::Reflexive) {
            const double r = std::max(pos.getX(), eps_pos);
            if (r < r_cut_mur) {
                p->ajouterForce(+calcule_force_mur(r), 0.0, 0.0);
            }
        }

        // Right wall: x = Lx
        if (condl_xmax == ConditionLimite::Reflexive) {
            const double r = std::max(Lds[0] - pos.getX(), eps_pos);
            if (r < r_cut_mur) {
                p->ajouterForce(-calcule_force_mur(r), 0.0, 0.0);
            }
        }

        // z walls only in 3D
        if (dim == 3 && condl_zmin == ConditionLimite::Reflexive) {
            const double r = std::max(pos.getZ(), eps_pos);
            if (r < r_cut_mur) {
                p->ajouterForce(0.0, 0.0, +calcule_force_mur(r));
            }
        }

        if (dim == 3 && condl_zmax == ConditionLimite::Reflexive) {
            const double r = std::max(Lds[2] - pos.getZ(), eps_pos);
            if (r < r_cut_mur) {
                p->ajouterForce(0.0, 0.0, -calcule_force_mur(r));
            }
        }
    }
}


/** @brief Applique la gravité à toutes les particules.
 * @param u Univers auquel appliquer la gravité.
 */
void univers::applique_gravite() {
    if (G == 0.0) return;

    for (particule* p : particules) {
        if (dim == 1) {
            // Pas de direction verticale en 1D dans ce modèle.
            continue;
        }

        if (dim == 2) {
            p->ajouterForce(0.0, p->getMasse() * G, 0.0);
        }

        if (dim == 3) {
            p->ajouterForce(0.0, 0.0, p->getMasse() * G);
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
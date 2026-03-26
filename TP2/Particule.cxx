#include <string>
#include <iostream>
#include <set>
#include <list>
#include <vector>
#include <deque>
#include <random>
#include <chrono>
#include <math.h>

using std::string;
using std::list;
using std::set;
using std::vector;
using std::deque;

/**
 * @mainpage Simulation de particules
 *
 * @section intro Introduction
 *
 * Ce projet simule un système de particules en interaction.
 *
 * @section contenu Contenu
 *
 * - Classe vecteur : manipulation des vecteurs
 * - Classe particule : représente une particule
 * - Classe univers : gère les interactions
 *
 * @section algo Algorithme
 *
 * Utilisation de la méthode de Störmer-Verlet pour l'évolution temporelle.
 * On remplit un fichier .txt et on simule grâce à un fichier python.
 */


/**
 * @brief Représente une particule en 2 dimensions.
 *
 * Une particule possède une position, une vitesse, une force,
 * une masse, un identifiant et un type.
 */
class Particule{
  private:
    double x;      ///< Coordonnée x
    double y;      ///< Coordonnée y
    double vx;     ///< Vitesse selon x
    double vy;     ///< Vitesse selon y
    double fx;     ///< Force selon x
    double fy;     ///< Force selon y
    double masse;  ///< Masse de la particule
    int id;        ///< Identifiant unique
    string type;   ///< Type ou nom de la particule

  public:
    /**
     * @brief Constructeur complet d'une particule.
     *
     * @param x position initiale selon x
     * @param y position initiale selon y
     * @param vx vitesse initiale selon x
     * @param vy vitesse initiale selon y
     * @param fx force initiale selon x
     * @param fy force initiale selon y
     * @param masse masse de la particule
     * @param type type de la particule
     * @param id identifiant de la particule
     */
    Particule(double x,double y,double vx,double vy,double fx,
              double fy,double masse,string type,int id){
       this->x = x;
       this->y = y;
       this->vx = vx;
       this->vy = vy;
       this->fx = fx;
       this->fy = fy;
       this->masse = masse;
       this->type = type;
       this->id = id;
    }

    /**
     * @brief Compare deux particules selon leur identifiant.
     *
     * Cet opérateur permet d'insérer les particules dans un `set`.
     *
     * @param other autre particule à comparer
     * @return true si l'identifiant courant est plus petit
     * @return false sinon
     */
    bool operator<(const Particule& other) const {
      return id < other.id;
    }

    /// @brief Retourne la masse.
    double getMasse(){
      return this->masse;
    }

    /// @brief Retourne la coordonnée x.
    double getX(){
      return this->x;
    }

    /// @brief Retourne la coordonnée y.
    double getY(){
      return this->y;
    }

    /// @brief Retourne la vitesse selon x.
    double getVx(){
      return this->vx;
    }

    /// @brief Retourne la vitesse selon y.
    double getVy(){
      return this->vy;
    }

    /// @brief Retourne le type de la particule.
    string getType(){
      return this->type;
    }

    /// @brief Retourne la force selon x.
    double getFx(){
      return this->fx;
    }

    /// @brief Retourne la force selon y.
    double getFy(){
      return this->fy;
    }

    /// @brief Retourne l'identifiant.
    int getId(){
      return this->id;
    }

    /**
     * @brief Calcule la norme de la vitesse.
     *
     * @return vitesse scalaire de la particule
     */
    double getVitesse(){
      return sqrt(this->vx*this->vx + this->vy*this->vy);
    }
};

/**
 * @brief Calcule les forces d'interaction entre les particules.
 *
 * La force de chaque particule est stockée dans le tableau `F`,
 * avec deux composantes : x et y.
 *
 * @param F tableau des forces calculées
 * @param liste_particules ensemble des particules du système
 */
void calcul_force_F(vector<vector<double>>&,vector<Particule>&);

/**
 * @brief Applique l'algorithme de Störmer-Verlet au système.
 *
 * Met à jour les positions et les vitesses des particules,
 * puis calcule les énergies cinétiques et potentielles.
 *
 * @param liste_particules système de particules à faire évoluer
 * @param delta_t pas de temps
 * @param t_end temps final de la simulation
 */
void algo_Stromer(vector<Particule>& liste_particules,double delta_t,double t_end);

/**
 * @brief Programme principal.
 *
 * Compare les temps d'insertion dans plusieurs structures de données
 * (`list`, `set`, `vector`), puis lance une simulation sur un système
 * solaire simplifié.
 *
 * @return 0 en fin d'exécution
 */
int main() {
  std::random_device rd;
  list<Particule> liste_particules;
  set<Particule> set_particules;
  vector<Particule> vector_particules;

  int taille_list = 400000;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  std::uniform_real_distribution<double> vitesse(5.0, 20.0);

  auto start = std::chrono::steady_clock::now();
  for (int i = 0; i < taille_list; ++i){
    double x = dist(mt);
    double y = dist(mt);
    double fx = dist(mt);
    double fy = dist(mt);
    double masse = dist(mt);
    double vx = vitesse(mt);
    double vy = vitesse(mt);
    int id = i;
    Particule p(x,y,vx,vy,fx,fy,masse,"default",id);
    liste_particules.push_back(p);
  }
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "elapsed time for list: " << elapsed_seconds.count() << "s\n";

  auto start_2 = std::chrono::steady_clock::now();
  for (int i = 0; i < taille_list; ++i){
    double x = dist(mt);
    double y = dist(mt);
    double fx = dist(mt);
    double fy = dist(mt);
    double masse = dist(mt);
    double vx = vitesse(mt);
    double vy = vitesse(mt);
    int id = i;
    Particule p(x,y,vx,vy,fx,fy,masse,"default",id);
    set_particules.insert(p);
  }
  auto end_2 = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds_2 = end_2-start_2;
  std::cout << "elapsed time for set: " << elapsed_seconds_2.count() << "s\n";

  auto start_3 = std::chrono::steady_clock::now();
  for (int i = 0; i < taille_list; ++i){
    double x = dist(mt);
    double y = dist(mt);
    double fx = dist(mt);
    double fy = dist(mt);
    double masse = dist(mt);
    double vx = vitesse(mt);
    double vy = vitesse(mt);
    int id = i;
    Particule p(x,y,vx,vy,fx,fy,masse,"default",id);
    vector_particules.push_back(p);
  }
  auto end_3 = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds_3 = end_3-start_3;
  std::cout << "elapsed time for vector: " << elapsed_seconds_3.count() << "s\n";

  vector<Particule> systeme_solaire;
  Particule soleil(0,0,0,0,0,0,1,"soleil",0);
  Particule terre(0,1,-1,0,0,0,3.0e-6,"terre",1);
  Particule jupiter(0,5.36,-0.425,0,0,0,9.55e-4,"jupiter",2);
  Particule haley(34.75,0,0,0.0296,0,0,1e-14,"haley",3);

  systeme_solaire.push_back(soleil);
  systeme_solaire.push_back(terre);
  systeme_solaire.push_back(jupiter);
  systeme_solaire.push_back(haley);

  double delta_t = 0.015;
  double t_end = 468.5;
  algo_Stromer(systeme_solaire,delta_t,t_end);
}

/**
 * @brief Calcule les forces exercées sur chaque particule.
 *
 * Pour chaque particule i, on somme les contributions des autres particules j.
 * La force est stockée dans un vecteur à deux composantes : x et y.
 *
 * @param F tableau de sortie contenant les forces
 * @param liste_particules particules du système
 */
void calcul_force_F(vector<vector<double>>& F,vector<Particule>& liste_particules){
  int len_liste = liste_particules.size();
  F.clear();

  for (int i = 0; i < len_liste; i++){
    double F_i_x = 0;
    double F_i_y = 0;
    vector<double> F_i;
    Particule p_i = liste_particules[i];
    double x_i = p_i.getX();
    double y_i = p_i.getY();

    for (int j = 0; j < len_liste; j++){
      if (j == i){
        continue;
      }

      Particule p_j = liste_particules[j];
      double x_j = p_j.getX();
      double y_j = p_j.getY();
      double r_i_j_x = x_j - x_i;
      double r_i_j_y = y_j - y_i;

      double norme_r_i_j = sqrt((p_i.getX() - p_j.getX())*(p_i.getX() - p_j.getX())
                              + (p_i.getY() - p_j.getY())*(p_i.getY() - p_j.getY()));

      F_i_x += p_i.getMasse()*p_j.getMasse()*r_i_j_x/(norme_r_i_j*norme_r_i_j*norme_r_i_j);
      F_i_y += p_i.getMasse()*p_j.getMasse()*r_i_j_y/(norme_r_i_j*norme_r_i_j*norme_r_i_j);
    }

    F_i.push_back(F_i_x);
    F_i.push_back(F_i_y);
    F.push_back(F_i);
  }
}

/**
 * @brief Exécute la simulation temporelle par la méthode de Störmer-Verlet.
 *
 * À chaque pas de temps :
 * - met à jour les positions
 * - recalcule les forces
 * - met à jour les vitesses
 * - calcule les énergies cinétique et potentielle
 * - affiche les positions des particules
 *
 * @param liste_particules système étudié
 * @param delta_t pas de temps
 * @param t_end temps final de simulation
 */
void algo_Stromer(vector<Particule>& liste_particules,double delta_t,double t_end){
  double t = 0;
  vector<vector<double>> F_old;
  vector<vector<double>> F;
  vector<double> Ec;
  vector<double> Ep;
  int len_liste = liste_particules.size();

  calcul_force_F(F,liste_particules);

  while (t < t_end){
    t += delta_t;

    for (int i = 0; i < len_liste; i++ ){
      Particule& p = liste_particules[i];
      vector<double> F_i = F[i];

      double x = delta_t*(p.getVx()+0.5/p.getMasse()*F_i[0]*delta_t) + p.getX();
      double y = delta_t*(p.getVy()+0.5/p.getMasse()*F_i[1]*delta_t) + p.getY();

      liste_particules[i] = Particule(x,y,p.getVx(),p.getVy(),p.getFx(),p.getFy(),
                                      p.getMasse(),p.getType(),p.getId());
      F_old.push_back(F_i);
    }

    calcul_force_F(F,liste_particules);

    for (int i = 0; i < len_liste; i++ ){
      Particule& p = liste_particules[i];
      vector<double> F_i = F[i];
      vector<double> F_old_i = F_old[i];

      double vx = delta_t*0.5/p.getMasse()*(F_i[0] + F_old_i[0]) + p.getVx();
      double vy = delta_t*0.5/p.getMasse()*(F_i[1] + F_old_i[1]) + p.getVy();

      liste_particules[i] = Particule(p.getX(),p.getY(),vx,vy,p.getFx(),p.getFy(),
                                      p.getMasse(),p.getType(),p.getId());
    }

    for (int i = 0; i < len_liste; i++){
      Particule p_i = liste_particules[i];
      double Ep_p = 0;
      double m_i = p_i.getMasse();

      for (int j = i+1; j < len_liste; j++){
        Particule p_j = liste_particules[j];
        double m_j = p_j.getMasse();

        double norme_r_i_j = sqrt((p_i.getX() - p_j.getX())*(p_i.getX() - p_j.getX())
                                + (p_i.getY() - p_j.getY())*(p_i.getY() - p_j.getY()));

        Ep_p += - m_i*m_j/norme_r_i_j;
      }

      double v = p_i.getVitesse();
      double Ec_p = 1/2*p_i.getMasse()*v*v;

      Ec.push_back(Ec_p);
      Ep.push_back(Ep_p);
    }

    std::cout << "t = " << t << "\n";
    for (int i = 0; i < len_liste; i++){
      Particule p_i = liste_particules[i];
      std::cout << "nom = " << p_i.getType() << "\n";
      std::cout << "x = " << p_i.getX() << "\n";
      std::cout << "y = " << p_i.getY() << "\n";
    }

    F_old.clear();
  }
}
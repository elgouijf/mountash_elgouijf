#ifndef IO_HXX
#define IO_HXX

#include <fstream>
#include <string>
#include "univers.hxx"

/* IO n'est pas une classe mais un ensemble de fonctions utilitaires pour la sauvegarde et la visualisation des données de simulation.
 Ainsi on définit pas un namespace juste pour pouvoir faire io:: dans le cxx, on garde les choses simples*/

/**
 * @brief Sauvegarde une frame dans un fichier texte.
 *
 * Chaque ligne contient :
 * frame_id x y z type
 *
 * @param file Flux de sortie.
 * @param u Univers contenant les particules.
 * @param frame_id Identifiant de la frame.
 */
void sauvegarde_frame_txt(std::ofstream& file, const univers& u, int frame_id);

/**
 * @brief Sauvegarde une frame au format VTK legacy (.vtk).
 *
 * Ce format est compatible avec ParaView et correspond à
 * l'ancien export utilisé dans le projet.
 *
 * @param u Univers contenant les particules.
 * @param frame_id Identifiant de la frame.
 * @param dossier Dossier de sortie.
 */
void sauvegarde_frame_vtk(const univers& u, int frame_id, const std::string& dossier);

/**
 * @brief Sauvegarde une frame au format VTK XML (.vtu).
 *
 * Ce format correspond à l'export demandé dans le TP5.
 * Le fichier contient les positions, vitesses et masses des particules.
 *
 * @param u Univers contenant les particules.
 * @param frame_id Identifiant de la frame.
 * @param dossier Dossier de sortie.
 */
void sauvegarde_frame_vtu(const univers& u, int frame_id, const std::string& dossier);



/**
 * @brief Écrit un fichier JSON pour une série de frames.
 *
 * @param nb_frames Nombre de frames.
 * @param dt Pas de temps.
 * @param save_every Fréquence de sauvegarde.
 * @param dossier Dossier de sortie.
 * @param extension Extension des fichiers.
 */
void ecrire_fichier_series_json(int nb_frames, double dt, int save_every,
                                const std::string& dossier,
                                const std::string& extension);


/**
 * @brief Sauvegarde le cadre de l'univers au format VTK.
 *
 * @param Lds Dimensions de l'univers.
 * @param dossier Dossier de sortie.
 */
void sauvegarde_cadre_vtk(const std::vector<double>& Lds, const std::string& dossier);


/**
 * @brief Écrit l'entête du fichier d'énergie.
 *
 * @param file Flux de sortie.
 */
void ecrire_entete_energie(std::ofstream& file);

/**
 * @brief Écrit les valeurs d'énergie dans le fichier.
 *
 * @param file Flux de sortie.
 * @param frame Identifiant de la frame.
 * @param time Temps écoulé.
 * @param Ec Énergie cinétique.
 * @param Ep Énergie potentielle.
 * @param Em Énergie mécanique totale.
 */
void ecrire_energie(std::ofstream& file,
                    int frame,
                    double time,
                    double Ec,
                    double Ep,
                    double Em);

#endif
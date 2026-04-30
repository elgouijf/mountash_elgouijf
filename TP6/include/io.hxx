#ifndef IO_HXX
#define IO_HXX

#include <fstream>
#include <string>
#include "univers.hxx"

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
void sauvegarde_frame_vtk(const univers& u, int frame_id, const std::string& dossier = "vtk_frames");

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
void sauvegarde_frame_vtu(const univers& u, int frame_id, const std::string& dossier = "vtu_frames");

/**
 * @brief Génère un fichier animation.vtk.series pour une suite de fichiers .vtk ou .vtu.
 *
 * @param nb_frames Nombre de frames.
 * @param dt Pas de temps.
 * @param save_every Intervalle de sauvegarde.
 * @param dossier Dossier contenant les fichiers.
 * @param extension Extension des fichiers de frame ("vtk" ou "vtu").
 */
void ecrire_fichier_series_json(int nb_frames, double dt, int save_every,
                                const std::string& dossier,
                                const std::string& extension);



void sauvegarde_cadre_vtk(const std::vector<double>& Lds, const std::string& dossier);
#endif
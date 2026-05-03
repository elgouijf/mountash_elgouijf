#include "io.hxx"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

/**
 * @brief Sauvegarde une frame au format VTK XML (.vtu).
 *
 * Le fichier décrit un nuage de points non connectés.
 * Les champs enregistrés sont :
 * - Position
 * - Velocity
 * - Masse
 *
 * @param u Univers contenant les particules.
 * @param frame_id Identifiant de la frame.
 * @param dossier Dossier de sortie.
 */

void sauvegarde_frame_txt(std::ofstream& file, const univers& u, int frame_id) {
    for (particule* p : u.getParticules()) {
        const vecteur& pos = p->getPosition();
        file << frame_id << " "
             << pos[0] << " "
             << pos[1] << " "
             << pos[2] << " "
             << p->getType() << "\n";
    }
}


/** @brief Sauvegarde une frame au format VTK.
 * @param u Univers contenant les particules.
 * @param frame_id Identifiant de la frame.
 * @param dossier Dossier de sortie.
 */
void sauvegarde_frame_vtk(const univers& u, int frame_id, const std::string& dossier) {
    std::ostringstream nom;
    nom << dossier << "/frame_"
        << std::setw(6) << std::setfill('0') << frame_id
        << ".vtk";

    std::ofstream file(nom.str());
    if (!file.is_open()) {
        std::cerr << "Impossible d'ouvrir " << nom.str() << "\n";
        return;
    }

    const std::vector<particule*>& particules = u.getParticules();
    int N = static_cast<int>(particules.size());

    file << "# vtk DataFile Version 3.0\n";
    file << "Frame " << frame_id << "\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    file << "POINTS " << N << " float\n";
    for (particule* p : particules) {
        const vecteur& pos = p->getPosition();
        file << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
    }

    file << "\nVERTICES " << N << " " << 2 * N << "\n";
    for (int i = 0; i < N; ++i) {
        file << "1 " << i << "\n";
    }

    file << "\nPOINT_DATA " << N << "\n";
    file << "SCALARS type int 1\n";
    file << "LOOKUP_TABLE default\n";
    for (particule* p : particules) {
        file << p->getType() << "\n";
    }
}


/** @brief Sauvegarde une frame au format VTK XML (.vtu).
 * @param u Univers contenant les particules.
 * @param frame_id Identifiant de la frame.
 * @param dossier Dossier de sortie.
 */
void sauvegarde_frame_vtu(const univers& u, int frame_id, const std::string& dossier) {
    std::ostringstream nom;
    nom << dossier << "/frame_"
        << std::setw(6) << std::setfill('0') << frame_id
        << ".vtu";

    std::ofstream file(nom.str());
    if (!file.is_open()) {
        std::cerr << "Impossible d'ouvrir " << nom.str() << "\n";
        return;
    }

    const std::vector<particule*>& particules = u.getParticules();
    const int N = static_cast<int>(particules.size());

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    // file << "    <Piece NumberOfPoints=\"" << N << "\" NumberOfCells=\"0\">\n";
    file << "    <Piece NumberOfPoints=\"" << N << "\" NumberOfCells=\"" << N << "\">\n";
    file << "      <Points>\n";
    file << "        <DataArray name=\"Position\" type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    file << "          ";
    for (particule* p : particules) {
        const vecteur& pos = p->getPosition();
        file << pos[0] << " " << pos[1] << " " << pos[2] << " ";
    }
    file << "\n";
    file << "        </DataArray>\n";
    file << "      </Points>\n";

    file << "      <PointData Vectors=\"vector\">\n";

    file << "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    file << "          ";
    for (particule* p : particules) {
        const vecteur& vit = p->getVitesse();
        file << vit[0] << " " << vit[1] << " " << vit[2] << " ";
    }
    file << "\n";
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"Float32\" Name=\"Masse\" format=\"ascii\">\n";
    file << "          ";
    for (particule* p : particules) {
        file << p->getMasse() << " ";
    }
    file << "\n";
    file << "        </DataArray>\n";

    file << "      </PointData>\n";

    file << "      <Cells>\n";

    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    file << "          ";
    for (int i = 0; i < N; ++i) {
        file << i << " ";
    }
    file << "\n";
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    file << "          ";
    for (int i = 0; i < N; ++i) {
        file << i + 1 << " ";
    }
    file << "\n";
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    file << "          ";
    for (int i = 0; i < N; ++i) {
        file << "1 ";
    }
    file << "\n";
    file << "        </DataArray>\n";

    file << "      </Cells>\n";
    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";
}

/**
 * @brief Génère un fichier animation.vtk.series pour ParaView.
 *
 * Ce fichier permet de charger une séquence temporelle de fichiers .vtu.
 *
 * @param nb_frames Nombre total de frames.
 * @param dt Pas de temps.
 * @param save_every Intervalle de sauvegarde.
 * @param dossier Dossier contenant les .vtu.
 */
void ecrire_fichier_series_json(int nb_frames, double dt, int save_every,
                                const std::string& dossier,
                                const std::string& extension) {
    const std::string nom_series = dossier + "/animation." + extension + ".series";

    std::ofstream file(nom_series);
    if (!file.is_open()) {
        std::cerr << "Impossible d'ouvrir " << nom_series << "\n";
        return;
    }

    file << "{\n";
    file << "  \"file-series-version\" : \"1.0\",\n";
    file << "  \"files\" : [\n";

    for (int i = 0; i < nb_frames; ++i) {
        std::ostringstream nom;
        nom << "frame_" << std::setw(6) << std::setfill('0') << i
            << "." << extension;

        double temps = i * save_every * dt;

        file << "    { \"name\" : \"" << nom.str() << "\", \"time\" : " << temps << " }";
        if (i != nb_frames - 1) file << ",";
        file << "\n";
    }

    file << "  ]\n";
    file << "}\n";
}


/** @brief Écrit l'en-tête du fichier d'énergie.
 * @param file Référence au fichier d'énergie.
 */
void ecrire_entete_energie(std::ofstream& file) {
    file << "frame,time,Ec,Ep,Em\n";
}


/** @brief Écrit les valeurs d'énergie dans le fichier.
 * @param file Référence au fichier d'énergie
* @param frame Identifiant de la frame.
* @param time Temps écoulé.
* @param Ec Énergie cinétique.
* @param Ep Énergie potentielle.
* @param Em Énergie mécanique (Ec + Ep).
*/
void ecrire_energie(std::ofstream& file,
                    int frame,
                    double time,
                    double Ec,
                    double Ep,
                    double Em) {
    file << frame << ","
         << time << ","
         << Ec << ","
         << Ep << ","
         << Em << "\n";
}
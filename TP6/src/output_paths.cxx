#include "output_paths.hxx"
#include <stdexcept>

/** @brief Retourne le répertoire de sortie pour un type donné.
 * @param type Type de fichier de sortie.
 * @return Chemin vers le répertoire de sortie.
 */
std::filesystem::path get_output_dir(OutputType type) {
    const std::filesystem::path root = PROJECT_SOURCE_DIR;

    switch (type) {
        case OutputType::FramesTxt:
            return root / "frames";

        case OutputType::FramesVTK:
            return root / "vtk_frames";

        case OutputType::FramesVTU:
            return root / "vtu_frames";
        
        case OutputType::FramesTxt3D:
            return root / "frames_3d";

        case OutputType::FramesVTK3D:
            return root / "vtk_frames_3d";

        case OutputType::FramesVTU3D:
            return root / "vtu_frames_3d";

        case OutputType::Energy:
            return root / "energy";

        case OutputType::Plots:
            return root / "plots";
    }

    throw std::invalid_argument("Type de sortie inconnu");
}


/** @brief Assure que le répertoire de sortie pour un type donné existe.
 * @param type Type de fichier de sortie.
 * @return Chemin vers le répertoire de sortie.
 */
std::filesystem::path ensure_output_dir(OutputType type) {
    auto dir = get_output_dir(type);
    std::filesystem::create_directories(dir);
    return dir;
}
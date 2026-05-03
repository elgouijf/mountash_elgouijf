#ifndef OUTPUT_PATHS_HXX
#define OUTPUT_PATHS_HXX

#include <filesystem>
/** @brief Enumération des types de fichiers de sortie. */
enum class OutputType {
    FramesTxt,
    FramesVTK,
    FramesVTU,
    FramesTxt3D,
    FramesVTK3D,
    FramesVTU3D,
    Energy,
    Plots
};

/** @brief Retourne le répertoire de sortie pour un type donné.
 * @param type Type de fichier de sortie.
 * @return Chemin vers le répertoire de sortie.
 */
std::filesystem::path get_output_dir(OutputType type);

/** @brief Assure que le répertoire de sortie pour un type donné existe.
 * @param type Type de fichier de sortie.
 * @return Chemin vers le répertoire de sortie.
 */
std::filesystem::path ensure_output_dir(OutputType type);

#endif 
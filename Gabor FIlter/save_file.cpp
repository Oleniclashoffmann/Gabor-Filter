#include "save_file.h"

// Constructor for the save_file class
// It takes in a file path, a Point representing the tip, a duration representing time, and an angle as parameters
save_file::save_file(std::string path, cv::Point2f tip, std::chrono::duration<double> time, double angle)
{
    // Declare two output file streams for writing data to files
    std::ofstream outfile;
    std::ofstream file2;

    // Open the 'Tip' file in append mode (so as not to overwrite existing content)
    outfile.open("C:/Users/ohoff/Documents/Bachelorarbeit/Schreiben/Erro analysis/Angle/Gabor/Fulneedle-Judith/Tip/0grad.txt", std::ios_base::app);

    // Open the 'Angle' file in append mode
    file2.open("C:/Users/ohoff/Documents/Bachelorarbeit/Schreiben/Erro analysis/Angle/Gabor/Fulneedle-Judith/Angle/0grad.txt", std::ios_base::app);

    // Check if the 'Tip' file was opened successfully
    if (!outfile) {
        std::cerr << "Could not open the file!" << std::endl;
    }

    // Write the angle value to the 'Angle' file followed by a newline
    file2 << angle << std::endl;

    // Write the file path, tip coordinates, and the time value to the 'Tip' file in a comma-separated format
    outfile << path << ", " << tip.x << ", " << tip.y << ", " << time.count() << std::endl;

    // Close both files
    outfile.close();
    file2.close();
}

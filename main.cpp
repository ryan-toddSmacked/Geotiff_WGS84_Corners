#include <iostream>
#include <string>
#include <vector>

#include "gtiff_corners.hpp"
#include "gtiff_write.hpp"


int main(int argc, char** argv)
{

    // 1 argument, the name of the file
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <geotiff_file>" << std::endl;
        return 1;
    }

    // Load the geotiff file
    std::string filename_on_disk = argv[1];

    GCP upperLeft, upperRight, lowerLeft, lowerRight;
    double dx, dy;

    GTIFF_CORNERS_ERROR error = gtiff_loadCorners(filename_on_disk, upperLeft, upperRight, lowerLeft, lowerRight, &dx, &dy);
    if (error != GTIFF_CORNERS_SUCCESS)
    {
        std::cerr << "Error: " << error << std::endl;
        return 1;
    }

    fprintf(stdout, " Upper Left  | %.3e , %.3e | --> | %.6e , %.6e |\n", upperLeft.x, upperLeft.y, upperLeft.lon, upperLeft.lat);
    fprintf(stdout, " Upper Right | %.3e , %.3e | --> | %.6e , %.6e |\n", upperRight.x, upperRight.y, upperRight.lon, upperRight.lat);
    fprintf(stdout, " Lower Right | %.3e , %.3e | --> | %.6e , %.6e |\n", lowerRight.x, lowerRight.y, lowerRight.lon, lowerRight.lat);
    fprintf(stdout, " Lower Left  | %.3e , %.3e | --> | %.6e , %.6e |\n", lowerLeft.x, lowerLeft.y, lowerLeft.lon, lowerLeft.lat);

    fprintf(stdout, "dx ~ %lf [m]\ndy ~ %lf [m]\n", dx, dy);


    uint8_t* data = new uint8_t[100*100];
    for (int i = 0; i < 100*100; i++)
    {
        data[i] = i % 256;
    }

    upperLeft.x = 0;
    upperLeft.y = 0;
    upperRight.x = 100;
    upperRight.y = 0;
    lowerLeft.x = 0;
    lowerLeft.y = 100;
    lowerRight.x = 100;
    lowerRight.y = 100;

    GTIFF_WRITE_ERROR write_error = gtiff_write("output.tif", data, 100, 100, upperLeft, upperRight, lowerLeft, lowerRight);

    if (write_error != GTIFF_WRITE_SUCCESS)
    {
        std::cerr << "Error: " << write_error << std::endl;
        return 1;
    }


    return 0;
}


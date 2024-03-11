#pragma once
#ifndef __GTIFF_CORNERS_HPP__
#define __GTIFF_CORNERS_HPP__


#include <string>

struct GCP
{
    double x, y;
    double lon, lat;    // Always in WGS84
};


// Error codes
enum GTIFF_CORNERS_ERROR
{
    GTIFF_CORNERS_SUCCESS = 0,                      // Success
    GTIFF_CORNERS_ERROR_OPEN_FILE,                  // Could not open the file
    GTIFF_CORNERS_ERROR_READ_TAGS,                  // Could not read the geotiff tags
    GTIFF_CORNERS_ERROR_READ_DEFINITION,            // Could not read the geotiff definition
    GTIFF_CORNERS_ERROR_PCS_CONVERT,                // Could not convert to PCS
    GTIFF_CORNERS_ERROR_PROJ4_CONVERT,              // Could not convert WGS84 to lon/lat
    GTIFF_CORNERS_ERROR_UNKOWN_COORDINATE_SYSTEM,   // Unknown coordinate system
    GTIFF_CORNERS_ERROR_UNKNOWN                     // Unknown error
};
/**
 * @brief Load the corners of the geotiff file. Stores the corners in the GCP structs.
 * 
 * @param filename Filename on disk
 * @param upperLeft Upper left corner of the image.
 * @param upperRight Upper right corner of the image.
 * @param lowerLeft Lower left corner of the image.
 * @param lowerRight Lower right corner of the image.
 * @param dx Optional, the estimated pixel size in the x direction. In meters.
 * @param dy Optional, the estimated pixel size in the y direction. In meters.
 * @return (0) on sucess, else see the error codes above.
 */
GTIFF_CORNERS_ERROR gtiff_loadCorners(const std::string &filename, GCP& upperLeft, GCP& upperRight, GCP& lowerLeft, GCP& lowerRight, double* dx=nullptr, double* dy=nullptr);

/**
 * @brief Return a string representation of the error code.
 * 
 * @param error Error code
 * @return const char* NULL terminated string
 */
const char* gtiff_stringError(GTIFF_CORNERS_ERROR error);


#endif // __GTIFF_CORNERS_HPP__

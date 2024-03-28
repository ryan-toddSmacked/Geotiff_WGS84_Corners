#pragma once
#ifndef __GTIFF_WRITE_HPP__
#define __GTIFF_WRITE_HPP__

#include <string>
#include <cstdint>
#include "gtiff_GCP.hpp"



enum GTIFF_WRITE_ERROR
{
    GTIFF_WRITE_SUCCESS = 0,            // Success
    GTIFF_WRITE_ERROR_OPEN_FILE,        // Could not open the file
    GTIFF_WRITE_ERROR_WRITE_TAGS,       // Could not write the geotiff tags
    GTIFF_WRITE_ERROR_WRITE_DEFINITION, // Could not write the geotiff definition
    GTIFF_WRITE_ERROR_UNKNOWN           // Unknown error
};


GTIFF_WRITE_ERROR gtiff_write(const std::string &filename, const uint8_t *data, int32_t width, int32_t height, const GCP& upperLeft, const GCP& upperRight, const GCP& lowerLeft, const GCP& lowerRight);





#endif // __GTIFF_WRITE_HPP__

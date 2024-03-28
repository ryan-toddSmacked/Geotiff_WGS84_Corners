
#include "gtiff_corners.hpp"
#include "haversine.hpp"

#include "geotiff.h"
#include "xtiffio.h"
#include "geo_normalize.h"
#include "geo_simpletags.h"
#include "geovalues.h"
#include "tiffio.h"

#include <cstdint>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static inline int readPixel(GTIF* gtif, GTIFDefn* defn, double x, double y, double& lon, double& lat)
{

    // First conver to PCS
    if (!GTIFImageToPCS(gtif, &x, &y))
    {
        return GTIFF_CORNERS_ERROR_PCS_CONVERT;
    }

    if (defn->Model == ModelTypeGeographic)
    {
        // Already in lon/lat, No-op
    }
    else if(defn->Model == ModelTypeProjected)
    {
        if (!GTIFProj4ToLatLong(defn, 1, &x, &y))
        {
            return GTIFF_CORNERS_ERROR_PROJ4_CONVERT;
        }
    }
    else
    {
        return GTIFF_CORNERS_ERROR_UNKOWN_COORDINATE_SYSTEM;
    }

    lon = x;
    lat = y;

    return GTIFF_CORNERS_SUCCESS;
}



GTIFF_CORNERS_ERROR gtiff_loadCorners(const std::string &filename, GCP& upperLeft, GCP& upperRight, GCP& lowerLeft, GCP& lowerRight, double* dx, double* dy)
{
    TIFF* tif = XTIFFOpen(filename.c_str(), "r");
    if (!tif)
    {
        return GTIFF_CORNERS_ERROR_OPEN_FILE;
    }

    GTIF* gtif = GTIFNew(tif);
    if (!gtif)
    {
        XTIFFClose(tif);
        return GTIFF_CORNERS_ERROR_READ_TAGS;
    }

    GTIFDefn defn;
    if (!GTIFGetDefn(gtif, &defn))
    {
        GTIFFree(gtif);
        XTIFFClose(tif);
        return GTIFF_CORNERS_ERROR_READ_DEFINITION;
    }

    // Read image width and height
    uint32_t width, height;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);

    unsigned short raster_type = RasterPixelIsArea;
    GTIFKeyGetSHORT(gtif, GTRasterTypeGeoKey, &raster_type, 0, 1);

    const double xmin = (raster_type == RasterPixelIsArea) ? 0 : -0.5;
    const double ymin = xmin;
    const double xmax = xmin + width;
    const double ymax = ymin + height;
    double lon, lat;

    int error = GTIFF_CORNERS_SUCCESS;

    error = readPixel(gtif, &defn, xmin, ymin, lon, lat);
    if (error != GTIFF_CORNERS_SUCCESS)
    {
        GTIFFree(gtif);
        XTIFFClose(tif);
        return (GTIFF_CORNERS_ERROR)error;
    }
    upperLeft.lon = lon;
    upperLeft.lat = lat;

    error = readPixel(gtif, &defn, xmax, ymin, lon, lat);
    if (error != GTIFF_CORNERS_SUCCESS)
    {
        GTIFFree(gtif);
        XTIFFClose(tif);
        return (GTIFF_CORNERS_ERROR)error;
    }
    upperRight.lon = lon;
    upperRight.lat = lat;

    error = readPixel(gtif, &defn, xmin, ymax, lon, lat);
    if (error != GTIFF_CORNERS_SUCCESS)
    {
        GTIFFree(gtif);
        XTIFFClose(tif);
        return (GTIFF_CORNERS_ERROR)error;
    }
    lowerLeft.lon = lon;
    lowerLeft.lat = lat;

    error = readPixel(gtif, &defn, xmax, ymax, lon, lat);
    if (error != GTIFF_CORNERS_SUCCESS)
    {
        GTIFFree(gtif);
        XTIFFClose(tif);
        return (GTIFF_CORNERS_ERROR)error;
    }
    lowerRight.lon = lon;
    lowerRight.lat = lat;

    GTIFFree(gtif);
    XTIFFClose(tif);

    upperLeft.x = xmin;
    upperLeft.y = ymin;

    upperRight.x = xmax;
    upperRight.y = ymin;

    lowerLeft.x = xmin;
    lowerLeft.y = ymax;

    lowerRight.x = xmax;
    lowerRight.y = ymax;

    // Calculate the pixel size in meters
    // The pixel size is the distance between the upper left and upper right corners divided by the width of the image
    // To estimate the distance between the corners, we use the haversine formula
    if (dx)
    {
        *dx = haversine(upperLeft.lon, upperLeft.lat, upperRight.lon, upperRight.lat) / (double)width;
    }
    if (dy)
    {
        *dy = haversine(upperLeft.lon, upperLeft.lat, lowerLeft.lon, lowerLeft.lat) / (double)height;
    }

    return GTIFF_CORNERS_SUCCESS;
}


const char* gtiff_stringError(GTIFF_CORNERS_ERROR error)
{
    switch (error)
    {
    case GTIFF_CORNERS_SUCCESS:
        return "GTIFF No error";
        break;
    case GTIFF_CORNERS_ERROR_OPEN_FILE:
        return "GTIFF Could not open the file";
        break;
    case GTIFF_CORNERS_ERROR_READ_TAGS:
        return "GTIFF Could not read the geotiff tags";
        break;
    case GTIFF_CORNERS_ERROR_READ_DEFINITION:
        return "GTIFF Could not read the geotiff definition";
        break;
    case GTIFF_CORNERS_ERROR_PCS_CONVERT:
        return "GTIFF Could not convert to PCS";
        break;
    case GTIFF_CORNERS_ERROR_PROJ4_CONVERT:
        return "GTIFF Could not convert projected coordinate system to WGS84 lon/lat";
        break;
    case GTIFF_CORNERS_ERROR_UNKOWN_COORDINATE_SYSTEM:
        return "GTIFF Unknown coordinate system embedded in the geotiff";
        break;
    case GTIFF_CORNERS_ERROR_UNKNOWN:
        return "GTIFF Unknown error";
        break;
    default:
        return "GTIFF Unknown error";
        break;
    }

    return "GTIFF Unknown error";
}



#include "gtiff_write.hpp"
#include "haversine.hpp"

#include "geotiff.h"
#include "xtiffio.h"
#include "geo_normalize.h"
#include "geo_simpletags.h"
#include "geovalues.h"
#include "tiffio.h"
#include "turku.hpp"

#include <cstdio>


inline void matmul_3x4_4x3(const double* A, const double* B, double* C)
{
    // C = A*B
    // A is 3x4
    // B is 4x3
    // C is 3x3

    // A = [a0, a1, a2, a3]
    //     [a4, a5, a6, a7]
    //     [a8, a9, a10, a11]

    // B = [b0, b1, b2]
    //     [b3, b4, b5]
    //     [b6, b7, b8]
    //     [b9, b10, b11]


    // C[0] = A[0]*B[0] + A[1]*B[3] + A[2]*B[6] + A[3]*B[9]
    C[0] = A[0]*B[0] + A[1]*B[3] + A[2]*B[6] + A[3]*B[9];

    // C[1] = A[0]*B[1] + A[1]*B[4] + A[2]*B[7] + A[3]*B[10]
    C[1] = A[0]*B[1] + A[1]*B[4] + A[2]*B[7] + A[3]*B[10];

    // C[2] = A[0]*B[2] + A[1]*B[5] + A[2]*B[8] + A[3]*B[11]
    C[2] = A[0]*B[2] + A[1]*B[5] + A[2]*B[8] + A[3]*B[11];

    // C[3] = A[4]*B[0] + A[5]*B[3] + A[6]*B[6] + A[7]*B[9]
    C[3] = A[4]*B[0] + A[5]*B[3] + A[6]*B[6] + A[7]*B[9];

    // C[4] = A[4]*B[1] + A[5]*B[4] + A[6]*B[7] + A[7]*B[10]
    C[4] = A[4]*B[1] + A[5]*B[4] + A[6]*B[7] + A[7]*B[10];

    // C[5] = A[4]*B[2] + A[5]*B[5] + A[6]*B[8] + A[7]*B[11]
    C[5] = A[4]*B[2] + A[5]*B[5] + A[6]*B[8] + A[7]*B[11];

    // C[6] = A[8]*B[0] + A[9]*B[3] + A[10]*B[6] + A[11]*B[9]
    C[6] = A[8]*B[0] + A[9]*B[3] + A[10]*B[6] + A[11]*B[9];

    // C[7] = A[8]*B[1] + A[9]*B[4] + A[10]*B[7] + A[11]*B[10]
    C[7] = A[8]*B[1] + A[9]*B[4] + A[10]*B[7] + A[11]*B[10];

    // C[8] = A[8]*B[2] + A[9]*B[5] + A[10]*B[8] + A[11]*B[11]
    C[8] = A[8]*B[2] + A[9]*B[5] + A[10]*B[8] + A[11]*B[11];
}

inline void matmul_3x4_4x1(const double* A, const double* B, double* C)
{
    // C = A*B
    // A is 3x4
    // B is 4x1

    // A = [a0, a1, a2, a3]
    //     [a4, a5, a6, a7]
    //     [a8, a9, a10, a11]

    // B = [b0]
    //     [b1]
    //     [b2]
    //     [b3]

    // C[0] = A[0]*B[0] + A[1]*B[1] + A[2]*B[2] + A[3]*B[3]
    C[0] = A[0]*B[0] + A[1]*B[1] + A[2]*B[2] + A[3]*B[3];

    // C[1] = A[4]*B[0] + A[5]*B[1] + A[6]*B[2] + A[7]*B[3]
    C[1] = A[4]*B[0] + A[5]*B[1] + A[6]*B[2] + A[7]*B[3];

    // C[2] = A[8]*B[0] + A[9]*B[1] + A[10]*B[2] + A[11]*B[3]
    C[2] = A[8]*B[0] + A[9]*B[1] + A[10]*B[2] + A[11]*B[3];
}

inline void matmul_3x3_3x3(const double* A, const double* B, double* C)
{
    // C = A*B
    // A is 3x3
    // B is 3x3
    // C is 3x3

    // A = [a0, a1, a2]
    //     [a3, a4, a5]
    //     [a6, a7, a8]

    // B = [b0, b1, b2]
    //     [b3, b4, b5]
    //     [b6, b7, b8]

    // Unrolled matrix multiplication
    C[0] = A[0]*B[0] + A[1]*B[3] + A[2]*B[6];
    C[1] = A[0]*B[1] + A[1]*B[4] + A[2]*B[7];
    C[2] = A[0]*B[2] + A[1]*B[5] + A[2]*B[8];

    C[3] = A[3]*B[0] + A[4]*B[3] + A[5]*B[6];
    C[4] = A[3]*B[1] + A[4]*B[4] + A[5]*B[7];
    C[5] = A[3]*B[2] + A[4]*B[5] + A[5]*B[8];

    C[6] = A[6]*B[0] + A[7]*B[3] + A[8]*B[6];
    C[7] = A[6]*B[1] + A[7]*B[4] + A[8]*B[7];
    C[8] = A[6]*B[2] + A[7]*B[5] + A[8]*B[8];
}

inline void matmul_3x3_3x1(const double* A, const double* B, double* C)
{
    // C = A*B
    // A is 3x3
    // B is 3x1

    // Unrolled matrix multiplication
    C[0] = A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
    C[1] = A[3]*B[0] + A[4]*B[1] + A[5]*B[2];
    C[2] = A[6]*B[0] + A[7]*B[1] + A[8]*B[2];
}

inline void matinv_3x3(const double* A, double* A_inv)
{
    // 3x3 matrix inversion
    double det = A[0] * (A[4] * A[8] - A[7] * A[5]) -
                 A[1] * (A[3] * A[8] - A[5] * A[6]) +
                 A[2] * (A[3] * A[7] - A[4] * A[6]);

    double inv_det = 1.0 / det;

    A_inv[0] = (A[4] * A[8] - A[7] * A[5]) * inv_det;
    A_inv[1] = (A[2] * A[7] - A[1] * A[8]) * inv_det;
    A_inv[2] = (A[1] * A[5] - A[2] * A[4]) * inv_det;
    A_inv[3] = (A[5] * A[6] - A[3] * A[8]) * inv_det;
    A_inv[4] = (A[0] * A[8] - A[2] * A[6]) * inv_det;
    A_inv[5] = (A[2] * A[3] - A[0] * A[5]) * inv_det;
    A_inv[6] = (A[3] * A[7] - A[4] * A[6]) * inv_det;
    A_inv[7] = (A[1] * A[6] - A[0] * A[7]) * inv_det;
    A_inv[8] = (A[0] * A[4] - A[1] * A[3]) * inv_det;
}



void createTransformation(const GCP& upperLeft, const GCP& upperRight, const GCP& lowerLeft, const GCP& lowerRight, double* adfGeoTransform)
{
    // Custom transformation
    // We will be solving the following equations:
    // loni = a0*xi + a1*yi + a2
    // lati = a3*xi + a4*yi + a5

    double A[] = {
        upperLeft.x, upperLeft.y, 1,
        upperRight.x, upperRight.y, 1,
        lowerLeft.x, lowerLeft.y, 1,
        lowerRight.x, lowerRight.y, 1
    };

    double AT[] = {
        upperLeft.x, upperRight.x, lowerLeft.x, lowerRight.x,
        upperLeft.y, upperRight.y, lowerLeft.y, lowerRight.y,
        1, 1, 1, 1
    };

    double b[] = {
        upperLeft.lon,
        upperRight.lon,
        lowerLeft.lon,
        lowerRight.lon
    };

    double ATA[9];
    matmul_3x4_4x3(AT, A, ATA);

    double ATb[4];
    matmul_3x4_4x1(AT, b, ATb);

    double ATA_inv[9];
    matinv_3x3(ATA, ATA_inv);

    double lon_params[3];
    matmul_3x3_3x1(ATA_inv, ATb, lon_params);


    // Now the same for latitude
    b[0] = upperLeft.lat;
    b[1] = upperRight.lat;
    b[2] = lowerLeft.lat;
    b[3] = lowerRight.lat;

    matmul_3x4_4x1(AT, b, ATb);

    double lat_params[3];
    matmul_3x3_3x1(ATA_inv, ATb, lat_params);

    // Set the transformation
    // 16 doubles
    // lon_params[0] = a0
    // lon_params[1] = a1
    // 0.0 = a2
    // lon_params[2] = a3
    // lat_params[0] = a4
    // lat_params[1] = a5
    // 0.0 = a6
    // lat_params[2] = a7
    // 0.0 = a8
    // 0.0 = a9
    // 0.0 = a10
    // 0.0 = a11
    // 0.0 = a12
    // 0.0 = a13
    // 0.0 = a14
    // 1.0 = a15
    adfGeoTransform[0] = lon_params[0];
    adfGeoTransform[1] = lon_params[1];
    adfGeoTransform[2] = 0.0;
    adfGeoTransform[3] = lon_params[2];
    adfGeoTransform[4] = lat_params[0];
    adfGeoTransform[5] = lat_params[1];
    adfGeoTransform[6] = 0.0;
    adfGeoTransform[7] = lat_params[2];
    adfGeoTransform[8] = 0.0;
    adfGeoTransform[9] = 0.0;
    adfGeoTransform[10] = 0.0;
    adfGeoTransform[11] = 0.0;
    adfGeoTransform[12] = 0.0;
    adfGeoTransform[13] = 0.0;
    adfGeoTransform[14] = 0.0;
    adfGeoTransform[15] = 1.0;
}


GTIFF_WRITE_ERROR gtiff_write(const std::string &filename, const uint8_t *data, int32_t width, int32_t height, const GCP& upperLeft, const GCP& upperRight, const GCP& lowerLeft, const GCP& lowerRight)
{
    TIFF* tif = XTIFFOpen(filename.c_str(), "w");
    if (!tif)
    {
        return GTIFF_WRITE_ERROR_OPEN_FILE;
    }

    GTIF* gtif = GTIFNew(tif);
    if (!gtif)
    {
        XTIFFClose(tif);
        return GTIFF_WRITE_ERROR_WRITE_TAGS;
    }

    GTIFKeySet(gtif, GTModelTypeGeoKey, TYPE_SHORT, 1, ModelTypeGeographic);
    GTIFKeySet(gtif, GTRasterTypeGeoKey, TYPE_SHORT, 1, RasterPixelIsArea);
    GTIFKeySet(gtif, GeographicTypeGeoKey, TYPE_SHORT, 1, GCS_WGS_84);
    GTIFKeySet(gtif, GeogCitationGeoKey, TYPE_ASCII, 7, "WGS 84");
    GTIFKeySet(gtif, GeogAngularUnitsGeoKey, TYPE_SHORT, 1, 9102);
    GTIFKeySet(gtif, GeogSemiMajorAxisGeoKey, TYPE_DOUBLE, 1, 6378137.0);
    GTIFKeySet(gtif, GeogInvFlatteningGeoKey, TYPE_DOUBLE, 1, 298.257223563);

    double adfGeoTransform[16];
    createTransformation(upperLeft, upperRight, lowerLeft, lowerRight, adfGeoTransform);

    TIFFSetField(tif, TIFFTAG_GEOTRANSMATRIX, 16, adfGeoTransform);

    // Get the resolution, by dividing the distance between the upper left and upper right corners by the width
    double xresolution = haversine(upperLeft.lon, upperLeft.lat, upperRight.lon, upperRight.lat) / width;
    double yresolution = haversine(upperLeft.lon, upperLeft.lat, lowerLeft.lon, lowerLeft.lat) / height;

    TIFFSetField(tif, TIFFTAG_XRESOLUTION, xresolution);
    TIFFSetField(tif, TIFFTAG_YRESOLUTION, yresolution);

    // Set the image dimensions
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);

    // Set the bits per sample
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);

    // Set the samples per pixel
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);

    // Set the photometric interpretation, we will use the turku color table
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_PALETTE);

    // Set the color map
    // The color map is a 256x3 array of uint16_t, stored in turku.hpp
    TIFFSetField(tif, TIFFTAG_COLORMAP, crameri_turku_colormap, crameri_turku_colormap + 256, crameri_turku_colormap + 512);

    // Write the image data
    for (int32_t y = 0; y < height; y++)
    {
        TIFFWriteScanline(tif, (void*)(data + y*width), y, 0);
    }

    GTIFWriteKeys(gtif);

    GTIFFree(gtif);
    XTIFFClose(tif);

    return GTIFF_WRITE_SUCCESS;
}

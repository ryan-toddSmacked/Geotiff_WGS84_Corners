#pragma once
#ifndef __GTIFF_GCP_HPP__
#define __GTIFF_GCP_HPP__

struct GCP
{
    double x, y;
    double lon, lat;    // Always in WGS84
};

#endif // __GTIFF_GCP_HPP__

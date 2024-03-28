
#include "haversine.hpp"
#include <cmath>


double haversine(double lon1, double lat1, double lon2, double lat2)
{
    const double R = 6371000; // metres
    const double phi1 = lat1 * M_PI / 180; // φ, λ in radians
    const double phi2 = lat2 * M_PI / 180;
    const double deltaPhi = (lat2 - lat1) * M_PI / 180;
    const double deltaLambda = (lon2 - lon1) * M_PI / 180;

    const double a = sin(deltaPhi / 2) * sin(deltaPhi / 2) +
        cos(phi1) * cos(phi2) *
        sin(deltaLambda / 2) * sin(deltaLambda / 2);
    const double c = 2 * atan2(sqrt(a), sqrt(1 - a));

    return R * c; // in metres
}

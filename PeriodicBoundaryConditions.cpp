#include "tools.h"

Coordinates periodicBoundaryCondidions(Coordinates coordinates){ // works if disk isn t moving further then 1/2 of box length
    Coordinates coordinatesNew;
    coordinatesNew.x = 0.5 * putInside(2 * coordinates.x);
    coordinatesNew.y = 0.5 * putInside(2 * coordinates.y);
    return coordinatesNew;
}

double putInside(double x){ // box [-1, 1]
        return x - 2 * (int) x;
}
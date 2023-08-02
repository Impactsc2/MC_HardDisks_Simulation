#include "tools.h"
//#include "storage.h"

Coordinates& toUnitBox(Coordinates& coordinates, Box box){
    double multiplicant = 1 / (box.dimensions[0][0] * box.dimensions[1][1] -  box.dimensions[0][1] * box.dimensions[1][0]); // 1/(DET) of box matrix

    coordinates.x = (coordinates.x * box.dimensions[1][1] - coordinates.y * box.dimensions[1][0]) * multiplicant; // to [0, 1] box "vector * matrix"
    coordinates.y = (coordinates.y * box.dimensions[0][0] - coordinates.x * box.dimensions[0][1]) * multiplicant;

    coordinates.x -= 0.5 ; // to [-0.5, 0.5] box
    coordinates.y -= 0.5;

    return coordinates;
}

Coordinates toBox(Coordinates coordinates, Box box){
    Coordinates coordinatesNew;
        
    coordinates.x += 0.5 ; // to [0, 1] box
    coordinates.y += 0.5;

    coordinatesNew.x = coordinates.x * box.dimensions[0][0] + coordinates.y * box.dimensions[1][0];
    coordinatesNew.y = coordinates.x * box.dimensions[0][1] + coordinates.y * box.dimensions[1][1];
    return coordinatesNew;
}

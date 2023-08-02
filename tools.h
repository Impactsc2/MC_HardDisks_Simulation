#pragma once
#include "storage.h"

#ifndef TOOLS_H_
#define TOOLS_H_

Coordinates& toUnitBox(Coordinates& coordinates, Box box); // changing Coordinates from Box passed in parameters to Unit Box(x and y in range (-0.5, 0,5))

Coordinates toBox(Coordinates coordinates, Box box); // changing Coordinates from Unit Box to Box passed in parameters


double generateDisplacement(double maxMove); // generate number between in range (-maxMove, maxMove)

double calculateDistance(Coordinates coordinates1, Coordinates coordinates2); // returns distance between two points

Coordinates periodicBoundaryCondidions(Coordinates coordinates); // implements periodic boundary conditions of the system

double putInside(double x);

#endif
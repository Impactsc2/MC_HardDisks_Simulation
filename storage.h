#pragma once

#ifndef STORAGE_H_
#define STORAGE_H_

struct Coordinates{
    double x;
    double y;
};

class Box{

    public:
    double dimensions[2][2];

    void setDimensions(double xx, double xy, double yx, double yy);

};

#endif
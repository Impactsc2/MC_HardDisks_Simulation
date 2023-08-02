#include "storage.h"

void Box::setDimensions(double xx, double xy, double yx, double yy){
    dimensions[0][0] = xx;
    dimensions[0][1] = xy;    
    dimensions[1][0] = yx;
    dimensions[1][1] = xx;
}
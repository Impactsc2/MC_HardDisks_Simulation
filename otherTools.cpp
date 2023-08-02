#include "tools.h"

#include <math.h>

double generateDisplacement(double maxMove){
    double x = (2 * drand48() - 1) * maxMove;
    //std::cout << x<< std::endl;
    return x; //generate displacement between (- Move) and  move
}
        
double calculateDistance(Coordinates coordinates1, Coordinates coordinates2){ // Pythagoras theorem
    return sqrt((coordinates1.x - coordinates2.x)*(coordinates1.x - coordinates2.x) + (coordinates1.y - coordinates2.y)*(coordinates1.y - coordinates2.y)); 
}
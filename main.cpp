#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <fstream>

#include "storage.h"
#include "tools.h"
#include "Disk.h"

using namespace std;

/*class Coordinates{

    public:
        double x;
        double y;
};

class Box{

    public:
    double dimensions[2][2];

    void setDimensions(double xx, double xy, double yx, double yy){
        dimensions[0][0] = xx;
        dimensions[0][1] = xy;
        dimensions[1][0] = yx;
        dimensions[1][1] = xx;
    }

};

class tools{

    public:
        double squere(double x){
            return x * x;
        }
        double generateDisplacement(double maxMove){
            double x = (2 * drand48() - 1) * maxMove;
            //std::cout << x<< std::endl;
            return x; //generate displacement between (- Move) and  move
        }
        
        double calculateDistance(Coordinates coordinates1, Coordinates coordinates2){
            return sqrt(squere(coordinates1.x - coordinates2.x) + squere(coordinates1.y - coordinates2.y)); // Pythagoras theorem
        }

        Coordinates periodicBoundaryCondidions(Coordinates coordinates){ // works if disk isn t moving further then 1/2 of box length
            Coordinates coordinatesNew;
            coordinatesNew.x = 0.5 * puttingInside(2 * coordinates.x);
            coordinatesNew.y = 0.5 * puttingInside(2 * coordinates.y);
            return coordinatesNew;
        }

        double puttingInside(double x){ // box [-1, 1]
            return x - 2 * (int) x;
        }

        Coordinates toUnitBox(Coordinates coordinates, Box box){
            double multiplicant = 1 / (box.dimensions[0][0] * box.dimensions[1][1] -  box.dimensions[0][1] * box.dimensions[1][0]); // 1/(DET) of box matrix

            Coordinates coordinatesNew;

            coordinatesNew.x = (coordinates.x * box.dimensions[1][1] - coordinates.y * box.dimensions[1][0]) * multiplicant; // to [0, 1] box "vector * matrix"
            coordinatesNew.y = (coordinates.y * box.dimensions[0][0] - coordinates.x * box.dimensions[0][1]) * multiplicant;

            coordinatesNew.x -= 0.5 ; // to [-0.5, 0.5] box
            coordinatesNew.y -= 0.5;

            return coordinatesNew;
        }

        Coordinates toBox(Coordinates coordinates, Box box){
            Coordinates coordinatesNew;
        
            coordinates.x += 0.5 ; // to [0, 1] box
            coordinates.y += 0.5;

            coordinatesNew.x = coordinates.x * box.dimensions[0][0] + coordinates.y * box.dimensions[1][0];
            coordinatesNew.y = coordinates.x * box.dimensions[0][1] + coordinates.y * box.dimensions[1][1];
            return coordinatesNew;
        }

        Coordinates full(Coordinates coordinates, Box box){
            Coordinates coordinatesNew = coordinates;
            coordinatesNew = toUnitBox(coordinatesNew, box);
            coordinatesNew = periodicBoundaryCondidions(coordinatesNew);
            coordinatesNew = toBox(coordinatesNew, box);
            return coordinatesNew;
        }


};


class Disk{

    public:
        double radius;
        Coordinates coordinates;

        void setRadius(double radiusSet){
            radius = radiusSet;
        }

};
*/
class Structure{

    public: 
        int numberOfDiscs;
        int numberOfDiscsX;
        int numberOfDiscsY;
        vector<Disk> disks;
        vector<Disk> disksUnit;


        void generateStructure(int nX, int nY, vector<vector<double>> startingPositions, vector<double> radiuses){
            for(int i = 0; i < nX * nY; i++){
                Disk diskTemp;
                diskTemp.radius = radiuses[i];
                diskTemp.coordinates.x = startingPositions[i][0];
                diskTemp.coordinates.y = startingPositions[i][1];
                disks.push_back(diskTemp);
                disksUnit.push_back(diskTemp);                    
                //std:: cout << i << "  " << diskTemp.coordinates.x << "  " << diskTemp.coordinates.y  << "    "<< disks[i].coordinates.x << "  " << disks[i].coordinates.y << std::endl;
                

            }
        }

        void initializeStructure(int nX, int nY, vector<vector<double>> startingPositions, vector<double> radiuses){
            numberOfDiscs = nX * nY; //Total number of discs
            numberOfDiscsX = nX;
            numberOfDiscsY = nY;
            generateStructure(nX, nY, startingPositions, radiuses);
        }

};
/*class MC{
    public:
        void MCSimulation(int numberOfCycles, int balancing, int numberOfParticles){
            for(int i = 0; i < numberOfCycles; i++){
                MCCycle(numberOfParticles);
            }
        }
        void MCCycle(int numberOfParticles){
            for(int i; i < numberOfParticles; i++){
                MCStep();
            }
        }
        void MCStep(){

        }
};*/

class Model : public Box{
    public:
    int acceptedMoves;
    double maxMove;
    Box box;
    Structure structure;

    Model(double maxMoveSet){
        maxMove = maxMoveSet;
        acceptedMoves = 0;
    }

    void diskUnitUpdate(){ // updating every disk position i unit box
        for(int i = 0; i < structure.numberOfDiscs; i++){
            structure.disksUnit[i].coordinates = toUnitBox(structure.disks[i].coordinates, box);
            //std::cout << structure.disksUnit[i].coordinates.x << "  " << structure.disksUnit[i].coordinates.y << std::endl;
        }
    }

    void MCSimulation(int balancing, int numberOfCycles){
        diskUnitUpdate();
        for(int i = 1; i <= numberOfCycles; i++){
            if(! (i% 1000)){
                double acceptedMovesRatio = acceptedMoves / (10.0 * structure.numberOfDiscs); // in percents
                std::cout << i << "  " << maxMove  << "  " << acceptedMovesRatio << "  " << acceptedMoves << std::endl;
                if(i < balancing){
                    if(acceptedMovesRatio > 25) maxMove *= 1.1;
                    else maxMove *= 0.9;
                }
                acceptedMoves = 0;
            }
            MCCycle();
        }
    }

    void MCCycle(){ // moving all disks one by one
        for(int i = 0; i < structure.numberOfDiscs; i++){
            MCStep(i);
        }
    }

    void MCStep(int i){ // moving single disk
        Disk diskNew = generateNewDisk(i);
        Disk diskNewUnit;

        diskNew.coordinates = toUnitBox(diskNew.coordinates, box);
        diskNew.coordinates = periodicBoundaryCondidions(diskNew.coordinates);

        diskNewUnit.coordinates = diskNew.coordinates; 

        diskNew.coordinates = toBox(diskNew.coordinates, box);

        if(isOk(diskNew, i, diskNewUnit)){ // if move is valiable update disk positions in normal box and unit box
            structure.disks[i] = diskNew;
            structure.disksUnit[i].coordinates = diskNewUnit.coordinates;
            acceptedMoves ++;
        }
        output();
    }

    vector<double> generateDeltas(){
        double deltaX = generateDisplacement(maxMove);
        double deltaY = generateDisplacement(maxMove);

        vector <double> result;

        result.push_back(deltaX);
        result.push_back(deltaY);

        return result;
    }

    Disk generateNewDisk(int i){
        vector<double> deltas = generateDeltas();
        Disk diskNew;
        diskNew.radius = structure.disks[i].radius;
        diskNew.coordinates.x = structure.disks[i].coordinates.x + deltas[0];
        diskNew.coordinates.y = structure.disks[i].coordinates.y + deltas[1];
        return diskNew;
    }

    bool isOk(Disk diskNew, int position, Disk diskNewUnit){
        for(int i = 0; i < structure.numberOfDiscs; i++){
            if(i == position) continue;
            
            if(isOverlapedV2(diskNew, diskNewUnit, i)) return false;
        }
        return true;
    }

    bool isOverlaped(Disk disk1, Disk disk2){
        //std::cout << calculateDistance(disk1.coordinates, disk2.coordinates) << "  " << disk2.radius + disk2.radius << std::endl;
        if(calculateDistance(disk1.coordinates, disk2.coordinates) < disk1.radius + disk2.radius) return true;
        
        return false;
    }

    bool isOverlapedV2(Disk disk1, Disk disk1Unit, int d2){
        Disk diskTemp = structure.disksUnit[d2];
        diskTemp.coordinates.x -= disk1Unit.coordinates.x;
        diskTemp.coordinates.y -= disk1Unit.coordinates.y;
        
        diskTemp.coordinates = periodicBoundaryCondidions(diskTemp.coordinates);

        diskTemp.coordinates.x += disk1Unit.coordinates.x;
        diskTemp.coordinates.y += disk1Unit.coordinates.y;
        diskTemp.coordinates = toBox(diskTemp.coordinates, box);
        //std::cout << diskTemp.coordinates.x << "  " << diskTemp.coordinates.y << std::endl;


        return isOverlaped(disk1, diskTemp);
    }
    
    void print(){
        for(int i = 0; i < structure.numberOfDiscs; i++){
            std::cout << structure.disks[i].coordinates.x << "  " << structure.disks[i].coordinates.y << std::endl;
        }
    }

    void output(){
        ofstream file;
        file.open("results.txt", std::ios_base::app);
        for(int i = 0; i < structure.numberOfDiscs; i++){
            file << structure.disks[i].coordinates.x << "  " << structure.disks[i].coordinates.y << std::endl;
        }
        file << std::endl;
        file.close();
    }
    

};

int main()
{
    ofstream file;
    file.open("results.txt");
    file.close();
    srand(10000);
    Model model1(1.5);
    model1.box.setDimensions(10, 0, 0, 10);
    vector<vector<double>> pos = {{1, 1}, {1, 3.5}, {1, 6}, {1, 8.5}, 
                                {3.5, 1}, {3.5, 3.5}, {3.5, 6}, {3.5, 8.5}, 
                                {6, 1}, {6, 3.5}, {6, 6}, {6, 8.5}, 
                                {8.5, 1}, {8.5, 3.5}, {8.5, 6}, {8.5, 8.5}};

    vector<vector<double>> pos1 = {{1,1} , {3,3} , {5,5} , {7,7}};

    vector<double> R = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    vector<double> R1 = {1,1,1,1};
    model1.structure.initializeStructure(4, 4, pos, R);
    model1.MCSimulation(5000, 10000);
    //model1.print();
    //std::cout << model1.isOverlapedV2(model1.structure.disksUnit[0], 0, 1) << std::endl;









    /*model1.structure.initializeStructure(2, 2, {{1, 1}, {1, 2.9}, {1, 1}, {1, 0.2}}, {1, 1, 1, 1});
    model1.box.setDimensions(10, 0, 0, 10);
    model1.diskUnitUpdate();
    std::cout << model1.structure.disksUnit[0].coordinates.y << std::endl;*/
    //std::cout << model1.isOverlapedV2(0, 1) << std::endl;
    /*Structure structure;
    structure.initializeStructure(2, 2, {{1,1}, {0,1.5}, {1, 1}, {1, 3}}, {1, 1, 1, 1});
    
    Box box = model1.box;
    Coordinates coordinates1 = structure.disks[3].coordinates;
    std::cout << coordinates1.x << "  " << coordinates1.y << std::endl;

    coordinates1 = model1.toUnitBox(coordinates1, box);
    std::cout << coordinates1.x << "  " << coordinates1.y << std::endl;

    coordinates1 = model1.periodicBoundaryCondidions(coordinates1);
    std::cout << coordinates1.x << "  " << coordinates1.y << std::endl;

    coordinates1 = model1.toBox(coordinates1, box);
    std::cout << coordinates1.x << "  " << coordinates1.y << std::endl;

    coordinates1 = model1.full(structure.disks[3].coordinates, box);
    std::cout << coordinates1.x << "  " << coordinates1.y << std::endl;
    Coordinates coordinates2;
    coordinates2.x = 0;
    coordinates2.y = 0.2;
    model1.periodicBoundaryCondidions(coordinates2);
    std::cout << coordinates2.x << "  " << coordinates2.y << std::endl;    
    */
    
}

// first you must update all Unit box coordinates
// use arithmetic rather then if's while updationg acceptedMovesRatio
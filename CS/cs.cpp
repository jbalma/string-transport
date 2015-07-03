#include <iostream>
#include "Distribution.h"


using namespace std;

int main(int argc, char **argv)
{
    int
    luser,  //the order of the legendre polynomial
    precision, //the number of points used l<precision
    g; //precision used in loops

    string kernel_name;


    bool aboutX=false;
    bool aboutY=false;
    bool aboutZ=false;
    char choice;

    Points Mapper;
    Points RotatedX, RotatedY, RotatedZ;


    Distribution Simulate;

    cout << "Preparing simulation for run..." << endl;
    cout << endl;
    cout << "Enter l:";
    cin >> luser;
    cout << endl;

    cout << "Choices for point numbers are (Be Careful):" << endl;
    cout << "6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230" << endl;
    cout << "266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730" << endl;
    cout << "2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810" << endl;
    cout << endl;

    cout << "Enter number of lebedev points (length of file): ";
    cin >> precision;
    cout << endl;





//------------------------declare objects-----------------------------//

    g=precision;

    Mapper.SetPrecisions(g); // set up point object

   
    Mapper.GenerateKernel(luser, g, "kernel"); 
	//generates a kernel file of order l 
	//and size g x g (after you generate a kernel one time, 
	//you can comment this line out and continue 
	//using it for trial runs at the same l and number of lebedev points)

    kernel_name = "kernel"; //whatever the name of the kernel file you generated with GenerateKernel() was

    Mapper.Rotate('Z', g, kernel_name); //rotates a kernel file in direction 'Z' this is the initial points object about z

    vector<double> f(g), df(g), th(g), ph(g), theta(g), phi(g), weight(g);


    //initialize rotated lists
    RotatedX.SetPrecisions(g);
    RotatedY.SetPrecisions(g);
    RotatedZ.SetPrecisions(g);

    //build rotated lits based on kernel file produced by Mapper
    RotatedX.Rotate('X',g,kernel_name); //sets up a list of rotated points about X in memory for Simulate object to use
    RotatedY.Rotate('Y',g,kernel_name); //sets up a list of rotated points about Y in memory for Simulate object to use
    RotatedZ.Rotate('Z',g,kernel_name); //no rotation occurs here


    Simulate.GenerateInitialDistribution(g);  //sets up the initial conditions and rotated points for use by simulate object
    Simulate.EvolveDistribution(g, 1000000); //runs the simulation using EvolveDistribution(precision, length of time in steps) in Distribution object



    return 0;
}

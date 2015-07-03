
#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H
#define PI 3.1415926535897932384626433832795028841971693993751058
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
//#include "gnuplot_i.hpp"
#include "Vex.h"
#include "Points.h"



class Distribution
{
    public:

        Distribution();
        void Setup(int);
        void GenerateInitialDistribution(int);

        double GetVMF(Vex, Vex, Vex, Vex, double, double);

        std::vector< std::vector<double> > GetEqDist();


        void EvolveDistribution(int, int);


        virtual ~Distribution();


    protected:


    private:

    std::vector< std::vector<double> > init_dist;

    std::vector<double> theta, phi, weight;

    double p_size;


};

#endif // DISTRIBUTION_H

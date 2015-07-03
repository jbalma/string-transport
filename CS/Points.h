#ifndef POINTS_H
#define POINTS_H
#include <vector>
#include "Vex.h"
#include <iostream>


class Points
{
    public:
        Points();

        std::string GetPrecision();

        void Setup(std::vector<double>, std::vector<double>);

        void GenerateKernel(int, int, std::string);

        std::vector<std::vector<double> > KernelMap(char rotation); //returns 2d vector of kernel for selected number of points

        std::vector<std::vector<double> > PointMap(char rotation); //returns 2d vector of lebedev points that correspond to kernel values

        int SetPrecisions(int);

        void Rotate(char, int, std::string);

        double GetWeight(int);
        double GetWeightMap(int, int);
        double GetKernel(int, int);
        double GetFkernel(int, int);
        double GetThetaPrime(int);
        double GetPhiPrime(int);

        double GetTheta(int);
        double GetPhi(int);

        double GetTheta_RotatedX(int);
        double GetPhi_RotatedX(int);

        double GetTheta_RotatedY(int);
        double GetPhi_RotatedY(int);

        std::vector<Vex> GetRotatedPoint();


        bool Check(double, double, double, std::vector<double>, std::vector<double>, std::vector<double>); //has a point been found before?

        virtual ~Points();




    protected:

    private:

    std::string precision2;

    int point_count;

    std::vector<Vex> TtoZ;

    std::vector<double> list;
    std::vector<double> Weight;

    std::vector< std::vector<double> > kernel_map;
    std::vector< std::vector<double> > weight_map;
    std::vector< std::vector<double> > fkernel_map;

    std::vector<double> Theta;
    std::vector<double> Phi;

    std::vector<double> theta_prime;
    std::vector<double> phi_prime;

    std::vector<double> theta_rotatedX;
    std::vector<double> theta_rotatedY;
    std::vector<double> theta_rotatedZ;
    std::vector<double> phi_rotatedX;
    std::vector<double> phi_rotatedY;
    std::vector<double> phi_rotatedZ;


};

#endif // POINTS_H

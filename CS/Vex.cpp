#include "Vex.h"
#include <cmath>
#include <iostream>
#define PI 3.1415926535897932384626433832795028841971693993751058


Vex::Vex()
{
    //ctor
}

void Vex::Setup(double my_theta, double my_phi, double my_r)
{

    Vex::theta=my_theta;
    Vex::phi=my_phi;
    Vex::r=my_r;

    Vex::x = my_r*cos(my_phi)*sin(my_theta);
    Vex::y = my_r*sin(my_phi)*sin(my_theta);
    Vex::z = my_r*cos(my_theta);



}

void Vex::SetXYZ(double my_x, double my_y, double my_z)
{
    Vex::x=my_x;
    Vex::y=my_y;
    Vex::z=my_z;

    Vex::r = std::sqrt( pow(my_x,2) + pow(my_y,2) + pow(my_z,2) );
}

void Vex::SetSpherical(double my_theta, double my_phi, double my_r)
{
    Vex::theta=my_theta;
    Vex::phi=my_phi;
    Vex::r=my_r;
}

void Vex::XYZtoSpherical(double my_x, double my_y, double my_z)
{


    Vex::r = std::sqrt( pow(my_x,2) + pow(my_y,2) + pow(my_z,2) );

    //Vex::r=1;


    Vex::theta = acos(my_z / Vex::r);
    Vex::phi = atan2(my_y,my_x);

}


void Vex::SphericaltoXYZ(double my_theta, double my_phi, double my_r)
{

    Vex::x = my_r*cos(my_phi)*sin(my_theta);
    Vex::y = my_r*sin(my_phi)*sin(my_theta);
    Vex::z = my_r*cos(my_theta);

}


double Vex::GetPhi()
{
    return Vex::phi;
}

double Vex::GetTheta()
{
    return Vex::theta;
}

double Vex::GetR()
{
    return Vex::r;
}

double Vex::GetX()
{
    return Vex::x;
}

double Vex::GetY()
{
    return Vex::y;
}

double Vex::GetZ()
{
    return Vex::z;
}



Vex Vex::operator^(Vex myvex) //cross product returns vector
{
    double my_x, my_y, my_z;
    Vex NewVector;

    my_x= Vex::y * myvex.GetZ() - Vex::z * myvex.GetY();
    my_y= Vex::z * myvex.GetX() - Vex::x * myvex.GetZ();
    my_z= Vex::x * myvex.GetY() - Vex::y * myvex.GetX();



    NewVector.SetXYZ(my_x, my_y, my_z);

    return(NewVector);

}

double Vex::operator<(Vex myvex) //dot product returns scalar
{
    double dot_product;
    double my_x=0, my_y=0, my_z=0;

    my_x=Vex::x*myvex.GetX();
    my_y=Vex::y*myvex.GetY();
    my_z=Vex::z*myvex.GetZ();

    dot_product=my_x+my_y+my_z;


    return dot_product;
}

Vex Vex::operator*(double constant)
{
    Vex NewVector;
    double my_x, my_y, my_z;

    my_x=Vex::x*constant;
    my_y=Vex::y*constant;
    my_z=Vex::z*constant;


    NewVector.SetXYZ(my_x, my_y, my_z);

    return(NewVector);

}

Vex Vex::operator-(Vex myvex)
{
    Vex NewVector;
    double my_x, my_y, my_z;

    my_x=Vex::x-myvex.GetX();
    my_y=Vex::y-myvex.GetY();
    my_z=Vex::z-myvex.GetZ();

    NewVector.SetXYZ(my_x, my_y, my_z);

    return(NewVector);
}

Vex Vex::operator+(Vex myvex)
{
    Vex NewVector;
    double my_x, my_y, my_z;

    my_x=Vex::x+myvex.GetX();
    my_y=Vex::y+myvex.GetY();
    my_z=Vex::z+myvex.GetZ();

    NewVector.SetXYZ(my_x, my_y, my_z);

    return(NewVector);
}

Vex Vex::operator/(double constant)
{
    Vex NewVector;
    double my_x, my_y, my_z;

    my_x=Vex::x/constant;
    my_y=Vex::y/constant;
    my_z=Vex::z/constant;

    NewVector.SetXYZ(my_x, my_y, my_z);

    return(NewVector);
}


Vex Vex::operator-(double constant) //subtract scalar
{
    Vex NewVector;
    double my_x, my_y, my_z;

    my_x=Vex::x-constant;
    my_y=Vex::y-constant;
    my_z=Vex::z-constant;

    NewVector.SetXYZ(my_x, my_y, my_z);

    return(NewVector);
}

Vex Vex::operator+(double constant) //add scalar
{
    Vex NewVector;
    double my_x, my_y, my_z;

    my_x=Vex::x+constant;
    my_y=Vex::y+constant;
    my_z=Vex::z+constant;

    NewVector.SetXYZ(my_x, my_y, my_z);

    return(NewVector);
}

bool Vex::operator==(Vex myvex)
{
    double TOL=1e-17;
    if( (std::abs(Vex::x-myvex.GetX())<=TOL) && (std::abs(Vex::y-myvex.GetY())<=TOL) && (std::abs(Vex::z-myvex.GetZ())<=TOL) )
    {
        return true;
    }else{return false;}
}

Vex::~Vex()
{
    //dtor
}

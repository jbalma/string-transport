#include "Points.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include "Vex.h"
#include <vector>
#include <iomanip>
#define PI 3.1415926535897932384626433832795028841971693993751058
#include "Legendre.h"


Points::Points()
{
    //ctor
}

int Points::SetPrecisions(int precision)
{
    int g=precision;
    double con=PI/180.0;


    ///-----------build arrays based on precision of selected lebedev points
    //std::vector<double> phi(g), theta(g), weight(g);
    //std::vector<double> phi_rotatedX(g), theta_rotatedX(g), theta_rotatedY(g), phi_rotatedY(g), theta_rotatedZ(g), phi_rotatedZ(g);

    std::vector<double> phi_clean(g), theta_clean(g), weight_clean(g);
    //std::vector<double> thh(g), phh(g);

    Points::Theta.resize(g);
    Points::Phi.resize(g);
    Points::Weight.resize(g);

    Points::theta_rotatedX.resize(g);
    Points::phi_rotatedX.resize(g);

    Points::theta_rotatedY.resize(g);
    Points::phi_rotatedY.resize(g);

    Points::theta_rotatedZ.resize(g);
    Points::phi_rotatedZ.resize(g);


    Points::point_count=precision;



    ///---------get lebedev point file based on precision ------------///

    if (precision == 6)
    {
        Points::precision2="003";
    }
    else if (precision == 14)
    {
        Points::precision2="005";
    }
    else if (precision == 26)
    {
        Points::precision2="007";
    }
    else if (precision == 38)
    {
        Points::precision2="009";
    }
    else if (precision == 50)
    {
        Points::precision2="011";
    }
    else if (precision == 74)
    {
        Points::precision2="013";
    }
    else if (precision == 86)
    {
        Points::precision2="015";
    }
    else if (precision == 110)
    {
        Points::precision2="017";
    }
    else if (precision == 146)
    {
        Points::precision2="019";
    }
    else if (precision == 170)
    {
        Points::precision2="021";
    }
    else if (precision == 194)
    {
        Points::precision2="023";
    }
    else if (precision == 230)
    {
        Points::precision2="025";
    }
    else if (precision == 266)
    {
        Points::precision2="027";
    }
    else if (precision == 302)
    {
        Points::precision2="029";
    }
    else if (precision == 350)
    {
        Points::precision2="031";
    }
    else if (precision == 434)
    {
        Points::precision2="035";
    }
    else if (precision == 590)
    {
        Points::precision2="041";
    }
    else if (precision == 770)
    {
        Points::precision2="047";
    }
    else if (precision == 974)
    {
        Points::precision2="053";
    }
    else if (precision == 1202)
    {
        Points::precision2="059";
    }
    else if (precision == 1454)
    {
        Points::precision2="065";
    }
    else if (precision == 1730)
    {
        Points::precision2="071";
    }
    else if (precision == 2030)
    {
        Points::precision2="077";
    }
    else if (precision == 2354)
    {
        Points::precision2="083";
    }
    else if (precision == 2702)
    {
        Points::precision2="089";
    }
    else if (precision == 3074)
    {
        Points::precision2="095";
    }
    else if (precision == 3470)
    {
        Points::precision2="101";
    }
    else if (precision == 3890)
    {
        Points::precision2="107";
    }
    else if (precision == 4334)
    {
        Points::precision2="113";
    }
    else if (precision == 4802)
    {
        Points::precision2="119";
    }
    else if (precision == 5294)
    {
        Points::precision2="125";
    }
    else if (precision == 5810)
    {
        Points::precision2="131";
    }
    else
    {
        Points::precision2="003";
        std::cout << precision <<" is not an option! Setting default (003)!" << std::endl;
    }



//-----------------------get lebedev points from files------------------------------------//


    std::ifstream infile;
    std::string filename1, filename2, filename3;
    filename1="lebedev_";
    filename2=precision2;
    filename3=".txt";


    int num = 0; // num must start at 0

    infile.open(filename1+filename2+filename3);// file containing numbers in 3 columns
    infile.precision(15);

    if(infile.fail()) // checks to see if file opended
    {
        std::cout << "error, no file!" << std::endl;
        return 1; // no point continuing if the file didn't open...
    }


    while(!infile.eof() && num < precision) // reads file to end of *file*, not line
    {

        infile >> phi_clean[num]; // from -pi to pi
        infile >> theta_clean[num]; // from 0 to pi
        infile >> weight_clean[num]; // weight column


        //-------arrays of lebedev points--------//

        Points::Phi[num]=phi_clean[num]*con;
        Points::Theta[num]=theta_clean[num]*con;
        Points::Weight[num] = weight_clean[num];


        ///setup lists of rotated points for use in kernel rotation
        //-----------------ThetaX / PhiX -------------//

        Vex VexX;
        VexX.SphericaltoXYZ(Points::Theta[num], Points::Phi[num], 1.0);
        VexX.SetXYZ(VexX.GetZ(), VexX.GetY(), -1*VexX.GetX());
        VexX.XYZtoSpherical(VexX.GetX(), VexX.GetY(), VexX.GetZ());

        Points::theta_rotatedX[num] = VexX.GetTheta();
        Points::phi_rotatedX[num] = VexX.GetPhi();


        //--------------ThetaZ / PhiZ ----------------//

        Vex VexZ;
        VexZ.SetSpherical(Points::Theta[num], Points::Phi[num], 1.0);

        Points::theta_rotatedZ[num] = VexZ.GetTheta();
        Points::phi_rotatedZ[num] = VexZ.GetPhi();


        //------------- ThetaY / PhiY ---------------//

        Vex VexY;
        VexY.SphericaltoXYZ(Points::Theta[num], Points::Phi[num], 1.0);
        VexY.SetXYZ(VexY.GetX(), VexY.GetZ(), -1*VexY.GetY());
        VexY.XYZtoSpherical(VexY.GetX(), VexY.GetY(), VexY.GetZ());


        theta_rotatedY[num] = VexY.GetTheta();
        phi_rotatedY[num] = VexY.GetPhi();


        //std::cout << "object rotated: " << Points::theta_rotatedZ[num] << ",    " << Points::phi_rotatedZ[num]<< ",   " << std::endl;
        //std::cout << "object original: " << Points::Theta[num] << ",    " << Points::Phi[num] << std::endl;

        ++num;
    }
    infile.close();

    return 0;

}

void Points::GenerateKernel(int L, int point_count, std::string title)
{

    int g=point_count, m=0;

    LEGENDRE P_lm;
    LEGENDRE Y_P;
    LEGENDRE dP_lm;

    std::vector<std::vector<double> > d_kern(g); //kernel for derivative reconstruction
    std::vector<std::vector<double> > f_kern(g); //kernel for function (test) reconstruction
    std::vector<std::vector<double> > WT(g);

    std::complex<double> Y(0,0), Ylm(0,0), dYlm(0,0), Ymp1(0,0), ej(0,0), function(0,0), derivative(0,0);
    std::complex<double> im(0,1);

    double th1=0, ph1=0, sign=0;



    std::cout << title << std::endl;
    std::ofstream kernel(title);
    kernel.precision(15);




    for(int i=0; i<g; i++)
    {
        d_kern[i].resize(g);
        f_kern[i].resize(g);
        WT[i].resize(g);

        for(int j=0; j<g; j++)
        {
            for(double l=0; l<=L; l++)
            {
                for(double m_it=0; m_it<=(2*l); m_it++)
                {

                    m=0;

                    m = l-m_it;

                    //std::cout << "m = " <<  m << ", l = " << l <<  std::endl;

                    ej = m;
                    sign=pow(-1.0,m_it);

                    std::complex<double> exponential_prime(cos( Points::Phi[i]), (-1)*sin(Points::Phi[i]));
                    std::complex<double> exponential(cos(m*Points::Phi[j]), sin(m*Points::Phi[j]));


                    Ylm = P_lm.Yml(m, l, Points::Theta[i], Points::Phi[i]);
                    Y = Y_P.Yml(m, l, Points::Theta[j], Points::Phi[j]);

                    if( Theta[i] != 0 && ((m+1)<=l) )
                    {
                        Ymp1 = m * (1.0/tan(Points::Theta[i])) * dP_lm.Yml(m, l, Points::Theta[i], Points::Phi[i]) + sqrt( (l-m)*(l+m+1) ) * exponential_prime * dP_lm.Yml(m+1, l, Points::Theta[i], Points::Phi[i]);

                    }

                    ///fill arrays with f=Y*Y for the function kernel and derivative kernel

                    f_kern[i][j] += (conj(Y)*Ylm).real();//Y_real*Y_prime_real;
                    d_kern[i][j] += (conj(Y)*Ymp1).real();

                }

            }

            ///absorb weights into kernel

            WT[i][j] = Points::Weight[j]*4.0*PI;
            kernel << d_kern[i][j]*Points::Weight[j]*4.0*PI  << "      " << f_kern[i][j]*Points::Weight[j]*4.0*PI << "      " << WT[i][j] << std::endl;

        }
    }

    kernel.close();

}





void Points::Rotate(char axis, int point_count, std::string title)
{
    int g=point_count;

    std::vector<int> index(g*g);
    std::vector<int> rindex(g*g);
    std::vector<int> It(g);
    std::vector<double> theta(g);
    std::vector<double> phi(g);
    std::vector<double> theta_rotated(g);
    std::vector<double> phi_rotated(g);



    std::vector<double> theta_index(g*g);
    std::vector<double> phi_index(g*g);
    Points::kernel_map.resize(g);
    Points::weight_map.resize(g);
    Points::fkernel_map.resize(g);
    Points::theta_prime.resize(g*g);
    Points:: phi_prime.resize(g*g);

    std::vector<std::vector<double> > kern(g);
    std::vector<std::vector<double> > fkern(g);
    std::vector<std::vector<double> > weight(g);
    std::vector<std::vector<double> > rkern(g);
    std::vector<std::vector<double> > rkern2(g);
    std::vector<std::vector<double> > rotated_kern(g);
    std::vector<double> r_weight(g*g);
    std::vector<std::vector<double> > rfkernel_value(g);

    std::vector<Vex> Original(g);
    std::vector<Vex> Rotated(g);
    Points::TtoZ.resize(g);

    theta=Points::Theta;
    phi=Points::Phi;






    if(axis=='X')
    {
        theta_rotated=Points::theta_rotatedX;
        phi_rotated=Points::phi_rotatedX;

    }
    else if(axis=='Y')
    {

        theta_rotated=Points::theta_rotatedY;
        phi_rotated=Points::phi_rotatedY;

    }
    else if(axis=='Z')
    {

        theta_rotated=Points::theta_rotatedZ;
        phi_rotated=Points::phi_rotatedZ;

    }
    else
    {

    }


    std::vector<int> i_it(g*g);
    std::vector<int> j_it(g*g);


      for(int i=0; i<g; i++)
    {

        Original[i].SphericaltoXYZ(theta[i],phi[i],1.0);
        Rotated[i].SphericaltoXYZ(theta_rotated[i], phi_rotated[i], 1.0);

       // std::cout <<"x, y, z =" << Original[i].GetX() << ", " << Original[i].GetY() << ", "<< Original[i].GetZ() << ", "<< std::endl;
    }


    for(int i=0; i<g; i++)
    {
        kern[i].resize(g);
        fkern[i].resize(g);
        weight[i].resize(g);

        rkern[i].resize(g);
        rkern2[i].resize(g);

        //std::cout << "original xyz = ("  << Original[i].GetX() << ", " << Original[i].GetY() << ", " << Original[i].GetZ() << ")" << std::endl;


    }



//------------ load unroated kernel matrix from file ------------------//


    std::cout <<axis+title << std::endl;
    std::ifstream kernel_read(title);
    kernel_read.precision(15);
    kernel_read.clear();                 // clear fail and eof bits
    kernel_read.seekg(0, std::ios::beg); // back to the start!

    int u=0;

    for(int i=0; i<g; i++)
    {
        kernel_map[i].resize(g);
        fkernel_map[i].resize(g);
        weight_map[i].resize(g);

        for(int j=0; j<g; j++)
        {
            //kernel_read >> theta_index[u] >> phi_index[u] >> kern[i][j] >> fkern[i][j] >> weight[i][j];
            kernel_read >> kern[i][j] >> fkern[i][j] >> weight[i][j];
            u++;
        }

    }

    kernel_read.close();


//--------------find index that matches rotation -----------------------//


    for(int i=0; i<g; i++)
    {
        for(int j=0; j<g; j++)
        {
            if(Original[i]==Rotated[j])
            {
                Points::TtoZ[i] = Rotated[j];
                It[i]=j;
            }
        }
    }





//-------------------use rotated index to map kernel ------------------------//

    for(int i=0; i<g; i++)
    {
        Points::theta_prime[i]=theta[ It[i] ];
        Points::phi_prime[i]=phi[ It[i] ];

        for(int j=0; j<g; j++)
        {
            Points::kernel_map[i][j]=kern[ It[i] ][ It[j] ];
            Points::fkernel_map[i][j]=fkern[ It[i] ][ It[j] ];
            Points::weight_map[i][j]=weight[ It[i] ][ It[j] ];

        }

    }




    std::cout << std::endl;
    std::cout << "Rotation Passed! " <<  std::endl;
    std::cout << std::endl;

}


//-----------------Rotated Point object-----------------------------//

std::vector<Vex> Points::GetRotatedPoint()
{
    return Points::TtoZ;
}




//----------------kernel and weight accessors------------------------//

double Points::GetWeightMap(int i, int j)
{
    return Points::weight_map[i][j];
}

double Points::GetKernel(int i, int j)
{
    return Points::kernel_map[i][j];
}

double Points::GetFkernel(int i, int j)
{
    return Points::fkernel_map[i][j];
}







//----------------Points Accessors----------------//

double Points::GetThetaPrime(int i)
{
    return Points::theta_prime[i];
}

double Points::GetPhiPrime(int i)
{
    return Points::phi_prime[i];
}

double Points::GetTheta(int i)
{
    return Points::Theta[i];
}

double Points::GetPhi(int i)
{
    return Points::Phi[i];
}

double Points::GetWeight(int i)
{
    return Points::Weight[i];
}

double Points::GetTheta_RotatedX(int i)
{
    return Points::theta_rotatedX[i];
}

double Points::GetPhi_RotatedX(int i)
{
    return Points::phi_rotatedX[i];
}

double Points::GetTheta_RotatedY(int i)
{
    return Points::theta_rotatedY[i];
}
double Points::GetPhi_RotatedY(int i)
{
    return Points::phi_rotatedY[i];
}

std::string Points::GetPrecision()
{
    return Points::precision2;
}
Points::~Points()
{
    //dtor
}

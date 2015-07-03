#include "Distribution.h"
#include "gnuplot_i.hpp"
#include <random>
#include <iomanip>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#define PI 3.1415926535897932384626433832795028841971693993751058
#include <sys/resource.h>

Distribution::Distribution()
{
    //ctor
}

void Distribution::Setup(int point_count)
{


}

void Distribution::GenerateInitialDistribution(int g)
{
    //Generate Mise-Fisher Distribution of p_size gxg
    Distribution::init_dist.resize(g);
    Distribution::p_size=g;



    for(int i=0; i<g; i++)
    {
        Distribution::init_dist[i].resize(g);

        for(int j=0; j<g; j++)
        {
            Distribution::init_dist[i][j] = 0; //( (1.0/(sigma*sqrt(2.0*PI))) * exp( -1*pow((j/gg-mean), 2) / (2*pow(sigma,2)) ));
            //std::cout << init_dist[i][j] << ", ";
        }
        // std::cout << std::endl;
    }



}


//GenDist takes in A(theta,phi) and uses the resulting xyz coords to create dstribution
double Distribution::GetVMF(Vex A, Vex B, Vex a, Vex b, double alpha, double beta)
{


    double lower;
    double upper;

    double p;
    double P;

    //p(A)=(1 / (4pi * alpha * sinh(1/alpha))) * e^(u*A/alpha)

    P = std::abs( ((a+b)-(A+B))<( (b-a)^(B-A) ) );


    //std::cout << P << std::endl;

    ///----------------P(A)----------------------///

    //upper = exp(((a<A) / alpha));
    //lower=(PI * 4.0 * alpha * sinh(1.0/alpha));

    ///----------------P(B)----------------------///

    //upper = exp(((b<B) / beta));
    //lower=(PI * 4.0 * beta * sinh(1.0/beta));

    ///-------------Mixed AB--------------------///

    upper = exp( -( ((a<A)/alpha) + ((b<B) / beta)  + P));

    lower = (PI * 4.0 * alpha * sinh(1.0/alpha)) * (PI * 4.0 * beta * sinh(1.0/beta));





    //p = upper/lower; //probablity density of element of von-misFisher distribution
    p=exp(0); //flat distribution



    return p;
}




void Distribution::EvolveDistribution(int point_count, int t) //This function takes in size of arrays, and total simulation time
{

    bool IncludeP=false, plot_ph=false, plot_sphere=false, save_animation=false, save_entropy=false, leave=false;
    bool solve_P=false, plot=false, sphere=false, animation=false, ent=false;
    double dt;



        ///----------------Choose what to plot--------------------///

    while(!leave)
    {
        leave=false;


        std::cout << "step size: ";
        std::cin >> dt;
        std::cout << std::endl;

        std::cout << "Plot P(A) and P(B) on sphere (1 = yes, 0 = no): ";
        std::cin >> plot_sphere;
        std::cout << std::endl;

        std::cout << "Show entropy and energy density? (1 = yes, 0 = no): ";
        std::cin >> save_entropy;
        std::cout << std::endl;






        if (plot_sphere == 1)
        {
            sphere=true;
            leave=true;
        }
        else
        {
            sphere=false;
            leave=true;
        }


        if (save_entropy==1)
        {

            ent=true;
            leave=true;
        }
        else
        {
            ent=false;
            leave=true;
        }

    }

//----------------------------declare variables------------------------------------------//

    int g=point_count;


    Gnuplot g0, g1, g2("lines"), g3, g4, g5, g6, g7, g8, g9, g11("lines"); //gnuplot objects


// arrays for gravitational term
    std::vector<double> dAX(g);
    std::vector<double> dAY(g);
    std::vector<double> dAZ(g);

    std::vector<double> dBX(g);
    std::vector<double> dBY(g);
    std::vector<double> dBZ(g);



    double energy=1;  //energy density
    double E0=0; //initial energy density
    double slope = 0; //slope of the energy density
    double timer=0; //the simulation time
    double P=0; //wedge product terms


//arrays for holding data to plot
    std::vector<double> Entropy(t);
    std::vector<double> EntropyLimit(t);
    std::vector<double> Time(t);
    std::vector<double> Energy(t);
    std::vector<double> Temperature(t);
    std::vector<double> Falloff(t);
    std::vector<double> Slope(t);

    std::vector<double> w(g);



    Points PointSet;
    PointSet.SetPrecisions(g);

    Vex First, Second, Third, Fourth;
    Vex a, b;
    Vex V_last;

    std::vector<Vex> Point(g);


    Distribution::theta.resize(g);
    Distribution::phi.resize(g);
    Distribution::weight.resize(g);

    std::ofstream myfile;
    std::ofstream myfile2;

    Distribution::p_size=g;

    //--------------------end variable declaration----------------------//



///-----------------------------------------------------------------------------///




    std::vector< std::vector<double> > f(g);
    std::vector< std::vector<double> > df(g);
    std::vector< std::vector<double> > df_small(g);
    std::vector< std::vector<double> > df_large(g);
    std::vector< std::vector<double> > df_long(g);
    std::vector< std::vector<double> > dfg(g);
    std::vector< std::vector<double> > dA(g);
    std::vector< std::vector<double> > dB(g);




    for(int i=0; i<g; i++)
    {
        dA[i].resize(g);
        dB[i].resize(g);
        f[i].resize(g);
        df[i].resize(g);
        dfg[i].resize(g);
        df_small[i].resize(g);
        df_large[i].resize(g);
        df_long[i].resize(g);

    }


    //----------------------setup lebedev points and vectors--------------//

    for(int i=0; i<g; i++)
{
        Distribution::theta[i] = PointSet.GetTheta(i);
        Distribution::phi[i] = PointSet.GetPhi(i);
        Distribution::weight[i] = PointSet.GetWeight(i);
        w[i] = PointSet.GetWeight(i);
        Point[i].Setup(PointSet.GetTheta(i), PointSet.GetPhi(i), 1.0);

        //std::cout << "A(x,y,z) = " << "(" << Point[i].GetTheta() << ", " << Point[i].GetPhi() << ", " << Point[i].GetZ() << ")" << std::endl;
    }









///---------------------------------Setup derivative terms ---------------------------------------//
    Points RotatedX;
    Points RotatedY;
    Points RotatedZ;

    RotatedX.SetPrecisions(g);
    RotatedY.SetPrecisions(g);
    RotatedZ.SetPrecisions(g);

    RotatedX.Rotate('X', g, "kernel");
    RotatedY.Rotate('Y', g, "kernel");
    RotatedZ.Rotate('Z', g, "kernel");


///------------------------------------------------------------------------------///













///--------------------------SETUP VMF DISTRIBUTION-------------------------------///

    std::cout << std::endl;
    std::cout << std::endl;

    double alpha = 0.8;
    double beta = 0.5;

    //a.Setup(Point[10].GetTheta(), Point[10].GetPhi(), 1.0);
    //b.Setup(Point[10].GetTheta(), Point[10].GetPhi(), 1.0);
    //b.Setup(Point[3].GetTheta(), Point[3].GetPhi(), 1.0);
    a.Setup(3*PI/4.0, PI/4.0, 1.0);
    b.Setup(PI/3.0, 5*PI/4.0, 1.0);


    for (int i=0; i<g; i++)
    {
        for(int j=0; j<g; j++)
        {
            Distribution::init_dist[i][j] = Distribution::GetVMF(Point[i], Point[j], a, b, alpha, beta);
            f[i][j] = Distribution::GetVMF(Point[i], Point[j], a, b, alpha, beta); //fill the distribution with VMF function output

        }
    }

    for (int i=0; i<g; i++)
    {
        for(int j=0; j<g; j++)
        {
            E0 += f[i][j]*w[i]*w[j];  //calculate initial energy density of distribution for normalization
        }
    }



    std::cout << "Initial energy rho: " << E0 << std::endl;
    std::cout << std::endl;

    for (int i=0; i<g; i++)
    {
        for(int j=0; j<g; j++)
        {
            f[i][j] = f[i][j]/E0; //normalize it.

        }
    }


///--------------------------------------------------------------------------------///






///-----------------------Generate P---------------------------------------///

	///This term only needs to be computed once. So do it here before the time loop

    std::vector< std::vector< std::vector< std::vector<double> > > > P_list(g);

    for(int i=0; i<g; i++)
    {
        P_list[i].resize(g);

        for(int j=0; j<g; j++)
        {
            P_list[i][j].resize(g);

            for(int k=0; k<g; k++)
            {
                P_list[i][j][k].resize(g);
            }
        }
    }



    for(int A=0; A<g; A++)
    {
        for(int B=0; B<g; B++)
        {
            for(int Ap=0; Ap<g; Ap++)
            {
                for(int Bp=0; Bp<g; Bp++)
                {
                    P_list[A][B][Ap][Bp] = (1.0/4.0)*std::abs( ((Point[Ap]+Point[Bp])-(Point[A]+Point[B]))<( (Point[Bp]-Point[Ap])^(Point[B]-Point[A]) ) );

                }
            }
        }
    }
    std::cout << std::endl;
    std::cout << std::endl;
///---------------------------------------------------------------------------------------------///


    double result=0;
    double avg_slope=0;

///some tests to verify vector math works, as well as lebedev quadrature

double X=0;
double Y=0;
double Z=0;

    for(int i=0; i<g; i++)
    {
        X=sin(Point[i].GetTheta())*cos(Point[i].GetPhi());
        Y=sin(Point[i].GetTheta())*sin(Point[i].GetPhi());
        Z=cos(Point[i].GetTheta());

        result += 4.0*PI * pow(cos(Point[i].GetTheta()),4)*pow(sin(Point[i].GetPhi()),4) * w[i];
    }

    std::cout << "Here is the integral f(theta,phi): " << result << std::endl;


    double mm=0;
    double counter=0;

/*
    for (int A=0; A<g; A++)
    {
        for(int B=0; B<g; B++)
        {
            for(int Ap=0; Ap<g; Ap++)
            {
                for(int Bp=0; Bp<g; Bp++)
                {
                    //std::cout << std::abs(((Point[B]+Point[Bp]-Point[A]-Point[Ap]).GetR())) << std::endl;

                    if(((Point[B]+Point[Bp]-Point[A]-Point[Ap]).GetR())==0)
                    {
                        // std::cout << std::cout << std::abs(((Point[B]+Point[Bp]-Point[A]-Point[Ap]).GetR())) << std::endl;
                        counter++;
                    }

                }
            }
        }
    }

*/


    mm=counter/(g*g*g*g);

    std::cout << "counter: " << mm << std::endl;
    std::cout << counter << std::endl;
    std::cout << g*g*g*g << std::endl;


    Vex R,O;

    R.SetXYZ(0.1,0.2,0);
    O.SetXYZ(1,0.5,0.2);

    std::cout << (R^O).GetX() << ", " << (R^O).GetY() << ", "<< (R^O).GetZ() << std::endl;


    std::cout << "Press Enter to start..." << std::endl;
    std::cin.get();
    std::cin.get();


///-------------start timer------------------------///
    std::clock_t start1 = clock();
///------------------------------------------------///

    double P_A=0;
    double P_B=0;
    double U2=0, V2=0;


    std::vector<double> myRow(g), myCol(g), myRowsum(g), myColsum(g), pA(g), pB(g);
    Vex U, V, Abar, Bbar;



    double V2bar = 0;
    double U2bar=0;
    Vex Ubar;
    Vex Vbar;

    double t0=10.0;

    double large_loops=0;
    double small_loops=0;
    double redshift=0;

    double c=1;
    double d=1;

    double cl=1;

    int count=0;



///----------------start time loop------------------------///

    for(int s=1; s<t; s++)
    {

        timer = t0 + (s-1.0)*dt;

        double inter_prob = 1; //intercommutation probabiliy
        double rho0=0;
        double CUTOFF=1e-12; //size of small loops to remove


//initial energy density

        for(int i=0; i<g; i++)
        {
            for(int j=0; j<g; j++)
            {
                rho0 +=f[i][j]*w[i]*w[j];
            }
        }


        // Falloff[s] = 1.0/(timer);//-2.1*(-1+E0)/(E0*timer);

///----------------------------Terminal Output ----------------------------///

        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "//////////////////////////////////////////" << std::endl;
        std::cout << "//" << std::endl;
        std::cout << "// time = " << timer << " / " << t0 + dt*t << std::endl;
        std::cout << "//" << std::endl;
        std::cout << "// Energy Density: " << Energy[s-1] << std::endl;
        std::cout << "// Entropy: " << Entropy[s-1] << std::endl;
        std::cout << "// Entropy Limit: " << EntropyLimit[s-1] << std::endl;
        std::cout << "// Slope of log(rho) vs log(t): " << Slope[s-1] << std::endl;
        std::cout << "// Average slope: " << avg_slope << std::endl;
        std::cout << "// log(rho): " << log(Energy[s-1]) << std::endl;
        std::cout << "//" << std::endl;
        std::cout << "//" << std::endl;
        //std::cout << "// P(A), P(B) : " << P_A << ", " << P_B << std::endl;
        std::cout << "//" << std::endl;
        std::cout << "// <U>^2, <V>^2 : " << (U<U) << ", " << (V<V) << ")" << std::endl;
        std::cout << "//" << std::endl;
        std::cout << "// <U^2>, <V^2> : (" << U2 << ", " << V2 << ")" <<  " <u^2> + <v^2> = " << U2+V2 << std::endl;
        std::cout << "//" << std::endl;
        std::cout << "// Average A and B: " << Abar.GetR() << ", " << Bbar.GetR() << std::endl;
        std::cout << "//" << std::endl;
        std::cout << "// loops ratio (small/large) :" << small_loops/large_loops << ", " << (double)d/(double)c << std::endl;
        std::cout << "//" << std::endl;
        std::cout << "//" << std::endl;
        std::cout << "// t^2 rho: " << rho0*pow(timer,2) << std::endl;
        std::cout << "// d = " << d << std::endl;
        std::cout << "//////////////////////////////////////////" << std::endl;
        std::cout << std::endl;



///-------------start timer------------------------///
        std::clock_t start2 = clock(); //time integral computation
///------------------------------------------------///


        double term=0;
        double random1=(double)rand()/RAND_MAX; //normalized random number





        for (int A=0; A<g; A++)
        {
            for(int B=0; B<g; B++)
            {
                df[A][B] = 0;
                df_large[A][B]=0;
                df_small[A][B]=0;
                df_long[A][B]=0;

                for(int Ap=0; Ap<g; Ap++)
                {
                    for(int Bp=0; Bp<g; Bp++)
                    {
                        P = P_list[A][B][Ap][Bp];

                        //df[A][B] += c * sqrt(rho0) * (1.0/rho0)*( (f[Ap][B])*(f[A][Bp]) - (f[A][B])*(f[Ap][Bp]) ) * (1.0 + rho0*inter_prob*P) * w[Ap] * w[Bp] * dt; //no loop removal, comment out everything after this line in inner loop to use this

                        ///remove small loops

                        if((Point[B]+Point[Bp]-Point[A]-Point[Ap]).GetR()<=CUTOFF)
                        {
                            term=0;
                            df_long[A][B] += (1.0/rho0)*( (f[Ap][B])*(f[A][Bp]) ) * (1.0/cl) * w[Ap] * w[Bp] * dt;

                        }
                        else
                        {
                            term=f[Ap][B]*f[A][Bp];

                        }

                        df_small[A][B] += (1.0/rho0)*( term - (f[A][B]) * (f[Ap][Bp])) * c * (1.0/cl) * w[Ap] * w[Bp] * dt; //df/dt = dA'dB' P(A,B,A'B') * [ f(A',B) f(A,B') - f(A,B) f(A',B') ]

                        ///remove large loops


                        if((Point[B]+Point[Bp]-Point[A]-Point[Ap]).GetR()>CUTOFF)
                        {
                            term=0;
                            df_long[A][B] += ( (f[Ap][B])*(f[A][Bp]) ) * (inter_prob*P)* (cl) * w[Ap] * w[Bp] * dt;

                        }
                        else
                        {
                            term=f[Ap][B]*f[A][Bp];
                        }

                        df_large[A][B] += (1.0/rho0)*( term - (f[A][B]) * (f[Ap][Bp])) *pow(inter_prob,1) * P *  d * (rho0*cl) * w[Ap] * w[Bp] * dt; //df/dt = dA'dB' P(A,B,A'B') * [ f(A',B) f(A,B') - f(A,B) f(A',B') ]


                    } //end j'
                } // end i'
            } //end j
        } //end i


///------------end timer-------------------------///
        printf("// Time elapsed for integration: %f\n", ((double)clock() - start2) / CLOCKS_PER_SEC);
///-----------------------------------///



///--------------------------------------Calculate Derivative Terms---------------------------------------///


        ///dA components

        for(int A=0; A<g; A++)
        {
            dAX[A]=0;
            dAY[A]=0;
            dAZ[A]=0;
        }



        for(int B=0; B<g; B++)
        {

            for(int k=0; k<g; k++)
            {
                dAX[B] +=   (f[k][B]) * RotatedX.GetKernel(B,k);
                dAY[B] +=   (f[k][B]) * RotatedY.GetKernel(B,k);
                dAZ[B] +=   (f[k][B]) * RotatedZ.GetKernel(B,k);
            }

        }


        for(int A=0; A<g; A++)
        {
            for(int B=0; B<g; B++)
            {
                dA[A][B] = - Point[B].GetX()*sin( RotatedX.GetThetaPrime(A) )*dAX[A]
                           - Point[B].GetY()*sin( RotatedY.GetThetaPrime(A) )*dAY[A]
                           - Point[B].GetZ()*sin( RotatedZ.GetThetaPrime(A) )*dAZ[A];
            }
        }


        ///dB components

        for(int A=0; A<g; A++)
        {
            dBX[A]=0;
            dBY[A]=0;
            dBZ[A]=0;
        }

        for(int A=0; A<g; A++)
        {

            for(int k=0; k<g; k++)
            {

                dBX[A] += (f[A][k]) * RotatedX.GetKernel(A,k); //[i][k]
                dBY[A] += (f[A][k]) * RotatedY.GetKernel(A,k);
                dBZ[A] += (f[A][k]) * RotatedZ.GetKernel(A,k);
            }


        }


        for(int A=0; A<g; A++)
        {
            for(int B=0; B<g; B++)
            {
                dB[A][B] = - Point[A].GetX()*sin( RotatedX.GetThetaPrime(B) )*dBX[B]
                           - Point[A].GetY()*sin( RotatedY.GetThetaPrime(B) )*dBY[B]
                           - Point[A].GetZ()*sin( RotatedZ.GetThetaPrime(B) )*dBZ[B];

            }
        }

///-----------------------------Hubble constant / scaling factor------------------------------------------///


        double Hubble=0; //Hubble rate
        //Hubble=0.001; //cosmological constant dominated
        //Hubble= 1.0/(2.0*timer); //radiation dominated
        //Hubble= 2.0/(3.0*timer); //matter dominated
        Hubble=(2.0)/(3*timer);

///-----------------------------------------------------------------------------------------------///






///-----------------------Update Gravitational Terms ----------------------------------------------------///

        std::vector< std::vector<double> > ddf(g);
        std::vector<double> fA(g), fB(g);



        for(int i=0; i<g; i++)
        {
            for(int j=0; j<g; j++)
            {
                dfg[i][j] = (dt)*Hubble*( dA[i][j]  + dB[i][j] -1.0 - 5.0*(Point[i]<Point[j]) )*f[i][j];
            }
        }



///-----------------------------------------Update Distribution--------------------------------------///


        large_loops = 0;
        small_loops =0;
        redshift=0;

        for(int A=0; A<g; A++)
        {
            for(int B=0; B<g; B++)
            {

                // f[A][B] = f[A][B] + df[A][B]; //no gravitational term or loop removal, must disable loop removal in integral term to use this

                //f[A][B] = f[A][B] + df[A][B] + dfg[A][B]; //update with gravitational term, must disable loop removal in integral term to use this

                f[A][B] = (f[A][B]  + df_small[A][B] + df_large[A][B] + dfg[A][B]); //remove small and large loops, adding in df_long should give same results as no loop removal

                large_loops += df_large[A][B]*w[A]*w[B];
                small_loops += df_small[A][B]*w[A]*w[B];
                redshift += dfg[A][B]*w[A]*w[B];

            }
        }

        double drho = small_loops + large_loops + redshift;


        double x=0, y=0, z=0;

        x=-1; //x<0 //small loops
        y=0.5;  //y>0 //large loops
        z=-1;  //z<0 //redshifted loops

        double gamma=Hubble*timer;

       // double fun= (-1.0)*((1+2*y + gamma*(1+2*z))/(1+2*x) );//(1+2*y)*(1+z)*(1+2*gamma);
        double fun= (1+2*y + 0.57*gamma*(z-y))/(-1-2*x+0.57*gamma*(x-z)) ;//(1+2*y)*(1+z)*(1+2*gamma);

        std::cout << "x: "<< x <<", |x|-y: "<< std::abs(x)-y <<", y+z: " << y+z<< std::endl;

        std::cout << "Predicted: " << fun << std::endl;

        //double AA = ( 1+2*y + (0.912786-0.584526*gamma)*gamma*(-y + z ) )/(2*y - 2*x);
        //double BB = ( 1+2*x + (0.912786-0.584526*gamma)*gamma*(-x + z ) )/(2*x - 2*y);
        double AA = ( 1+2*y + ((1 + z + (1.0/3.0)*x - (2.0/3.0)*y)/(-gamma+z))*gamma*(-y + z ) )/(2*y - 2*x);
        double BB = ( 1+2*x + ((1 + z + (1.0/3.0)*x - (2.0/3.0)*y)/(-gamma+z))*gamma*(-x + z ) )/(2*x - 2*y);
        //double CC = (0.467006)*gamma + (0.00051556)*z + (0.0385177)*gamma*z + (-0.255878)*gamma*gamma + (0.00498906)*z*z;
        double CC =1-AA-BB;

        //std::cout << "ll:" << large_loops << ", sm: " << small_loops << ", rs: " << redshift << std::endl;
        std::cout << "A: " << small_loops/drho << ", Predicted A:" << AA << std::endl;
        std::cout << "B: " << large_loops/drho << ", Predicted B:" << BB << std::endl;
        std::cout << "C: " << redshift/drho << ", Predicted C:" << CC << std::endl;
        std::cout << "A+B+C = " << AA + BB + CC << std::endl;

        cl = cl * (1+x*small_loops/rho0 + y*large_loops/rho0 + z*redshift/rho0);

        std::cout << "A/B= " << AA/BB << std::endl;

        std::cout << "correlation length / interstring distance: "<< cl*sqrt(rho0) << std::endl;

        std::cout << "x: " << x << std::endl;
        std::cout << "y: " << y << std::endl;
        std::cout << "z: " << z << std::endl;
        std::cout << "gamma: " << gamma << std::endl;




        //Energy density
        energy=0;

        for(int A=0; A<g; A++)
        {
            for(int B=0; B<g; B++)
            {
                energy += f[A][B] * w[A]*w[B];
            }
        }




        //calculate Entropy Limit

        for(int A=0; A<g; A++)
        {
            for(int B=0; B<g; B++)
            {
                fA[A] += f[A][B]*w[B]/energy;
            }
        }

        for(int A=0; A<g; A++)
        {
            for(int B=0; B<g; B++)
            {
                fB[A] += f[B][A]*w[B]/energy;
            }
        }

        for(int A=0; A<g; A++)
        {
            for(int B=0; B<g; B++)
            {
                EntropyLimit[s] += energy*fA[A]*fB[B]*log(energy*fA[A]*fB[B])*w[A]*w[B];
            }
        }

        //Entropy
        double S=0;

        for(int A=0; A<g; A++)
        {
            for(int B=0; B<g; B++)
            {
                S += f[A][B] * log(f[A][B]) * w[A]*w[B];
            }
        }




        Distribution::init_dist = f;

        Entropy[s] = S;
        Time[s] = timer;
        Energy[s] = energy;



///calculate slope of log(rho) vs log(time)

        double slope=0;
        double L=1.1;


        if(s>L)
        {
            slope = ( log(Energy[floor((double)s/L)]) - log(Energy[s]) )/( log(Time[floor((double)s/L)]) - log(Time[s]) );
        }



        std::cout << "Slope: " << slope << std::endl;


        avg_slope=0;
        Slope[s]=slope;


        for(int i=0; i<s; i++)
        {
            avg_slope += (double)Slope[i]/((double)(s));
        }






        ///----------------Get P(A) and P(B) ----------------///

//symmetry check for debugging

        for(int A=0; A<g; A++)
        {
            pA[A]=0;
            pB[A]=0;

            for(int B=0; B<g; B++)
            {
                pA[A] += f[A][B]*w[B]*w[A]/energy;
                pB[A] += f[B][A]*w[B]*w[A]/energy;

                myRowsum[A] += f[A][B]*w[B];
                myColsum[B] += f[B][A]*w[B];
            }
        }

        P_A=0;
        P_B=0;

        for(int i=0; i<g; i++)
        {
            P_A += pA[i];
            P_B += pB[i];
        }

        std::cout << P_A << ", " << P_B << std::endl;


        ///-------------------------------------------------------///


        ///-------------Get <U>^2, <V>^2, <U^2>, <V^2> -----------///


        Ubar.Setup(0,0,0);
        Vbar.Setup(0,0,0);


        U2=0;
        V2=0;


        Vex Vi;
        Vi=V;


        U.Setup(0,0,0);
        V.Setup(0,0,0);

        for(int A=0; A<g; A++)
        {
            for(int B=0; B<g; B++)
            {

                U = U + ((Point[B]-Point[A])/2.0)*f[A][B]*w[A]*w[B]/energy;
                V = V + ((Point[A]+Point[B])/2.0)*f[A][B]*w[A]*w[B]/energy;

                U2 += (((Point[B]-Point[A])/2.0)<((Point[B]-Point[A])/2.0))*f[A][B]*w[A]*w[B]/energy;
                V2 += (((Point[A]+Point[B])/2.0)<((Point[A]+Point[B])/2.0))*f[A][B]*w[A]*w[B]/energy;

            }
        }

        Falloff[s]=V2;







///-----------Calculate energy-momentum Tensory, pressure and equation of state w-----------///

        Vex zero;
        zero.Setup(0,0,0);

        Abar.Setup(0,0,0);
        Bbar.Setup(0,0,0);

        Abar=zero;
        Bbar=zero;

        double normA=0;
        double normB=0;

        for(int i=0; i<g; i++)
        {

            normA = normA + myRowsum[i]*w[i];
            normB = normB + myColsum[i]*w[i];

            Abar = Abar + Point[i]*myRowsum[i]*w[i];
            Bbar = Bbar + Point[i]*myColsum[i]*w[i];


        }

        Abar = Abar*(1.0/normA);
        Bbar = Bbar*(1.0/normB);





        std::vector<double> a1(3), b1(3);
        std::vector< std::vector<double> > AB(g);

        for(int i=0; i<g; i++)
        {
            AB[i].resize(g);
        }


        a1[0]=Abar.GetX();
        a1[1]=Abar.GetY();
        a1[2]=Abar.GetZ();

        b1[0]=Bbar.GetX();
        b1[1]=Bbar.GetY();
        b1[2]=Bbar.GetZ();




        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                AB[i][j]=0;

                for(int A=0; A<g; A++)
                {
                    for(int B=0; B<g; B++)
                    {
                        AB[i][j] = (a1[i]*b1[j]*f[A][B]*w[A]*w[B])/energy;
                    }
                }
            }
        }

        double omega=0;

        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                if(i==j)
                {
                    omega += (1.0/3.0)*AB[i][j]/(energy);
                }

            }
        }

        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "// energy momentum tensor:" << std::endl;
        std::cout << std::endl;

        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                std::cout << AB[i][j] << ", ";
            }
            std::cout << std::endl;
        }

        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "// equation of state w:" << omega << std::endl;
        std::cout << std::endl;


///-------------------------------------------------------///

///-----------------------------DONE WITH CALCULATIONS------------------------------------///





///----------------------------Build Output Files---------------------------------------///

        std::ofstream myfile4;
        std::ofstream myfile5;
        std::ofstream myfile6;

        myfile4.open("entropy.txt");
        myfile6.open("matrix");

        double r=0;


        for(int i=0; i<g; i++)
        {
            for(int j=0; j<g; j++)
            {
                myfile6 << f[i][j] << " ";
            }
            myfile6 << std::endl;
        }

        for (int i=1; i<100; i++)
        {

            r=i*std::abs(ceil(s/100));

            myfile4 << (Time[r]) << "   " <<  (Entropy[r]) << "    " << log(Energy[r]) << "    " << (EntropyLimit[r]) << "     " << Falloff[r] << "       " << Slope[r] << "       " << (log(Time[r])) << std::endl;


        }

        myfile4.close();


///-------------------------------------------------------------------------------///



///--------------------Do PLOTS-----------------------------------------------///

        if(sphere==true) //PLOT SPHERE
        {


            //g3.cmd("set palette defined ( 0 'black', 5 'red', 10 'green')");
            g3.cmd("set terminal x11");
            // g3.cmd("set size 1,1");
            // g3.cmd("set terminal x11 size 400,400");

            //g3.cmd("set terminal wxt size 600,600");
            g3.set_grid().set_xrange(-0.04,0.04).set_yrange(-0.04,0.04).set_zrange(-0.04,0.04);
            //g1.cmd("set pm3d depthorder hidden3d 2");
            //g3.cmd("set palette rgbformulae 8, 9, 7");


            //g1.cmd("set style line 2  linetype 1 linecolor palette linewidth 0.01"); // rgb \"#a0a0f0\"
            // g1.cmd("set style fill  transparent solid 0.90 border");

            g3.cmd("set style line 2 linetype 1 linecolor palette linewidth 1");
            g3.set_grid();

            //g1.cmd("set pm3d map");
            g3.cmd("set ticslevel 0");
            //g1.cmd("set arrow to $1");
            //g3.set_grid().set_xautoscale().set_yautoscale().set_zautoscale();
            // g3.cmd("splot 'dataA.txt' using 1:2:3:4 with points pt 5 ps 1 lc palette, \ 'dataB.txt' using 1:2:3:4 with points pt 7 ps 1 lc palette");
            g3.cmd("splot 'dataA.txt' using 1:2:3:4 with points pt 3 ps 1 lc palette");
            //g3.set_grid();
            //g3.reset_plot();

            //g1.cmd("");
            //g3.cmd("splot 'dataA.txt' u 1:2:3:4:5:6:7 with vectors lc palette lw 3, \ 'dataB.txt' u 1:2:3:4 with points pt 7 ps variable lc palette");
            g3.remove_tmpfiles();
            //std::system("clear");

        }
        if(ent==true)
        {





            g4.cmd("set terminal x11");
            g4.set_grid();
            g4.cmd("set size 1,1");
            g4.cmd("set terminal x11 size 400,400");

            //g4.set_xrange(0,s);
            //g4.set_yautoscale();
            //g4.set_yrange(-0.001, -0.002);
            //g4.set_xrange(0, t);
            g4.cmd("set ylabel \"H(t)\"");
            g4.cmd("set title \"Entropy\"");
            g4.cmd("plot 'entropy.txt' using 1:2 with lines lc rgb 'green', \ 'entropy.txt' using 1:4 with lines lc rgb 'red'"); //'entropy.txt' using 1:7 with points pt 6 lw 1 lc rgb 'blue'");
            std::system("clear");

            g8.cmd("set terminal x11");
            g8.set_grid();
            g8.cmd("set size 1,1");
            g8.cmd("set terminal x11 size 400,400");
            //g8.set_yrange(-0.3,0);
            //g6.set_xrange(0,(dt*t));
            g8.cmd("set ylabel \"E(t)\"");
            g8.cmd("set title \"Energy\"");
            //g6.cmd("plot 'entropy.txt' using 1:3 with points lt 6 lw 3 ps 2 lc rgb  'red', \ 'entropy.txt' using 1:6 with points lt 6 lw 3 lc rgb  'blue'");
            g8.cmd("plot 'entropy.txt' using 7:3 with lines lt 6 lw 1 lc rgb 'blue'");//, \ 'entropy.txt' using 7:5 with lines lc rgb 'purple'");
            std::system("clear");





            g6.cmd("set terminal x11");
            g6.set_grid();
            g6.cmd("set size 1,1");
            g6.cmd("set terminal x11 size 400,400");
            //g4.set_yrange(Energy[1], Energy[s]);
            //g6.set_xrange(0,(dt*t));
            g6.cmd("set ylabel \" <V^2> \"");
            g6.cmd("set title \" <V^2> vs t \"");
            g6.cmd("plot 'entropy.txt' using 7:5 with lines lt 6 lw 1 lc rgb  'red'");//", \ 'entropy.txt' using 1:6 with points lt 6 lw 3 lc rgb  'blue'");
            //g6.cmd("plot 'entropy.txt' using 1:10 with lines lt 6 lw 1 lc rgb  'black', \ 'entropy.txt' using 1:11 with lines lt 6 lw 1 lc rgb 'green' ");


            g7.cmd("set terminal x11");
            g7.set_grid();

            g7.set_xlabel("j");
            g7.set_ylabel("i");
            //g7.set_xrange(-20, size+20).set_yrange(-20,size+20);
            g7.cmd("set pm3d map");
            g7.set_cbrange(0,1);
            //g7.cmd("set palette model RGB defined ( 0 'white', 1 'black' ) ");
            g7.cmd("splot 'matrix' matrix with image lc palette");



            /*




                                        g7.cmd("set terminal x11");
                                        g7.set_grid();
                                        g7.cmd("set size 1,1");
                                        //g7.cmd("set terminal x11 size 400,400");
                                        // g7.set_yrange(-0.0000005, 0.0000005);
                                        //g7.set_yrange(0.261, 0.263);
                                        //g7.set_yrange(V_list[1].GetR(), V_list[s].GetR() );
                                        g7.set_xrange(0,t);
                                        g7.cmd("set ylabel \"V(t)\"");
                                        g7.cmd("set title \"Actual Average Velocity\"");
                                        g7.cmd("plot 'entropy.txt' using 1:4 with points lt 6 lw 1 lc rgb  'red'");


                        g8.cmd("set terminal x11");
                        g8.set_grid();
                        g8.cmd("set size 1,1");
                        //g8.set_yrange(dVTheory[2].GetR(),dVTheory[s].GetR());
                        //g8.cmd("set terminal x11 size 400,400");
                        //g8.set_yautoscale();
                        //g8.set_yrange(-0.000000005, 0.000000005);
                        g8.set_yrange(0, 0.01);
                        g8.set_xrange(0,(dt*t));
                        g8.cmd("set ylabel \"V(t)\"");
                        g8.cmd("set title \"Theoretical Average Velocity\"");
                        //g8.cmd("")
                        g8.cmd("plot 'entropy.txt' using 1:4 with points lt 6 ps 1 lc rgb  'red', \ 'entropy.txt' using 1:5 with points lt 6 ps 2 lc rgb  'green'");
                        // g8.cmd("plot 'entropy.txt' using 1:5 with points lt 6 lw 1 lc rgb  'green'");
                        */



            //g4.set_yrange(400, H_max);
            //g4.set_cbrange(0,time);
            //g4.cmd("set")
            //g4.set_grid();
            // g4.plot_xy(Time, Entropy, "H(t)");


            /*
                        std::string Result2;          // string which will contain the result
                        std::ostringstream convert2;   // stream used for the conversion
                        convert2 << std::setfill('0') << std::setw(2) << s;      // insert the textual representation of 'Number' in the characters in the stream
                        Result2 = convert2.str();

                        std::string mystring2;
                        mystring2 = "i = " + Result2;

                        std::cout << "mystring2= " << mystring2 << std::endl;

                        g4.cmd(mystring2);
                        g4.cmd("n=1");
                        g4.cmd("load \"loop2.plt\"");
                        */

        }


        std::system("clear");

    } //end time loop




///------------end timer-------------------------///
    printf("Time elapsed for full time loop: %f\n", ((double)clock() - start1) / CLOCKS_PER_SEC);
///-----------------------------------///



}


Distribution::~Distribution()
{
    //dtor
}



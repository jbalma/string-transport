//calculate associate legendre polynomial
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#define PI 3.1415926535897932384626433832795028841971693993751058



class LEGENDRE
{
public:

    double factorial(int n) //a factorial function for use in normalization function
    {

        int i=0;
        double fact=1.0;

        if(n<=1)
        {
            return(1);
        }

        else
        {
            for(i=1; i<=n; i++)
            {
                fact=fact*i;
            }
            //std::cout << fact << std::endl;
            return(fact);
        }


    }




    double Normalization(int m, int l) //calculates the normilzation factor for Y
    {
        double norm=1;


        norm = sqrt( ( (2.0*l+1.0) / (4.0*PI) ) * (factorial(l-m) / factorial(l+m)) );

        //std::cout << norm;


       // std::cout << norm << std::endl;


        return norm;
    }

    std::complex<double> Normalization_derivative(int mm, int ll) //calculates the normilzation factor for dY
    {
        std::complex<double> norm_c(0,0);

        std::complex<double> m, l;

        m=mm;
        l=ll;


        norm_c = sqrt(  (2.0*ll+1.0) / (4.0*PI) ) * sqrt( std::tgamma(ll-mm+1) ) / sqrt(std::tgamma(ll+mm+1));
        //std::cout << norm;


        //std::cout << norm << std::endl;


        return norm_c;
    }







    double P(int m, int l, double x)
    {

//takes a function of m, l and x, and returns a new function P(x)

        int i, ll;
        double fact, pll,pmm,pmmp1,somx2;





        if (m<0 || m>l || std::abs(x) > 1.0)
        {
            //std::cout << "Bad arguments";
            return 0;
        }
        pmm=1.0;

        if(m>0) //compute p(m,m) --->implies m=1.0
        {
            somx2=sqrt((1.0-x)*(1.0+x));
            fact=1.0;

            for(i=1; i<=m; i++)
            {
                pmm *= -fact*somx2;
                fact += 2.0;
            }
        }



        if (l==m) //compute p(m, m+1)
            return pmm;
        else
        {
            pmmp1=x*(2*m+1)*pmm;
            if(l==(m+1))
                return pmmp1;
            else  						//compute p(m,l > m+1)
            {
                for(ll=m+2; ll<=l; ll++)
                {
                    pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
                    pmm=pmmp1;
                    pmmp1=pll;

                }
                return pll;
            }
        }



    }

   double Pml(int m, int l, double x)
    {
        std::complex<double> P;

        P = boost::math::legendre_p(l, m, x);

        return P.real();

    }



    std::complex<double> Yml(int m, int l, double theta, double phi)
    {
        std::complex<double> Y;

        Y= boost::math::spherical_harmonic(l, m, theta, phi);

        return Y;
    }




};

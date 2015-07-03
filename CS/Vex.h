#ifndef VEX_H
#define VEX_H


class Vex
{
     private:

        double x, y, z;
        double phi, theta, r;



    public:

    //double TOL=1.0e-10;

        Vex();

        Vex operator^(Vex); //cross product

        double operator<(Vex); //dot product

        Vex operator*(double); //scalar mult

        Vex operator/(double); //scalar divide

        Vex operator-(Vex); //subtract vectors

        Vex operator+(Vex); //add vectors

        Vex operator-(double); //subtract vectors

        Vex operator+(double); //add vectors

        bool operator==(Vex);

        void Setup(double, double, double);

        void SetXYZ(double, double, double);

        void SetSpherical(double, double, double);

        void XYZtoSpherical(double, double, double);

        void SphericaltoXYZ(double, double, double);

        Vex GetVector();

        double GetX();
        double GetY();
        double GetZ();

        double GetTheta();
        double GetPhi();
        double GetR();

        int GetIt(int);


        virtual ~Vex();

    protected:



};

#endif // VEX_H

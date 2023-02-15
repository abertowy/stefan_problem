#include <vector>

class StefanProblemSolver {
public:
    StefanProblemSolver(){}

    void calculate(unsigned int val);
    void reinitPhi(std::vector<std::vector<double> > &re_val, double dx);
    void clean();

    std::vector<double> PTB {0.0};
    std::vector <double> phi;

    long long solveTime;
    std::vector<double> ts {0.0};

private:
    int spy_{60*60*24*365};              //Seconds per year
    // (Van der Veen 2013)
    // pg 144 Van der Veen - from Yen CRREL 1981
    double Lf_ {3.335e5};                 //Latent heat of fusion (J/kg)
    double rho_ {1000};                   //Bulk density of water (kg/m3), density changes are ignored
    double Ks_ {spy_*2.1};                 //Conductivity of ice (J/mKs)
    double cs_ {2009};                    //Heat capacity of ice (J/kgK) - ** Van der Veen uses 2097 but see Tr and Aschwanden 2012)
    double ks_ {Ks_/(rho_*cs_)};             //Cold ice diffusivity (m2/sec)
    // Engineering Toolbox
    double Kl_ {spy_*0.58};                //Conductivity of water (J/mKs)
    double cl_ {4217};                    //Heat capacity of water (J/kgK)
    double kl_ {Kl_/(rho_*cl_)};

    // Problem Constants
    //std::vector<double> ptb_ {0.0};      //phase transition boundary
    double t0_ {0.0};                    //start time
    double Tm_ {0.0};                    //melting point
    double T_ {-10.0};                  //left boundary condition
    double Tr_ {-16.33};


    double k {0.01};
    double kd {0.03};
};



#ifndef _FITGLE_H_
#define _FITGLE_H_

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <algorithm>

namespace FITGLE_NS {

struct InputParameters
{
    double start;  //start distance r0
    double end;    //end distance r1
    int    splineOrder;
    int    numSplines;
    FILE*  fileTraj;
};  //Structure to store input parameters

class FitGLE
{
public:
    FitGLE(int argc, char** argv);
    ~FitGLE();
    void exec();

    //helper functions
    void accumulateNormalEquation();
    void leastSquareSolver();
    void output();

private:
    Frame* trajFrame;
    struct InputParameters* info; 
    std::vector<std::vector<double> > normalMatrix;
    std::vector<double> normalVector;
    std::vector<double> splineCoefficients;
};

}

#endif

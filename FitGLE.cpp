#include <FitGLE.h>
#include <comm.h>
#include <cassert>
#include <gsl/gsl_bspline.h>

using namespace FITGLE_NS;

FitGLE::FitGLE(int argc, char** argv)
{
    assert(argc == 4);
    printf("Initializing FitGLE parameters...\n");

    // parsing the configuration parameters
    VAR_BEGIN
      GET_REAL(info->start)
      GET_REAL(info->end)
      GET_REAL(info->boxLength)
      GET_INT(info->splineOrder)
      GET_INT(info->numSplines)
      GET_INT(info->steps)
    VAR_END
           
    // Initialize the Normal Equation matrix and vectors
    // Set up the size of splines according to order and numbers
    int numBreakds = info->numSplines + 2 - info->splineOrder;
    normalVector.resize(numSplines);
    splineCoefficients.resize(numSplines);
    normalMatrix.resize(numSplines);
    for (auto&& i : normalMatrix)
    {
        i.resize(numSplines);
    }
}

FitGLE::~FitGLE()
{
    printf("Exiting the Fitting GLE process...\n");
}


// Accumulate the normal equation for this particular frame
void FitGLE::accumulateNormalEquation()
{
    int nall = trajFrame->numParticles;
    int nSplines = info->numSplines;

    std::vector<std::vector<double> > frameMatrix(nall, std::vector<double>(nSplines));
    
    for (int i=0; i<nall; i++)
    {
        for (int j=0; j<nall; j++)
        {
            double rij = distance(frameMatrix->positions[i], frameMatrix->positions[j])
        }
    }
    
}

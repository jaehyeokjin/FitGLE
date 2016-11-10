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
    int numBreaks = info->numSplines + 2 - info->splineOrder;
    normalVector.resize(numSplines);
    splineCoefficients.resize(numSplines);
    normalMatrix.resize(numSplines);
    for (auto&& i : normalMatrix)
    {
        i.resize(numSplines);
    }

    // Initialize the spline set up
    bw = gsl_bspline_alloc(info->splineOrder, numBreaks);
    splineValue = gsl_vector_alloc(info->numSplines);
    gsl_bspline_knots_uniform(info->start, info->end, bw);
}

FitGLE::~FitGLE()
{
    gsl_bspline_free(bw);
    gsl_vector_free(B);
    printf("Exiting the Fitting GLE process...\n");
}


// Accumulate the normal equation for this particular frame
void FitGLE::accumulateNormalEquation()
{
    int nall = trajFrame->numParticles;
    int nSplines = info->numSplines;

    std::vector<std::vector<double> > frameMatrix(nall, std::vector<double>(nSplines));
   
    // Computing Matrix F_km 
    for (int i=0; i<nall-1 ; i++)
    {
        for (int j = i + 1; j<nall; j++)
        {
            double rij = distance(trajFrame->positions[i], trajFrame->positions[j]);
            gsl_bspline_eval(rij, splineValue, bw);
            
            double dv = parallelVelocity(trajFrame->velocities[i], trajFrame->velocities[j]);
            
            for (int m=0; m<nSplines; m++)
            {
               double phim = gsl_vector_get(splineValue, m);
               frameMatrix[i][m] += phim * dv;
               frameMatrix[j][m] -= phim * dv;
            }
        }
    }
 
    // Constructing the normal Matrix and normal Vector
    for (int m=0; m<nSplines; m++)
    {
        for (int n=0; n<nSplines; n++)
        {
            double sum = 0.0;
            for (int k=0; k<nall; k++)
                sum += frameMatrix[k][m] * frameMatrix[k][n];
            normalMatrix[m][n] = sum;
        }
   
        double sum_b = 0.0; 
        for (int k=0; k<nall; k++)
            sum_b += frameMatrix[k][m] *   
    }

    
}

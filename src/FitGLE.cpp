#include "FitGLE.h"
#include "comm.h"
#include <cmath>
#include <numeric>
#include <functional>
#include <algorithm>
#include <cassert>
#include <gsl/gsl_bspline.h>
#include <lapacke.h>

using namespace FITGLE_NS;

FitGLE::FitGLE(int argc, char** argv)
{
    if (argc != 3)
    {
        printf("./FitGLE.x [trajectory Filename] [number of Particles]\n");
    }

    assert(argc == 3);
    printf("Initializing FitGLE parameters...\n");

    // parsing the configuration parameters
    info = std::make_shared<InputParameters>();
    VAR_BEGIN
      GET_REAL(info->start)
      GET_REAL(info->end)
      GET_REAL(info->boxLength)
      GET_REAL(info->outputPrecision)
      GET_INT(info->splineOrder)
      GET_INT(info->numSplines)
      GET_INT(info->steps)
    VAR_END
           
    printf("set up trajectory files\n");
    trajFrame = std::make_shared<Frame>(atoi(argv[2]), argv[1]);
    //trajFrame = new Frame(atoi(argv[2]), argv[1]);
    printf("before: %u\n", trajFrame);
    // Initialize the Normal Equation matrix and vectors
    // Set up the size of splines according to order and numbers
    printf("set up b-spline data structures\n");
    int numBreaks = info->numSplines + 2 - info->splineOrder;
    normalVector.resize(info->numSplines);
    splineCoefficients.resize(info->numSplines);
    normalMatrix.resize(info->numSplines);
    printf("set up containers\n");
    for (auto&& i : normalMatrix)
    {
        i.resize(info->numSplines);
    }

    // Initialize the spline set up
    bw = gsl_bspline_alloc(info->splineOrder, numBreaks);
    splineValue = gsl_vector_alloc(info->numSplines);
    gsl_bspline_knots_uniform(info->start, info->end, bw);
    printf("finishing configuration, entering normal equation accumulation\n");
}

FitGLE::~FitGLE()
{
    gsl_bspline_free(bw);
    gsl_vector_free(splineValue);
    printf("Exiting the Fitting GLE process...\n");
}

inline double FitGLE::distance(std::vector<double> & A, std::vector<double> & B)
{
    double dx = A[0] - B[0];
    if (dx > 0.5 * info->boxLength) dx = dx - info->boxLength;
    if (dx < -0.5 * info->boxLength) dx = dx + info->boxLength;
    double dy = A[1] - B[1];
    if (dy > 0.5 * info->boxLength) dy = dy - info->boxLength;
    if (dy < -0.5 * info->boxLength) dy = dy + info->boxLength;
    double dz = A[2] - B[2];
    if (dz > 0.5 * info->boxLength) dz = dz - info->boxLength;
    if (dz < -0.5 * info->boxLength) dz = dz + info->boxLength;
  
    return sqrt(dx*dx + dy*dy + dz*dz);
}

inline std::vector<double> FitGLE::parallelVelocity(int i, int j)
{
    double dx = trajFrame->positions[i][0] - trajFrame->positions[j][0];
    if (dx > 0.5 * info->boxLength) dx = dx - info->boxLength;
    if (dx < -0.5 * info->boxLength) dx = dx + info->boxLength;
    double dy = trajFrame->positions[i][1] - trajFrame->positions[j][1];
    if (dy > 0.5 * info->boxLength) dy = dy - info->boxLength;
    if (dy < -0.5 * info->boxLength) dy = dy + info->boxLength;
    double dz = trajFrame->positions[i][2] - trajFrame->positions[j][2];
    if (dz > 0.5 * info->boxLength) dz = dz - info->boxLength;
    if (dz < -0.5 * info->boxLength) dz = dz + info->boxLength;

    double rij = sqrt(dx*dx + dy*dy + dz*dz);
    double eij[] = {dx/rij, dy/rij, dz/rij};
    std::vector<double> vij;
    std::transform(trajFrame->velocities[i].begin(), trajFrame->velocities[i].end(), trajFrame->velocities[j].begin(), std::back_inserter(vij), std::minus<double>());
     
    double projection = vij[0] * eij[0] + vij[1] * eij[1] + vij[2] * eij[2];
    vij[0] = projection * eij[0];
    vij[1] = projection * eij[1];
    vij[2] = projection * eij[2];
    return vij;
}

// Accumulate the normal equation for this particular frame
void FitGLE::accumulateNormalEquation()
{
    int nall = trajFrame->numParticles;
    printf("np = %d\n", nall);
    int nSplines = info->numSplines;
    printf("nall %d splines %d\n", nall, nSplines);
    std::vector<std::vector<double> > frameMatrix(3 * nall, std::vector<double>(nSplines));
   
    // Computing Matrix F_km 
    for (int i=0; i<nall-1 ; i++)
    {
        for (int j = i + 1; j<nall; j++)
        {
            double rij = distance(trajFrame->positions[i], trajFrame->positions[j]);
            printf("rij = %lf, %d %d\n", rij, i, j);    
            if (rij < info->end)
            {
                gsl_bspline_eval(rij, splineValue, bw);
                printf("rij = %lf, %d %d\n", rij, i, j);    
                std::vector<double> dv = parallelVelocity(i, j);
            
                for (int m=0; m<nSplines; m++)
                {
                     double phim = gsl_vector_get(splineValue, m);
                     // For all three dimensions
                     frameMatrix[3*i][m] += phim * dv[0];
                     frameMatrix[3*i + 1][m] += phim * dv[1];
                     frameMatrix[3*i + 2][m] += phim * dv[2];
                     frameMatrix[3*j][m] -= phim * dv[0];
                     frameMatrix[3*j + 1][m] -= phim * dv[1];
                     frameMatrix[3*j + 2][m] -= phim * dv[2];
                }
            }  
        }
    }
    printf("finishing F\n");
 
    // Constructing the normal Matrix and normal Vector
    for (int m=0; m<nSplines; m++)
    {
        for (int n=0; n<nSplines; n++)
        {
            double sum = 0.0;
            for (int k=0; k<3 * nall; k++)
                sum += frameMatrix[k][m] * frameMatrix[k][n];
            normalMatrix[m][n] = sum;
        }
   
        double sum_b = 0.0; 
        for (int k=0; k<3 * nall; k++)
            sum_b += frameMatrix[k][m] * trajFrame->residualForces[k/3][k%3];
        normalVector[m] = sum_b;
    }
}

void FitGLE::leastSquareSolver()
{
    // Solving the least square normal equation G*phi = b
    double* G = new double[info->numSplines * info->numSplines];
    double* b = new double[info->numSplines];

    for (int m=0; m<info->numSplines; m++)
    {
        for (int n=0; n<info->numSplines; n++)
        {
            G[m*info->numSplines + n] = normalMatrix[m][n];
        }
        b[m] = normalVector[m];
    }
    
    // Solving this linear equation using Lapacke
    int N = info->numSplines;
    int NRHS = 1;
    int LDA = N;
    int ldb = 1;
    int* ipiv = new int[N];
    int solverInfo = LAPACKE_dgesv(LAPACK_ROW_MAJOR, N, NRHS, G, LDA, ipiv, b, ldb);
   
    printf("LSQ Solver Info: %d\n", solverInfo);

    for (int m=0; m<info->numSplines; m++)
        splineCoefficients[m] = b[m];

    delete[] G;
    delete[] b;
    delete[] ipiv;
}

// Output function for gamma(R)
void FitGLE::output()
{
    double start = info->start;
    double end = info->end;
    double precision = info->outputPrecision;

    FILE* fp = fopen("gamma_out.dat", "w");

    while (start < end)
    {
        double gamma_r = 0.0;
        gsl_bspline_eval(start, splineValue, bw);
        for (int m=0; m<info->numSplines; m++)
        {
           gamma_r += splineCoefficients[m] * gsl_vector_get(splineValue, m);
        }
        fprintf(fp, "%lf\t%lf\n", start, gamma_r);
        start = start + precision;
    }
    
    fclose(fp);
}

// Execution Process
void FitGLE::exec()
{
    printf("Accumulating the LSQ normal Matrix\n");
    printf("np = %d\n", trajFrame->numParticles);
    printf("after: %u\n", trajFrame);
    for (int i=0; i<info->steps; i++)
    {
        trajFrame->readFrame();
        printf("%u\n", trajFrame);
        printf("reading finished, np = %d\n", trajFrame->numParticles);
        accumulateNormalEquation();
        printf("finishing step %d (total %d)\r", i+1, info->steps);
    }
    leastSquareSolver();
    output();
}

#ifndef _FRAME_H_
#define _FRAME_H_

#include <cstdio>
#include <cstdlib>
#include <vector>

namespace FITGLE_NS {

class Frame
{
public:
    Frame(int n, char* filename);
    ~Frame();
    void readFrame();
private:
    FILE* trajectory;
    int numParticles;
    std::vector<std::vector<double> > positions;
    std::vector<std::vector<double> > residualForces;
    std::vector<std::vector<double> > velocities;
};
}

#endif

#include <cstdio>
#include <cstdlib>
#include <vector>

namespace FITGLE_NS {

class Frame
{
public:
    Frame();
    ~Frame();
    void readFrame();
private:
    FILE* trajectory;
    int numParticles;
    std::vector<double> positions;
    std::vector<double> residualForces;
    std::vector<double> velocities;
};
}

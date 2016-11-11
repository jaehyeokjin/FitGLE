#include <FitGLE.h>
#include <Frame.h>
#include <cstdlib>
#include <cstdio>

int main(int argc, char** argv)
{
   FitGLE* fg = new FitGLE(argc, argv);
   fg->exec();
   delete fg;
   return 0;
}

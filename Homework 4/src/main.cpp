#include "GlobalPlacer.h"
#include "Timer.hpp"
#include "Wrapper.hpp"
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        cerr << "Usage: " << argv[0] << " <.aux file> <.gp.pl file>\n";
        return 0;
    }

    // ##################################################
    //                    Parse Input
    // ##################################################
    cout << "##################### Parse Input ######################\n";
    wrapper::Placement placement;
    placement.readBookshelfFormat(argv[1], "");

    cout << "Benchmark: " << placement.name() << "\n";
    printf("Core region: (%.f, %.f) (%.f, %.f)\n",
           placement.boundryLeft(), placement.boundryBottom(),
           placement.boundryRight(), placement.boundryTop());
    printf("HPWL: %.0f\n", placement.computeHpwl());

    Timer timer(0);
    // ##################################################
    //                  Global Placement
    // ##################################################
    cout << "\n################### Global Placement ###################\n";
    timer.startTimer("global");

    // TODO: implement your own placer!

    GlobalPlacer globalPlacer(placement);
    globalPlacer.SAplace(99999);
    //globalPlacer.place();   // random
    //globalPlacer.place2();



    // TODO END
    timer.stopTimer("global");

    placement.outputBookshelfFormat(argv[2]);
    double gpWirelength = placement.computeHpwl();
    printf("\nHPWL: %.0f    Time: %.2f sec (%.2f min)\n",
           gpWirelength, timer.getTime("global"), timer.getTime("global") / 60.0);
    return 0;
}

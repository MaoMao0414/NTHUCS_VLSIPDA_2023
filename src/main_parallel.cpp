#include "GlobalPlacer.h"
#include "Timer.hpp"
#include <omp.h>
#include "Wrapper.hpp"
#include <iostream>
#include "bits/stdc++.h"

using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        cerr << "Usage: " << argv[0] << " <.aux file> <.gp.pl file>\n";
        return 0;
    }

    // thread num
    omp_set_num_threads(16);

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

    double bestHpwl = DBL_MAX;

    srand(0);

    #pragma omp parallel
    {

        wrapper::Placement localPlacement;
        localPlacement.readBookshelfFormat(argv[1], "");

        GlobalPlacer localGlobalPlacer(localPlacement);

        int threadId = omp_get_thread_num();
        int randSeed = threadId * (rand() % 100000);

        // 將其中一個thread設成99999(跟sequential version一樣以觀察)
        if(threadId == 0)
            randSeed = 99999;

        cout << "Thread " << threadId << " is running, rng seed value: " << randSeed << endl;

        localGlobalPlacer.SAplace(randSeed);  // SA

        double localHpwl = localPlacement.computeHpwl();  // 計算HPWL

        printf("Thread %d finishing SA, HPWL: %.f\n", threadId, localHpwl);

        // 更新最佳結果
        #pragma omp critical
        {
            if (localHpwl < bestHpwl) {

                cout << "Thread " << threadId << "updating best value!" << endl;

                bestHpwl = localHpwl;
                
                for (size_t i = 0; i < localPlacement.numModules(); ++i)
                {
                    double x = localPlacement.module(i).x();
                    double y = localPlacement.module(i).y();

                    placement.module(i).setPosition(x, y);
                }
            }
        }
    }

    // TODO END
    timer.stopTimer("global");

    placement.outputBookshelfFormat(argv[2]);
    double gpWirelength = placement.computeHpwl();
    printf("\nHPWL: %.0f    Time: %.2f sec (%.2f min)\n",
           gpWirelength, timer.getTime("global"), timer.getTime("global") / 60.0);
    return 0;
}

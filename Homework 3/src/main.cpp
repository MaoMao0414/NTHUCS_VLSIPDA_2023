#include "bits/stdc++.h"
#include "floorplan.hpp"

using namespace std;

int main(int argc, char** argv)
{
    start_time = chrono::high_resolution_clock::now();

    // 經典seed
    srand(1542955417);

    FloorPlan fp;

    fp.parseInputFile(argv[1]);

    int count = 0;
    do
    {
        //cout << "SA Round: " << count++ << endl;
        fp.SA(570000);
    }
    while(!fp.getTimeOut());


    if(!fp.getbestLegal())
        fp.refine();

    fp.writeOutputFile(argv[2]);

    auto end_time_parser = chrono::high_resolution_clock::now();
    auto duration_parser = chrono::duration_cast<chrono::milliseconds>(end_time_parser-start_time).count();

    //cout << "total run time(ms): " << duration_parser << endl;
    
    return 0;
}
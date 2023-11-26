#include "partition.hpp"
#include <chrono>
#include <bits/stdc++.h>

using namespace std;

int main(int argc, char** argv)
{
    start_time = chrono::high_resolution_clock::now();

    Partition myPartition;

    myPartition.parseInputFile(argv[1]);
    cout << "pass init" << endl;

    myPartition.initialPartition();
    myPartition.initialNetDistributionGain();
    myPartition.iterrativeRevise();

    cout << "pass" << endl;
    myPartition.writeOutputFile(argv[2]);


    return 0;
}
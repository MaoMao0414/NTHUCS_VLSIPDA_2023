#include <bits/stdc++.h>
#include <chrono>
#include "cell.h"
#include "cmpFunc.hpp"

#define ll long long

using namespace std;

std::chrono::time_point<std::chrono::high_resolution_clock> start_time;

class Partition
{
    public:
        void parseInputFile(string);
        bool cmpPartition(bool (*func)(Cell*, Cell*));
        void initialPartition();
        void initialNetDistributionGain();
        void iterrativeRevise();
        void writeOutputFile(string);
    private:
        ll techNum, netListcellNum, netNum;
        unordered_map<string, unordered_map<string, ll>> umap_cellNameToTechInfo;
        string dieTech[2];      // 0 是 dieA, 1是dieBf
        ll dieMaximumArea;
        ll dieUtilization[2];
        ll dieArea[2];
        unordered_map<string, string> umap_netlistCellToCellName;
        unordered_map<string, Cell*> umap_nameToNetlistCell;
        unordered_map<string, Net*> umap_nameToNet;
        vector<Cell*> vec_netlistCell;
        vector<Net*> vec_net;
        map<ll, unordered_set<Cell*>> map_bucketList;    // 維護兩個, 0 是 dieA, 1 是 dieB
};

void Partition::parseInputFile(string path)
{
    ifstream inFile;
    inFile.open(path);

    if(!inFile.is_open())
    {
        cerr << "Error in parseInputFile(), no input file" << endl;
        exit(1);
    }

    string temp;
    // first Read the number of technology file
    inFile >> temp >> techNum;

    for(int i = 0; i < techNum; i++)
    {   
        string curTechName;
        ll curTechCellCnt;

        inFile >> temp >> curTechName >> curTechCellCnt;
        for(int j = 0; j < curTechCellCnt; j++)
        {
            string cellName;
            ll width, height;

            inFile >> temp >> cellName >> width >> height;

            ll area = width*height;

            umap_cellNameToTechInfo[cellName][curTechName] = area;
        }
    }

    // read information of die
    ll dieWidth, dieHeight;
    inFile >> temp >> dieWidth >> dieHeight;
    dieMaximumArea = dieWidth*dieHeight;
    inFile >> temp >> dieTech[0] >> dieUtilization[0];
    inFile >> temp >> dieTech[1] >> dieUtilization[1];

    // read cell
    inFile >> temp >> netListcellNum;
    for(int i = 0; i < netListcellNum; i++)
    {
        string cellName, libCellName;

        inFile >> temp >> cellName >> libCellName;

        Cell* newCell = new Cell();
        newCell->name = cellName;
        for(auto psd : umap_cellNameToTechInfo[libCellName])
        {
            if(psd.first == dieTech[0])    // A的technology
                newCell->arr_cellArea[0] = psd.second;
            
            if(psd.first == dieTech[1])    // B的technology
                newCell->arr_cellArea[1] = psd.second;
        }
        newCell->cellAratio = (float)(newCell->arr_cellArea[0])/(newCell->arr_cellArea[0]+newCell->arr_cellArea[1]);
        newCell->isLock = false;

        umap_nameToNetlistCell[cellName] = newCell;
        vec_netlistCell.push_back(newCell);
    }    

    // read net
    inFile >> temp >> netNum;
    for(int i = 0; i < netNum; i++)
    {
        string netName;
        ll connectCellNum;

        inFile >> temp >> netName >> connectCellNum;

        Net* newNet = new Net();
        newNet->name = netName;

        for(int j = 0; j < connectCellNum; j++)
        {
            string cellName;

            inFile >> temp >> cellName;
            Cell* curCell = umap_nameToNetlistCell[cellName];

            newNet->vec_connectCell.push_back(curCell);
            curCell->vec_connectNet.push_back(newNet);
        }

        umap_nameToNet[netName] = newNet;
        vec_net.push_back(newNet);
    }

    inFile.close();

    return;
}

bool Partition::cmpPartition(bool (*func)(Cell*, Cell*))
{
    vector<Cell*> vec_dieACell, vec_dieBCell;

    dieArea[0] = 0, dieArea[1] = 0;

    sort(vec_netlistCell.begin(), vec_netlistCell.end(), func);

    // put all cell in dieA
    int cellIdx = 0;
    for(; cellIdx < vec_netlistCell.size()-1; cellIdx++)   // 留一個給dieB
    {
        if(((dieArea[0] + vec_netlistCell[cellIdx]->arr_cellArea[0])*100) > dieUtilization[0]*dieMaximumArea)
            break;

        dieArea[0] += vec_netlistCell[cellIdx]->arr_cellArea[0];
        vec_dieACell.push_back(vec_netlistCell[cellIdx]);
    }

    for(; cellIdx < vec_netlistCell.size(); cellIdx++)
    {
        dieArea[1] += vec_netlistCell[cellIdx]->arr_cellArea[1];
        vec_dieBCell.push_back(vec_netlistCell[cellIdx]);
    }

    // check utilization constraint
    if(((dieArea[0]*100) < dieUtilization[0]*dieMaximumArea) && ((dieArea[1]*100) < dieUtilization[1]*dieMaximumArea))
    {
        for(int i = 0; i < vec_dieACell.size(); i++)
            vec_dieACell[i]->curDie = 0;

        for(int i = 0; i < vec_dieBCell.size(); i++)
            vec_dieBCell[i]->curDie = 1;

        return true;
    }

    return false;
}

void Partition::initialPartition()
{
    bool initialGroup_ok = cmpPartition(&comparefun1);
    if(!initialGroup_ok)
    {
        initialGroup_ok = cmpPartition(&comparefun2);
    }
    if(!initialGroup_ok)
    {
        initialGroup_ok = cmpPartition(&comparefun3);
    }
    if(!initialGroup_ok)
    {
        initialGroup_ok = cmpPartition(&comparefun4);
    }
    if(!initialGroup_ok)
    {
        initialGroup_ok = cmpPartition(&comparefun5);
    }
    if(!initialGroup_ok)
    {
        initialGroup_ok = cmpPartition(&comparefun6);
    }
    if(!initialGroup_ok)
    {
        initialGroup_ok = cmpPartition(&comparefun7);
    }
    if(!initialGroup_ok)
    {
        initialGroup_ok = cmpPartition(&comparefun8);
    }

    return;
}

void Partition::initialNetDistributionGain()
{
    // 計算 distribution
    for(int i = 0; i < vec_net.size(); i++)
    {
        Net* curNet = vec_net[i];

        for(int j = 0; j < curNet->vec_connectCell.size(); j++)
            curNet->arr_netDistribution[curNet->vec_connectCell[j]->curDie]++;
    }

    // 計算gain
    for(int i = 0; i < vec_netlistCell.size(); i++)
    {
        Cell* curCell = vec_netlistCell[i];
        curCell->gain = 0;

        bool cellDie = curCell->curDie;
        if(cellDie == -1)
        {
            cerr << "Error in initialNetDistributionGain(), still ungroup cell" << endl;
            exit(1);
        }

        for(int j = 0; j < curCell->vec_connectNet.size(); j++)
        {
            Net* connectNet = curCell->vec_connectNet[j];

            if(connectNet->arr_netDistribution[cellDie] == 1)   // if F(n) = 1
                curCell->gain++;
            
            if(connectNet->arr_netDistribution[!cellDie] == 0)  // if T(n) = 0
                curCell->gain--;
        }

        map_bucketList[curCell->gain].insert(curCell);
    }

    return;
}

void Partition::iterrativeRevise()
{
    // count the cutSize
    ll cutSize = 0;
    for(Net* curNet : vec_net)
    {
        if(curNet->arr_netDistribution[0] > 0 && curNet->arr_netDistribution[1] > 0)
            cutSize++;
    }

    int passCnt = 0;
    while(1)
    {
        passCnt++;

        vector<pair<ll, Cell*>> vec_partitalSum;
        vector<Cell*> outofConstraint;

        // back up
        unordered_map<Cell*, ll> umap_originalGain;
        for(Cell* netlistCell : vec_netlistCell)
            umap_originalGain[netlistCell] = netlistCell->gain;

        unordered_map<Cell*, bool> umap_originalDie;
        for(Cell* netlistCell : vec_netlistCell)
            umap_originalDie[netlistCell] = netlistCell->curDie;

        unordered_map<Net*, vector<ll>> umap_originalDistribution;
        for(Net* curNet: vec_net)
        {
            vector<ll> copy(2);

            copy[0] = curNet->arr_netDistribution[0];
            copy[1] = curNet->arr_netDistribution[1];

            umap_originalDistribution[curNet] = copy;
        }

        ll oldArea_A = dieArea[0];
        ll oldArea_B = dieArea[1];
        
        while(!map_bucketList.empty())    // find a legal cell
        {
            auto end_time_parser = chrono::high_resolution_clock::now();
            auto duration_parser = chrono::duration_cast<chrono::milliseconds>(end_time_parser-start_time).count();

            if(duration_parser > 270000)
            {
                for(Cell* netlistCell : vec_netlistCell)
                {
                    netlistCell->isLock = false;
                    netlistCell->curDie = umap_originalDie[netlistCell];
                    netlistCell->gain = umap_originalGain[netlistCell];
                }

                for(Net* curNet : vec_net)
                {
                    curNet->arr_netDistribution[0] = umap_originalDistribution[curNet][0];
                    curNet->arr_netDistribution[1] = umap_originalDistribution[curNet][1];
                }

                return;
            }

            Cell* baseCell = *(map_bucketList.rbegin()->second.begin());
   
            ll cellGain = baseCell->gain;
            if(cellGain != map_bucketList.rbegin()->first)
            {
                cerr << "Error in gain ccal" << endl;
                exit(1);
            }

            bool cellDie = baseCell->curDie;
            ll moveArea = baseCell->arr_cellArea[!cellDie];
            ll removeArea = baseCell->arr_cellArea[cellDie];

            // check area constraint
            if(((dieArea[!cellDie] + moveArea)*100) <= dieUtilization[!cellDie]*dieMaximumArea)
            {
                vec_partitalSum.push_back({cellGain, baseCell});
                map_bucketList[cellGain].erase(baseCell);

                if(map_bucketList[cellGain].empty())
                    map_bucketList.erase(cellGain);
                
                // lock the base cell
                baseCell->isLock = true;

                // update connecting cell gain
                for(Net* connectNet : baseCell->vec_connectNet)
                {
                    for(Cell* connectCell : connectNet->vec_connectCell)
                    {
                        if(connectCell->isLock)
                            continue;

                        map_bucketList[connectCell->gain].erase(connectCell);

                        if(map_bucketList[connectCell->gain].empty())
                            map_bucketList.erase(connectCell->gain);

                        if(connectNet->arr_netDistribution[!cellDie] == 0)  // T(n) == 0，所有free cell gain+1
                            connectCell->gain++;

                        if(connectNet->arr_netDistribution[!cellDie] == 1 && connectCell->curDie == !cellDie)   // T(n) == 1，T內的cell gain-1
                            connectCell->gain--;
                    }

                    connectNet->arr_netDistribution[cellDie]--;
                    connectNet->arr_netDistribution[!cellDie]++;

                    for(Cell* connectCell : connectNet->vec_connectCell)
                    {
                        if(connectCell->isLock)
                            continue;

                        if(connectNet->arr_netDistribution[cellDie] == 0)  // F(n) == 0，所有free cell gain-1
                            connectCell->gain--;

                        if(connectNet->arr_netDistribution[cellDie] == 1 && connectCell->curDie == cellDie)   // F(n) == 1，F內的cell gain+1
                            connectCell->gain++;

                        map_bucketList[connectCell->gain].insert(connectCell);
                    }
                }

                dieArea[!cellDie] = dieArea[!cellDie] + moveArea;
                dieArea[cellDie] = dieArea[cellDie] - removeArea;
                baseCell->curDie = !baseCell->curDie;
            }
            else
            {
                baseCell->isLock = true;
                outofConstraint.push_back(baseCell);
                map_bucketList[cellGain].erase(baseCell);

                if(map_bucketList[cellGain].empty())
                    map_bucketList.erase(cellGain);
            }

            // count the cutSize
            ll cutSize = 0;
            for(Net* curNet : vec_net)
            {
                if(curNet->arr_netDistribution[0] > 0 && curNet->arr_netDistribution[1] > 0)
                    cutSize++;
            }
        }

        // fimd maximum partial sum
        for(int i = 1; i < vec_partitalSum.size(); i++)
            vec_partitalSum[i].first += vec_partitalSum[i-1].first; 

        ll maximumGain = LLONG_MIN;
        int maxIdx = -1;
        for(int i = 0; i < vec_partitalSum.size(); i++)
        {
            if(vec_partitalSum[i].first >= maximumGain)
            {
                maximumGain = vec_partitalSum[i].first;
                maxIdx = i;
            }
        } 

        for(Cell* netlistCell : vec_netlistCell)
        {
            netlistCell->isLock = false;
            netlistCell->curDie = umap_originalDie[netlistCell];
            netlistCell->gain = umap_originalGain[netlistCell];
        }

        for(Net* curNet : vec_net)
        {
            curNet->arr_netDistribution[0] = umap_originalDistribution[curNet][0];
            curNet->arr_netDistribution[1] = umap_originalDistribution[curNet][1];
        }

        auto end_time_parser = chrono::high_resolution_clock::now();
        auto duration_parser = chrono::duration_cast<chrono::milliseconds>(end_time_parser-start_time).count();

        if(duration_parser > 270000 || maximumGain <= 0)
        {
            for(Cell* netlistCell : vec_netlistCell)
            {
                netlistCell->isLock = false;
                netlistCell->curDie = umap_originalDie[netlistCell];
                netlistCell->gain = umap_originalGain[netlistCell];
            }

            for(Net* curNet : vec_net)
            {
                curNet->arr_netDistribution[0] = umap_originalDistribution[curNet][0];
                curNet->arr_netDistribution[1] = umap_originalDistribution[curNet][1];
            }

            return;
        }

        dieArea[0] = oldArea_A, dieArea[1] = oldArea_B;


        // 重新換那些cell
        for(int i = 0; i <= maxIdx; i++)
        {
            Cell* baseCell = vec_partitalSum[i].second;

            ll cellGain = baseCell->gain;
            bool cellDie = baseCell->curDie;
            ll moveArea = baseCell->arr_cellArea[!cellDie];
            ll removeArea = baseCell->arr_cellArea[cellDie];

            // update connecting cell gain
            for(Net* connectNet : baseCell->vec_connectNet)
            {
                connectNet->arr_netDistribution[cellDie]--;
                connectNet->arr_netDistribution[!cellDie]++;
            }

            dieArea[!cellDie] = dieArea[!cellDie] + moveArea;
            dieArea[cellDie] = dieArea[cellDie] - removeArea;
            baseCell->curDie = !baseCell->curDie;
        }

        // 計算gain
        for(int i = 0; i < vec_netlistCell.size(); i++)
        {
            Cell* curCell = vec_netlistCell[i];
            curCell->gain = 0;

            bool cellDie = curCell->curDie;

            for(int j = 0; j < curCell->vec_connectNet.size(); j++)
            {
                Net* connectNet = curCell->vec_connectNet[j];

                if(connectNet->arr_netDistribution[cellDie] == 1)   // if F(n) = 1
                    curCell->gain++;
                
                if(connectNet->arr_netDistribution[!cellDie] == 0)  // if T(n) = 0
                    curCell->gain--;
            }

            map_bucketList[curCell->gain].insert(curCell);
        }

        // 解鎖所有cell
        for(Cell* curCell : vec_netlistCell)
            curCell->isLock = false;

        // count the cutSize
        ll cutSize = 0;
        for(Net* curNet : vec_net)
        {
            if(curNet->arr_netDistribution[0] > 0 && curNet->arr_netDistribution[1] > 0)
                cutSize++;
        }
    }
}

void Partition::writeOutputFile(string path)
{
    ofstream outFile;

    outFile.open(path);
    if(!outFile.is_open())
    {
        cerr << "error in writeOutputFile, can't open file" << endl;
        exit(1);
    }

    // count the cutSize
    ll cutSize = 0;
    for(Net* curNet : vec_net)
    {
        if(curNet->arr_netDistribution[0] > 0 && curNet->arr_netDistribution[1] > 0)
            cutSize++;
    }
    outFile << "CutSize " << cutSize << endl;
    cout << "CutSize " << cutSize << endl;

    // count die A number
    ll dieANum = 0, dieBNum = 0;
    for(Cell* curCell : vec_netlistCell)
    {
        if(curCell->curDie == 0)
            dieANum++;
        else if(curCell->curDie == 1)
            dieBNum++;
    } 

    outFile << "DieA " << dieANum << endl;
    for(Cell* curCell : vec_netlistCell)
    {
        if(curCell->curDie == 0)
            outFile << curCell->name << endl;
    } 

    outFile << "DieB " << dieBNum << endl;
    for(Cell* curCell : vec_netlistCell)
    {
        if(curCell->curDie == 1)
            outFile << curCell->name << endl;
    } 

    outFile.close();
    return;
}
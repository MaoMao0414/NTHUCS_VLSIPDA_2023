#include "bits/stdc++.h"
#include <chrono>
#include "module.h"
#include "net.h"
#include "tree.h"
#include "cmpFun.hpp"

using namespace std;

#define ll long 

double LB = 0.5, UB = 2;   // soft modules的aspect ratio比例

std::chrono::time_point<std::chrono::high_resolution_clock> start_time;

struct ContourNode
{
    ll x1, x2, y;
};

class FloorPlan
{
    public:
        void initialBtree();    // 建一顆initial binary tree (complete)
        void BtreeTofp();

        double calCost();     // 計算當前b tree產生的floorplan的cost function 並 return cost (TODO: 調整cost function)
        long long calHPWL();    // 計算weighted HPWL

        void SA(int timeLimit);
        void perturb();

        // 備份與還原當前狀態
        void copy();
        void restore();

        // 更新最佳解的狀態
        void updateBest();
        void copyfromBestinitial();

        // 檢查此次擺放是否為feasible solution
        bool isLegal(Module*, ll, ll);
        bool checkNoOverlap();
        bool checkFeasible();

        // 讀寫檔相關
        void parseInputFile(string);
        void writeOutputFile(string);

        // Debug用
        void printCurrentState();

        bool getTimeOut(){return timeOut;}
        bool getbestLegal(){return bestLegal;}

        // if有module violate constraint，嘗試更改shape後放到每一個點
        void refine();
    
    private:
        ll numSoftModule, numFixModule, numNet;
        ll chipWidth, chipHeight;

        vector<Net*> vec_net;
        vector<Module*> vec_allModule;
        vector<Module*> vec_initialFixedModule;
        vector<Module*> vec_floatModule;
        vector<Module*> vec_copySoftModule;

        unordered_map<string, Module*> umap_nameToModule;
        unordered_map<string, Module*> umap_nameTobestResModule;
        unordered_map<Module*, Node*> umap_moduleToTreeNodeIdx;

        // copy 用
        unordered_map<Node*, Node*> umap_copyNodeToOldNode;
        unordered_map<Node*, Node*> umap_bestNodeToOldNode;

        // 這兩個變數用來記錄當前floorplan的最上和最右邊界
        ll placeTop, placeRight;

        // 這兩個變數用來normalize
        ll initialOutsum, initialHPWL;

        // 紀錄最佳解
        double bestCost;
        bool bestLegal;
        vector<Module*> vec_bestResult;

        // b Tree兄弟
        Node* bTreeRoot;
        Node* copyTreeRoot;
        Node* bestTreeRoot;

        bool timeOut;
};

void FloorPlan::parseInputFile(string filePath)
{
    //cout << "==============================" << endl;
    //cout << "Parse the input file" << endl;

    ifstream inFile;
    inFile.open(filePath);

    if(!inFile.is_open())
    {
        cerr << "Error in parseInputFile(), no input file" << endl;
        exit(1);
    }  

    string temp;
    // first Read the chipsize
    inFile >> temp >> chipWidth >> chipHeight;

    // then read the soft modules
    inFile >> temp >> numSoftModule;
    for(ll i = 0; i < numSoftModule; i++)
    {
        string moduleName;
        ll minArea;

        inFile >> temp >> moduleName >> minArea;

        // create a new module object
        Module* newModule = new Module();
        newModule->name = moduleName;
        newModule->minArea = minArea;
        newModule->isFixed = false;

        long long minHeight = sqrt(minArea*0.5) + 1;
        long long maxHeight = sqrt(minArea*2) - 1;

        long long interval = (maxHeight - minHeight)/100;
        if(interval == 0)
            interval = 1;
        
        for(ll h = minHeight; h < maxHeight; h += interval)
        {
            ll w = minArea/h + 1;

            // check
            if(w*h < minArea || h > w*2 || h*2 < w)
                cout << "Warning: Violate Soft Modules Contraint" << endl;
            else
                newModule->vec_allpossibleShape.push_back({w, h});
        }

        ll minWidth = minArea/maxHeight + 1;
        // check
        if(minWidth*maxHeight < minArea || maxHeight > minWidth*2 || maxHeight*2 < minWidth)
            cout << "Warning: Violate Soft Modules Contraint" << endl;
        else
            newModule->vec_allpossibleShape.push_back({minWidth, maxHeight});


        // random決定一個initial shape
        int shapeCnt = newModule->vec_allpossibleShape.size();
        int shapeIdx = rand()%shapeCnt;

        newModule->width = (newModule->vec_allpossibleShape[shapeIdx]).first;   
        newModule->height = (newModule->vec_allpossibleShape[shapeIdx]).second;  

        vec_allModule.push_back(newModule);
        vec_floatModule.push_back(newModule);
        umap_nameToModule[moduleName] = newModule;
    }

    // then read the fixed modules
    inFile >> temp >> numFixModule;
    for(ll i = 0; i < numFixModule; i++)
    {
        string moduleName;
        ll cordX, cordY, width, height;

        inFile >> temp >> moduleName >> cordX >> cordY >> width >> height;

        // create a new module object
        Module* newModule = new Module();
        newModule->name = moduleName;
        newModule->cordX = cordX;
        newModule->cordY = cordY;
        newModule->width = width;
        newModule->height = height;
        newModule->rightCordX = cordX + width;
        newModule->rightCordY = cordY + height;
        newModule->isFixed = true;

        vec_allModule.push_back(newModule);
        vec_initialFixedModule.push_back(newModule);
        umap_nameToModule[moduleName] = newModule;
    }

    // finally, read the nets
    inFile >> temp >> numNet;
    for(ll i = 0; i < numNet; i++)
    {
        string module1Name, module2Name;
        ll netWeight;

        inFile >> temp >> module1Name >> module2Name >> netWeight;

        // create a new net object
        Net* newNet = new Net();
        newNet->module1Name = module1Name;
        newNet->module2Name = module2Name;
        newNet->weight = netWeight;

        vec_net.push_back(newNet);
    }

    // initialize two vector for copy
    bestCost = 0;
    timeOut = false;
    for(int i = 0; i < numSoftModule; i++)
    {
        Module* m1 = new Module();
        Module* m2 = new Module();

        vec_copySoftModule.push_back(m1);
        vec_bestResult.push_back(m2);

        umap_nameTobestResModule[vec_floatModule[i]->name] = m2;
    }

    // 也要把fixed module的name丟進去umap_nameTobestResModule
    for(Module* fixedModule : vec_initialFixedModule)
        umap_nameTobestResModule[fixedModule->name] = fixedModule;

    inFile.close();

    //cout << endl;
    //cout << "=====================" << endl;
    //cout << " initial module info" << endl;
    for(int i = 0; i < vec_floatModule.size(); i++)
    {
        Module* softModule = vec_floatModule[i];

        //cout << "Name: " << softModule->name << " " << "width: " << softModule->width << " " << "height: " << softModule->height << endl;
    }
    //cout << "=====================" << endl;
    //cout << endl;

    auto end_time_parser = chrono::high_resolution_clock::now();
    auto duration_parser = chrono::duration_cast<chrono::milliseconds>(end_time_parser-start_time).count();
    //cout << "input file spend time: " << duration_parser << endl;
    //cout << "==============================" << endl;
    //cout << endl;

    return;
}

// initial : 建一顆complete B* tree
void FloorPlan::initialBtree()
{
   // cout << "==============================" << endl;
   // cout << "Build initial binary tree" << endl;

    vector<int> vec_initial(numSoftModule);
    for(int i = 0; i < vec_initial.size(); i++)
        vec_initial[i] = i;

    int idx = 0;
    Node* rootNode = new Node();
    rootNode->moduleIdx = vec_initial[idx++];
    rootNode->left = nullptr;
    rootNode->right = nullptr;

    umap_moduleToTreeNodeIdx[vec_floatModule[rootNode->moduleIdx]] = rootNode;

    bTreeRoot = rootNode;

    queue<Node*> queue_bfs;
    queue_bfs.push(rootNode);

    while(!queue_bfs.empty())
    {
        Node* curNode = queue_bfs.front();
        queue_bfs.pop();

        if(idx >= vec_floatModule.size())
            break;
        
        Node* leftChildNode = new Node();
        leftChildNode->moduleIdx = vec_initial[idx++];;
        leftChildNode->left = nullptr;
        leftChildNode->right = nullptr;

        umap_moduleToTreeNodeIdx[vec_floatModule[leftChildNode->moduleIdx]] = leftChildNode;
        curNode->left = leftChildNode;

        queue_bfs.push(leftChildNode);

        if(idx >= vec_floatModule.size())
            break;
            
        Node* rightChildNode = new Node();
        rightChildNode->moduleIdx = vec_initial[idx++];;
        rightChildNode->left = nullptr;
        rightChildNode->right = nullptr;

        umap_moduleToTreeNodeIdx[vec_floatModule[rightChildNode->moduleIdx]] = rightChildNode;
        curNode->right = rightChildNode;

        queue_bfs.push(rightChildNode);
    }

    initialHPWL = 0;
    initialOutsum = -1;

    auto end_time_parser = chrono::high_resolution_clock::now();
    auto duration_parser = chrono::duration_cast<chrono::milliseconds>(end_time_parser-start_time).count();
    //cout << "initialBtree spend time: " << duration_parser << endl;
    //cout << "==============================" << endl;
    //cout << endl;

    return;
}

// 從B* tree推出一個floorplanning
void FloorPlan::BtreeTofp()
{
    // reset所有soft module的座標
    for(Module* softModule : vec_floatModule)
    {
        softModule->cordX = 0;
        softModule->cordY = 0;
        softModule->isFixed = false;
    }

    // 計算新的邊界值是多少
    placeTop = 0, placeRight = 0;

    stack<Node*> stack_dfs;
    stack_dfs.push(bTreeRoot);

    // 先決定root的cord Y => x = 0，如果不行就往上擺
    // 代表root的module一定要合法
    Module* rootModule = vec_floatModule[bTreeRoot->moduleIdx];

    for(ll attempX = 0; attempX + rootModule->width <= chipWidth; attempX++)
    {
        bool find = false;
        for(ll attempY = 0; attempY + rootModule->height <= chipHeight; attempY++)
        {
            if(isLegal(rootModule, attempX, attempY))
            {
                rootModule->cordX = attempX;
                rootModule->cordY = attempY;
                rootModule->isFixed = true;
                
                find = true;

                break;
            }
        }

        if(find)
            break;
    }

    while(!stack_dfs.empty())
    {
        Node* curNode = stack_dfs.top();

        Module* curModule = vec_floatModule[curNode->moduleIdx];
        stack_dfs.pop();

        if(curModule->cordX + curModule->width > placeRight)
            placeRight = curModule->cordX + curModule->width;

        if(curModule->cordY + curModule->height > placeTop)
            placeTop = curModule->cordY + curModule->height;

        
        if(curNode->right)
        {
            Module* rightChildModule = vec_floatModule[curNode->right->moduleIdx];

             // ORIGIN VERSION
            ll attempX = curModule->cordX;

            for(ll attempY = curModule->cordY + curModule->height; attempY <= LLONG_MAX; attempY++)
            {
                bool find = false;
                if(isLegal(rightChildModule, attempX, attempY))
                {
                    rightChildModule->cordX = attempX;
                    rightChildModule->cordY = attempY;
                    rightChildModule->isFixed = true;
                    
                    find = true;

                    break;
                }  

                if(find)
                    break;              
            }

            stack_dfs.push(curNode->right);
        }

        // traverse left right child
        if(curNode->left)
        {
            Module* leftChildModule = vec_floatModule[curNode->left->moduleIdx];

            // ORIGIN VERSION
            ll attempX = curModule->cordX + curModule->width;
            
            for(ll attempY = 0; attempY <= LLONG_MAX; attempY++)
            {
                bool find = false;
                if(isLegal(leftChildModule, attempX, attempY))
                {
                    leftChildModule->cordX = attempX;
                    leftChildModule->cordY = attempY;
                    leftChildModule->isFixed = true;
                    
                    find = true;

                    break;
                }   

                if(find)
                    break;             
            }

            stack_dfs.push(curNode->left);
        }       
    }

    return;
} 

double FloorPlan::calCost()
{
    // 產生對應的floorplan
    BtreeTofp();

    // 計算它的weight HPWL
    ll weightHPWL = calHPWL();

    // 檢查是否outofbound
    ll outofBoundX = 0;
    if(placeRight-chipWidth > 0)
        outofBoundX = placeRight-chipWidth;

    ll outofBoundY = 0;
    if(placeTop-chipHeight > 0)
        outofBoundY = placeTop-chipHeight;

    ll outSum = outofBoundX + outofBoundY;

    if(initialHPWL == 0)
        initialHPWL = weightHPWL;
    
    if(initialOutsum == -1)
    {
        if(outSum == 0)
            initialOutsum = 1;
        else
            initialOutsum = outSum;
    }

    double cost = 50*((double)outSum/initialOutsum) + ((double)weightHPWL/initialHPWL);

    return cost;
}

long long FloorPlan::calHPWL()
{
    ll temp = 0;

    for(Net* n : vec_net)
    {
        Module* m1 = umap_nameToModule[n->module1Name];
        Module* m2 = umap_nameToModule[n->module2Name];

        ll weight = n->weight;

        ll center1X = (m1->cordX + (m1->cordX + m1->width))/2;
        ll center1Y = (m1->cordY + (m1->cordY + m1->height))/2;
        ll center2X = (m2->cordX + (m2->cordX + m2->width))/2;
        ll center2Y = (m2->cordY + (m2->cordY + m2->height))/2;

        temp += weight*(abs(center1X-center2X) + abs(center1Y-center2Y));
    }

    return temp;
}

void FloorPlan::updateBest()
{
    for(int i = 0; i < vec_floatModule.size(); i++)
    {
        vec_bestResult[i]->name = vec_floatModule[i]->name;
        vec_bestResult[i]->cordX = vec_floatModule[i]->cordX;
        vec_bestResult[i]->cordY = vec_floatModule[i]->cordY;
        vec_bestResult[i]->width = vec_floatModule[i]->width;
        vec_bestResult[i]->height = vec_floatModule[i]->height;
    }

    // copy B* tree
    umap_bestNodeToOldNode.clear();

    bestTreeRoot = nullptr;
    Node* copyRoot = new Node();
    bestTreeRoot = copyRoot;

    queue<Node*> queue_bfs;
    queue<Node*> queue_bfsCopy; 
    queue_bfs.push(bTreeRoot);
    queue_bfsCopy.push(bestTreeRoot);

    while(!queue_bfs.empty())
    {
        Node* frontNode = queue_bfs.front();
        Node* copyNode = queue_bfsCopy.front();

        umap_bestNodeToOldNode[copyNode] = frontNode;

        copyNode->left = nullptr;
        copyNode->right = nullptr;

        queue_bfs.pop();
        queue_bfsCopy.pop();
        
        if(frontNode->left)
        {
            Node* newCopyNode = new Node();
            copyNode->left = newCopyNode;

            queue_bfs.push(frontNode->left);
            queue_bfsCopy.push(copyNode->left);
        }

        if(frontNode->right)
        {
            Node* newCopyNode = new Node();
            copyNode->right = newCopyNode;

            queue_bfs.push(frontNode->right);
            queue_bfsCopy.push(copyNode->right);
        }
    }
}

void FloorPlan::SA(int timeLimit)
{
    //cout << "==================" << endl;
    //cout << "   start SA       " << endl;

    // hyperparameter
    double initialTemperature = 1000;
    double coolingRate = 0.85;           
    double terminateTemperature = 0.1;
    double maxRejectRate = 0.95;             // 每個溫度失敗上限
    double N = 8*numSoftModule;  // 每個溫度嘗試次數

    // 建一個initial solution
    if(bestCost == 0)
        initialBtree();
    else
        copyfromBestinitial();


    bestCost = calCost();
    bestLegal =  checkFeasible();
    updateBest();

    double curTemperature = initialTemperature;
    while(curTemperature >= terminateTemperature)
    {
        int totalAttempt = 0, reject = 0, upHill = 0;

        double curCost = calCost();
        while(upHill <= N && totalAttempt <= 2*N)
        {   
            auto end_time_parser = chrono::high_resolution_clock::now();
            auto duration_parser = chrono::duration_cast<chrono::milliseconds>(end_time_parser-start_time).count();

            // REMEMBER TO ADD BACK !!!
            if(duration_parser > timeLimit)
            {
                timeOut = true;
                //cout << duration_parser << ": " << endl;
                //cout << "TLE" << endl;
                return;
            }

            // 複製當前state
            copy();

            // 從四種move裡面選一種進行perturb
            perturb();

            // 計算perturb後的cost
            double newCost = calCost();
            if(newCost < curCost)   // down-hill move
            {
                curCost = newCost;
                if(newCost < bestCost)
                {
                    bestCost = newCost;
                    bestLegal = checkFeasible();

                    updateBest();
                }
            }
            else        // up-hill move
            {
                double changeProbability = exp(-(newCost-curCost)/curTemperature);
                double p = (double)rand()/(RAND_MAX+1.0);

                // accept
                if(p < changeProbability)
                    curCost = newCost;
                else
                {
                    restore();
                    reject++;
                }

                upHill++;
            }

            totalAttempt++;
        }

        //cout << "Temperature: " << curTemperature << " " << "best result: " << bestCost << endl;
        
        curTemperature = curTemperature*coolingRate;

        if((double)reject/totalAttempt > maxRejectRate)
        {
            //cout << "reject rate out of bound!!" << endl;
            break;
        }
    }

    auto end_time_parser = chrono::high_resolution_clock::now();
    auto duration_parser = chrono::duration_cast<chrono::milliseconds>(end_time_parser-start_time).count();
    //cout << "SA Time: " << duration_parser << endl;

    return;
}

void FloorPlan::perturb()
{
    // 從四種選擇中挑一個
    int perturbChoose = rand()%4;

    // rotate一個module
    if(perturbChoose == 0)
    {
        // 隨機選擇一個module
        int selectModuleIdx = rand()%numSoftModule;
        Module* selectModule = vec_floatModule[selectModuleIdx];

        swap(selectModule->width, selectModule->height);
    }
    else if(perturbChoose == 1)     // 搬動B* tree上的一個node到其他地方
    {
        // 決定要搬的node所存的module idx => 等同於挑選一個node

        // 隨機選擇一個module Idx
        int selectModuleIdx = rand()%numSoftModule;

        // bfs 找出這個node
        Node* selNode = nullptr;
        Node* selNodeParent = nullptr;

        queue<Node*> queue_bfs;
        queue_bfs.push(bTreeRoot);

        // TRAVERSE(bfs)
        while(!queue_bfs.empty())
        {
            Node* frontNode = queue_bfs.front();
            queue_bfs.pop();

            // 這個代表此node選到root
            if(frontNode->moduleIdx == selectModuleIdx)
            {
                selNode = frontNode;
                break;
            }

            if(frontNode->left)
            {
                if(frontNode->left->moduleIdx == selectModuleIdx)
                {
                    selNodeParent = frontNode;
                    selNode = frontNode->left;

                    break;
                }

                queue_bfs.push(frontNode->left);
            }

            if(frontNode->right)
            {
                if(frontNode->right->moduleIdx == selectModuleIdx)
                {
                    selNodeParent = frontNode;
                    selNode = frontNode->right;

                    break;
                }

                queue_bfs.push(frontNode->right);
            }
        }

        // 清空queue
        queue_bfs = queue<Node*>();

        // 刪除此node
        if(!selNode->left && !selNode->right)   // case 1: 沒有child
        {
            if(selNodeParent->left == selNode)
                selNodeParent->left = nullptr;
            else
                selNodeParent->right = nullptr;
        }
        else if(selNode->left && !selNode->right)  // case 2: 一個child
        {
            if(selNode == bTreeRoot)    // speical case
                bTreeRoot = selNode->left;
            else
            {
                if(selNodeParent->left == selNode)
                    selNodeParent->left = selNode->left;
                else
                    selNodeParent->right = selNode->left;                
            }

            selNode->left = nullptr;
        }
        else if(!selNode->left && selNode->right)  // case 2: 一個child
        {
            if(selNode == bTreeRoot)    // speical case
                bTreeRoot = selNode->right;
            else
            {
                if(selNodeParent->left == selNode)
                    selNodeParent->left = selNode->right;
                else
                    selNodeParent->right = selNode->right;                
            }

            selNode->right = nullptr;
        }
        else    // case 3: 二個child
        {
            int moveDirection = rand()%2;   // 決定要用left child或是right child往上取代，0 => right, 1 => left

            if(moveDirection == 0)  // right child
            {
                if(selNode == bTreeRoot)    // speical case
                    bTreeRoot = selNode->right;
                else
                {
                    if(selNodeParent->left == selNode)
                        selNodeParent->left = selNode->right;
                    else
                        selNodeParent->right = selNode->right;                
                }   

                vector<Node*> vec_traverseNode;
                vec_traverseNode.push_back(selNode);

                Node* traverseNode = selNode;
                while(traverseNode->right)
                {
                    vec_traverseNode.push_back(traverseNode->right);
                    traverseNode = traverseNode->right;
                }

                if(traverseNode->left)
                {
                    traverseNode->right = traverseNode->left;
                    traverseNode->left = nullptr;
                }

                for(int i = vec_traverseNode.size()-1; i >= 1; i--)
                {
                    Node* curNode = vec_traverseNode[i];
                    Node* preNode = vec_traverseNode[i-1];

                    curNode->left = preNode->left;
                }

                selNode->left = nullptr;
                selNode->right = nullptr;             
            }
            else        // left child
            {
                if(selNode == bTreeRoot)    // speical case
                    bTreeRoot = selNode->left;
                else
                {
                    if(selNodeParent->left == selNode)
                        selNodeParent->left = selNode->left;
                    else
                        selNodeParent->right = selNode->left;                
                }   

                vector<Node*> vec_traverseNode;
                vec_traverseNode.push_back(selNode);

                Node* traverseNode = selNode;
                while(traverseNode->left)
                {
                    vec_traverseNode.push_back(traverseNode->left);
                    traverseNode = traverseNode->left;
                }

                if(traverseNode->right)
                {
                    traverseNode->left = traverseNode->right;
                    traverseNode->right = nullptr;
                }

                for(int i = vec_traverseNode.size()-1; i >= 1; i--)
                {
                    Node* curNode = vec_traverseNode[i];
                    Node* preNode = vec_traverseNode[i-1];

                    curNode->right = preNode->right;
                }

                selNode->left = nullptr;
                selNode->right = nullptr;
            }
        }

        // 決定要插在哪個module的後面，if insertIdx == selNode->moduleIdx，則插在root
        int insertIdx = rand()%numSoftModule;

        if(insertIdx == selNode->moduleIdx)
        {
            int dir = rand()%2;

            if(dir == 0)    // 插在左邊
                selNode->left = bTreeRoot;
            else if(dir == 1)   // 插在右邊
                selNode->right = bTreeRoot;

            bTreeRoot = selNode;
        }
        else
        {
            queue_bfs.push(bTreeRoot);

            while(!queue_bfs.empty())
            {
                Node* frontNode = queue_bfs.front();
                queue_bfs.pop();

                if(frontNode->moduleIdx == insertIdx)
                {
                    int dir = rand()%2;

                    if(dir == 0)    // 插在左邊
                    {
                        if(frontNode->left)
                            selNode->left = frontNode->left;
                        
                        frontNode->left = selNode;
                    }
                    else if(dir == 1)   // 插在右邊
                    {
                        if(frontNode->right)
                            selNode->right = frontNode->right;
                        
                        frontNode->right = selNode;
                    }
                    
                    break;
                }

                if(frontNode->left)
                    queue_bfs.push(frontNode->left);
                
                if(frontNode->right)
                    queue_bfs.push(frontNode->right);
            }
        }
    }
    else if(perturbChoose == 2)     // swap B* tree上的兩個module
    {   
        int selectModuleIdx1 = rand()%numSoftModule;
        int selectModuleIdx2 = selectModuleIdx1;
        while(selectModuleIdx2 == selectModuleIdx1)
            selectModuleIdx2 = rand()%numSoftModule;

        // search 找出這兩個node
        Node *node1 = nullptr, *node1Parent = nullptr;
        Node *node2 = nullptr, *node2Parent = nullptr;

        queue<Node*> queue_bfs;
        queue_bfs.push(bTreeRoot);

        while(!queue_bfs.empty())
        {
            Node* frontNode = queue_bfs.front();
            queue_bfs.pop();

            if(frontNode->moduleIdx == selectModuleIdx1)
                node1 = frontNode;
            
            if(frontNode->moduleIdx == selectModuleIdx2)
                node2 = frontNode;
            
            if(frontNode->left)
            {
                if(frontNode->left->moduleIdx == selectModuleIdx1)
                    node1Parent = frontNode;
                
                if(frontNode->left->moduleIdx == selectModuleIdx2)
                    node2Parent = frontNode;
                
                queue_bfs.push(frontNode->left);
            }

            if(frontNode->right)
            {
                if(frontNode->right->moduleIdx == selectModuleIdx1)
                    node1Parent = frontNode;
                
                if(frontNode->right->moduleIdx == selectModuleIdx2)
                    node2Parent = frontNode;
                
                queue_bfs.push(frontNode->right);
            }            
        }

        Node* tempNode1 = node1, *tempNode2 = node2;

        // 處理parent
        if(node1 == bTreeRoot)
        {   
            // node2的parent連到node1
            if(node2Parent->left == tempNode2)
                node2Parent->left = tempNode1;
            else
                node2Parent->right = tempNode1;

            // node2 成為new treeRoot
            bTreeRoot = tempNode2;
        }
        else if(node2 == bTreeRoot)
        {
            // node1的parent連到node2
            if(node1Parent->left == tempNode1)
                node1Parent->left = tempNode2;
            else
                node1Parent->right = tempNode2;

            // node1 成為new treeRoot
            bTreeRoot = tempNode1;            
        }
        else
        {
            // node1的parent連到node2
            if(node1Parent->left == tempNode1)
                node1Parent->left = tempNode2;
            else
                node1Parent->right = tempNode2;   

            // node2的parent連到node1
            if(node2Parent->left == tempNode2)
                node2Parent->left = tempNode1;
            else
                node2Parent->right = tempNode1;         
        }

        // 處理child
        Node* node1Left = node1->left;
        Node* node1Right = node1->right;
        Node* node2Left = node2->left;
        Node* node2Right = node2->right;

        node1->left = node2Left;
        node1->right = node2Right;

        node2->left = node1Left;
        node2->right = node1Right;
    }
    else if(perturbChoose == 3)     // 將一個module變換另一種shape
    {
        // 隨機選擇一個module
        int selectModuleIdx = rand()%numSoftModule;

        Module* selectModule = vec_floatModule[selectModuleIdx];

        int shapeNum = selectModule->vec_allpossibleShape.size();
        int shapeIdx = rand()%shapeNum;

        selectModule->width = (selectModule->vec_allpossibleShape[shapeIdx]).first;
        selectModule->height = (selectModule->vec_allpossibleShape[shapeIdx]).second;
    }
}

void FloorPlan::copy()
{
    for(int i = 0; i < vec_floatModule.size(); i++)
    {
        vec_copySoftModule[i]->cordX = vec_floatModule[i]->cordX;
        vec_copySoftModule[i]->cordY = vec_floatModule[i]->cordY;
        vec_copySoftModule[i]->height = vec_floatModule[i]->height;
        vec_copySoftModule[i]->width = vec_floatModule[i]->width;
    }

    // copy B* tree
    umap_copyNodeToOldNode.clear();

    copyTreeRoot = nullptr;
    Node* copyRoot = new Node();
    copyTreeRoot = copyRoot;

    queue<Node*> queue_bfs;
    queue<Node*> queue_bfsCopy; 
    queue_bfs.push(bTreeRoot);
    queue_bfsCopy.push(copyRoot);

    while(!queue_bfs.empty())
    {
        Node* frontNode = queue_bfs.front();
        Node* copyNode = queue_bfsCopy.front();

        umap_copyNodeToOldNode[copyNode] = frontNode;

        copyNode->left = nullptr;
        copyNode->right = nullptr;

        queue_bfs.pop();
        queue_bfsCopy.pop();
        
        if(frontNode->left)
        {
            Node* newCopyNode = new Node();
            copyNode->left = newCopyNode;

            queue_bfs.push(frontNode->left);
            queue_bfsCopy.push(copyNode->left);
        }

        if(frontNode->right)
        {
            Node* newCopyNode = new Node();
            copyNode->right = newCopyNode;

            queue_bfs.push(frontNode->right);
            queue_bfsCopy.push(copyNode->right);
        }
    }

    return;
}

void FloorPlan::restore()
{
    for(int i = 0; i < vec_floatModule.size(); i++)
    {
        vec_floatModule[i]->cordX = vec_copySoftModule[i]->cordX;
        vec_floatModule[i]->cordY = vec_copySoftModule[i]->cordY;
        vec_floatModule[i]->height = vec_copySoftModule[i]->height;
        vec_floatModule[i]->width = vec_copySoftModule[i]->width;
    }

    // restore B* tree
    bTreeRoot = umap_copyNodeToOldNode[copyTreeRoot];

    queue<Node*> queue_bfscp;
    queue<Node*> queue_bfsori;
    queue_bfscp.push(copyTreeRoot);
    queue_bfsori.push(bTreeRoot);

    while(!queue_bfscp.empty())
    {
        Node* copyNode = queue_bfscp.front();
        Node* originNode = queue_bfsori.front();

        originNode->left = nullptr;
        originNode->right = nullptr;

        queue_bfscp.pop();
        queue_bfsori.pop();
        
        if(copyNode->left)
        {
            originNode->left = umap_copyNodeToOldNode[copyNode->left];

            queue_bfscp.push(copyNode->left);
            queue_bfsori.push(originNode->left);
        }

        if(copyNode->right)
        {
            originNode->right = umap_copyNodeToOldNode[copyNode->right];

            queue_bfscp.push(copyNode->right);
            queue_bfsori.push(originNode->right);
        }
    }

    return;   
}

void FloorPlan::copyfromBestinitial()
{
    for(int i = 0; i < vec_floatModule.size(); i++)
    {
        vec_floatModule[i]->cordX = vec_bestResult[i]->cordX;
        vec_floatModule[i]->cordY = vec_bestResult[i]->cordY;
        vec_floatModule[i]->height = vec_bestResult[i]->height;
        vec_floatModule[i]->width = vec_bestResult[i]->width;
    }

    // restore B* tree
    bTreeRoot = umap_bestNodeToOldNode[bestTreeRoot];

    queue<Node*> queue_bfscp;
    queue<Node*> queue_bfsori;
    queue_bfscp.push(bestTreeRoot);
    queue_bfsori.push(bTreeRoot);

    while(!queue_bfscp.empty())
    {
        Node* copyNode = queue_bfscp.front();
        Node* originNode = queue_bfsori.front();

        originNode->left = nullptr;
        originNode->right = nullptr;

        queue_bfscp.pop();
        queue_bfsori.pop();
        
        if(copyNode->left)
        {
            originNode->left = umap_bestNodeToOldNode[copyNode->left];

            queue_bfscp.push(copyNode->left);
            queue_bfsori.push(originNode->left);
        }

        if(copyNode->right)
        {
            originNode->right = umap_bestNodeToOldNode[copyNode->right];

            queue_bfscp.push(copyNode->right);
            queue_bfsori.push(originNode->right);
        }
    }

    initialHPWL = 0;
    initialOutsum = -1;

    return;     
}

bool FloorPlan::checkFeasible()
{
    // 檢查是否outofbound
    ll outofBoundX = 0;
    if(placeRight-chipWidth > 0)
        outofBoundX = placeRight-chipWidth;

    ll outofBoundY = 0;
    if(placeTop-chipHeight > 0)
        outofBoundY = placeTop-chipHeight;

    return (outofBoundX == 0) && (outofBoundY == 0) && checkNoOverlap();
}

bool FloorPlan::isLegal(Module* curModule, ll attempCordX, ll attempCordY)
{
    for(int i = 0; i < vec_allModule.size(); i++)
    {
        Module* testModule = vec_allModule[i];

        // 檢查這個位置擺下去會不會和其他已經擺好的module重疊
        // 以下四個條件滿足其一即不會重疊
        if(testModule->isFixed && testModule != curModule)
        {
            if(attempCordX + curModule->width <= testModule->cordX ||       // 條件1 : xi + wi <= xj
               attempCordY + curModule->height <= testModule->cordY ||      // 條件2 : yi + hi <= yj
               attempCordX - testModule->width >= testModule->cordX ||      // 條件3 : xi - wj >= xj
               attempCordY - testModule->height >= testModule->cordY)      // 條件4 : yi - hj >= yj
               continue;
            else
                return false;
        }
    }

    return true;
}

bool FloorPlan::checkNoOverlap()
{
    for(int i = 0; i < vec_allModule.size(); i++)
    {
        Module* testModule = vec_allModule[i];

        if(!isLegal(testModule, testModule->cordX, testModule->cordY))
        {
            //cout << testModule->name << endl;
            //cout << "QQQQQQQQQQQppp" << endl;
        }
    }

    return true;
}

void FloorPlan::writeOutputFile(string outPath)
{
    ofstream outFile;
    outFile.open(outPath);

    //cout << "===================================" << endl;
    //cout << "output .floorplan file" << endl;
    //cout << endl;

    // calculate wirelength again
    long long weightHPWL = 0;

    for(int i = 0; i < vec_net.size(); i++)
    {
        Net* n = vec_net[i];

        Module* m1 = umap_nameTobestResModule[n->module1Name];
        Module* m2 = umap_nameTobestResModule[n->module2Name];

        ll weight = n->weight;

        ll center1X = (m1->cordX + (m1->cordX + m1->width))/2;
        ll center1Y = (m1->cordY + (m1->cordY + m1->height))/2;
        ll center2X = (m2->cordX + (m2->cordX + m2->width))/2;
        ll center2Y = (m2->cordY + (m2->cordY + m2->height))/2;

        weightHPWL += weight*(abs(center1X-center2X) + abs(center1Y-center2Y));
    }

    //cout << "weight HPWL: " << weightHPWL << endl;
    outFile << "Wirelength " << weightHPWL << endl;
    outFile << endl;

    // module cord
    outFile << "NumSoftModules " << vec_bestResult.size() << endl;
    for(int i = 0; i < vec_bestResult.size(); i++)
    {
        Module* curModule = vec_bestResult[i];

        outFile << curModule->name << " " << curModule->cordX << " " << curModule->cordY << " " << curModule->width << " " << curModule->height << endl;
    }

    //cout << "NumSoftModules " << vec_bestResult.size() << endl;
    for(int i = 0; i < vec_bestResult.size(); i++)
    {
        Module* curModule = vec_bestResult[i];

        //cout << curModule->name << " " << curModule->cordX << " " << curModule->cordY << " " << curModule->width << " " << curModule->height << endl;
    }

    outFile.close();

    auto end_time_parser = chrono::high_resolution_clock::now();
    auto duration_parser = chrono::duration_cast<chrono::milliseconds>(end_time_parser-start_time).count();
    //cout << "output file spend time: " << duration_parser << endl;   
    //cout << "=========================================" << endl;
}

void FloorPlan::printCurrentState()
{
    //cout << "===================" << endl;
    //cout << "current state" << endl << endl;

    //cout << "cost: " << calCost() << endl;

    //cout << "NumSoftModules " << vec_floatModule.size() << endl;
    for(int i = 0; i < vec_floatModule.size(); i++)
    {
        Module* curModule = vec_floatModule[i];

        //cout << "Module " << i << " :" << curModule->name << " " << curModule->cordX << " " << curModule->cordY << " " << curModule->width << " " << curModule->height << endl;
    }
    
    //cout << "Traverse B* tree: " << endl;
    queue<Node*> queue_dfs;
    queue_dfs.push(bTreeRoot);

    while(!queue_dfs.empty())
    {
        Node* topNode = queue_dfs.front();
        //cout << topNode->moduleIdx << " ";

        queue_dfs.pop();

        if(topNode->left)
            queue_dfs.push(topNode->left);
        
        if(topNode->right)
            queue_dfs.push(topNode->right);
    }
    //cout << endl;

    //cout << "===================" << endl << endl;

    return;
}

void FloorPlan::refine()
{
    //cout << "=======================" << endl;
    //cout << "LET'SSSSSSSSSSSSSSS go!!!" << endl << endl;

    copyfromBestinitial();

    // 先找出violate的module
    for(Module* m : vec_floatModule)
        m->isFixed = true;

    vector<Module*> vec_moduleNeededFixed;
    for(Module* m : vec_floatModule)
    {
        if(m->cordX < 0 || m->cordY < 0 || m->cordX + m->width > chipWidth || m->cordY+m->height > chipHeight)
            vec_moduleNeededFixed.push_back(m);
    }

    sort(vec_moduleNeededFixed.begin(), vec_moduleNeededFixed.end(), cmpDecreaseArea);
    
    for(Module* m : vec_moduleNeededFixed)
    {
        bool find = false;

        for(int i = 0; i < m->vec_allpossibleShape.size(); i++)
        {    
            m->width = m->vec_allpossibleShape[i].first;
            m->height = m->vec_allpossibleShape[i].second;

            for(ll attempX = 0; attempX + m->width <= chipWidth; attempX++)
            {
                for(ll attempY = 0; attempY + m->height <= chipHeight; attempY++)
                {
                    if(isLegal(m, attempX, attempY))
                    {
                        m->cordX = attempX;
                        m->cordY = attempY;

                        find = true;

                        break;
                    }
                }

                if(find)
                    break;
            }

            if(find)
                break;
        }

        // 修復不了就擺爛
        if(!find)
        {
            //cout << "Wow..This is a hard question..." << endl;
            return;
        }
    }

    updateBest();

    //cout << "     GOOOOOOOOOOOOD!!!" << endl;

    auto end_time_parser = chrono::high_resolution_clock::now();
    auto duration_parser = chrono::duration_cast<chrono::milliseconds>(end_time_parser-start_time).count();
    //cout << "refine Time: " << duration_parser << endl;
    //cout << "=======================" << endl << endl;

    return;
}
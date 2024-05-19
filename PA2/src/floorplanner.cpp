#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <climits>
#include <map>
#include <random>
#include <chrono>
#include "module.h"
#include "floorplanner.h"
using namespace std;


void Floorplanner::parseBlk(fstream& blkFile)
{
    string str;
    blkFile >> str;//outline
    blkFile >> str;
    _xmax = stod(str);
    blkFile >> str;
    _ymax = stod(str);
    blkFile >> str;//numblocks
    blkFile >> str;
    _blockNum = stod(str);
    blkFile >> str;//numterminals
    blkFile >> str;
    _terminalNum = stod(str);

    // Set up whole circuit

    for(int i = 0;i<_blockNum;i++) {
        blkFile >> str;
        string blkname = str;
        blkFile >> str;
        double blkw = stod(str);
        blkFile >> str;
        double blkh = stod(str);
        _blockArray.push_back(new Block(blkname,blkw,blkh));
        //_blockArray[i]->setMaxX(_xmax);
        //_blockArray[i]->setMaxY(_ymax);
        _blockName2Id[blkname] = i;
    }
    for(int i = 0;i<_terminalNum;i++) {
        blkFile >> str;
        //cout <<str<<endl;
        string tername = str;
        blkFile >> str;
        blkFile >> str;
        double terx = stod(str);
        blkFile >> str;
        double tery = stod(str);
        _terminalArray.push_back(new Terminal(tername,terx,tery));
        _terminalName2Id[tername] = i;
    }
    
    return;
}


void Floorplanner::parseNet(fstream& netFile)
{
    string str;
    netFile >> str;
    netFile >> str;
    _netNum = stod(str);
    int netdegree;
    for(int i = 0;i<_netNum;i++){
        netFile >> str;
        netFile >> str;
        netdegree = stod(str);
        _netArray.push_back(new Net());
        for(int j = 0;j<netdegree;j++){
            netFile >> str;
            string tername = str;
            if(_terminalName2Id.count(tername)==0){
                _netArray[i]->addTerm(_blockArray[_blockName2Id[tername]]);
            }
            else{
                _netArray[i]->addTerm(_terminalArray[_terminalName2Id[tername]]);
            }
        }
    }
    
}

void Floorplanner::recur(Node* root)
{
    double y = 0;
    double x;
    //if encounter null node then return
    if(!root){
        return;
    }

    //the block corresponding to the node
    Block* rootBlk = _blockArray[_blockName2Id[root->getName()]];

    //if the node is not the root of the b tree
    if(root->getParent()){
        //the parent block
        Block* parentBlk = _blockArray[_blockName2Id[root->getParent()->getName()]];
        //if the node is the left child of its parent
        if(root->getParent()->getLeft() == root){
            x = parentBlk->getX2();
        }
        else{
            x = parentBlk->getX1();
        }
        int count = 0;
        for(int i = 0;i<doublyLinkedlist.size();i++){
            if(doublyLinkedlist[i][0] == x){
                //if there is value in doubly linked list encountering the x of block
                count++;
            }
            else if(count>0){
                //retrieve the latest y of the doubly linked list with the same x value
                y = doublyLinkedlist[i-1][1];
                count = 0;
            }
            if(doublyLinkedlist[i][0]>x && doublyLinkedlist[i][0]<(x+rootBlk->getWidth())){
                if(doublyLinkedlist[i][1]>y){
                    //get max value for y
                    y = doublyLinkedlist[i][1];
                }
            }
            if(doublyLinkedlist[i][0] == (x+rootBlk->getWidth())){
                if(doublyLinkedlist[i][1]>y){
                    //retreive the first y of the doubly linked list with x + w
                    y = doublyLinkedlist[i][1];
                }
                break;
            }
        }
        rootBlk->setPos(x,y,x+rootBlk->getWidth(),y+rootBlk->getHeight());
        vector<vector<double>> temp = {{rootBlk->getX1(),rootBlk->getY2()},{rootBlk->getX2(),rootBlk->getY2()},{rootBlk->getX2(),rootBlk->getY1()}};
        vector<vector<double>> bef = {};
        vector<vector<double>> aft = {};
        count = 0;
        //changing the doubly linked list
        for(int i = 0;i<doublyLinkedlist.size();i++){
            if(doublyLinkedlist[i][0]<temp[0][0]){
                bef.push_back(doublyLinkedlist[i]);
            }
            //only push the first element encountered because the latter elements are redundant
            if(doublyLinkedlist[i][0]==temp[0][0]){
                bef.push_back(doublyLinkedlist[i]);
                break;
            }
        }
        for(int i = 0;i<doublyLinkedlist.size();i++){
            if(doublyLinkedlist[i][0]>temp[2][0]){
                aft.push_back(doublyLinkedlist[i]);
            }
            //push back the redundant terms anyway because only the first will be used 
            //because of the rule of assignment of y(when encounter x+w assign and break instantly)
            if(doublyLinkedlist[i][0]==temp[2][0]){
                aft.push_back(doublyLinkedlist[i]);
            }
        }
        doublyLinkedlist.clear();
        doublyLinkedlist = bef;
        doublyLinkedlist.insert(doublyLinkedlist.end(),temp.begin(),temp.end());
        doublyLinkedlist.insert(doublyLinkedlist.end(),aft.begin(),aft.end());
    }
    else{
        //if the node is the root of the b tree
        vector<double> temp1 = {0,0};
        vector<double> temp2 = {INT_MAX,0};
        doublyLinkedlist.push_back(temp1);
        doublyLinkedlist.push_back(temp2);
        rootBlk->setPos(0,0,rootBlk->getWidth(),rootBlk->getHeight());
        vector<vector<double>> temp = {{rootBlk->getX1(),rootBlk->getY2()},{rootBlk->getX2(),rootBlk->getY2()},{rootBlk->getX2(),rootBlk->getY1()}};
        vector<vector<double>> bef = {};
        vector<vector<double>> aft = {};
        int count = 0;
        for(int i = 0;i<doublyLinkedlist.size();i++){
            if(doublyLinkedlist[i][0]<temp[0][0]){
                bef.push_back(doublyLinkedlist[i]);
            }
            if(doublyLinkedlist[i][0]==temp[0][0]){
                bef.push_back(doublyLinkedlist[i]);
                break;
            }
        }
        for(int i = 0;i<doublyLinkedlist.size();i++){
            if(doublyLinkedlist[i][0]>temp[2][0]){
                aft.push_back(doublyLinkedlist[i]);
            }
            if(doublyLinkedlist[i][0]==temp[2][0]){
                aft.push_back(doublyLinkedlist[i]);
            }
        }
        doublyLinkedlist.clear();
        doublyLinkedlist = bef;
        doublyLinkedlist.insert(doublyLinkedlist.end(),temp.begin(),temp.end());
        doublyLinkedlist.insert(doublyLinkedlist.end(),aft.begin(),aft.end());
    }
    //performing dfs on the b tree
    recur(root->getLeft());
    recur(root->getRight());
}

void Floorplanner::treetograph(Node* root)
{
    //for reconstructing the doublylinkedlist
    doublyLinkedlist.clear();
    recur(root);
}

//replaced by copytree
void Floorplanner::copyrecur(Node* &node, Node* &bestNode)
{
    if(!node){
        return;
    }
    bestNode = _blockArray[_blockName2Id[node->getName()]]->getBestNode();
    bestNode->setLeft(NULL);
    bestNode->setRight(NULL);
    bestNode->setParent(NULL);
    //cout << bestNode->getName() << endl;
    Node* parentnode = node->getParent();
    Node* leftnode = node->getLeft();
    Node* rightnode = node->getRight();

    if(parentnode){
        bestNode->setParent(_blockArray[_blockName2Id[parentnode->getName()]]->getBestNode());
    }
    if(leftnode){
        bestNode->setLeft(_blockArray[_blockName2Id[leftnode->getName()]]->getBestNode());
    }
    if(rightnode){
        bestNode->setRight(_blockArray[_blockName2Id[rightnode->getName()]]->getBestNode());
    }
    Node* bestleft;
    Node* bestright;
    copyrecur(leftnode,bestleft);
    copyrecur(rightnode,bestright);
}

//if true copy the current tree to the best tree
//else copy the best tree to the current tree
void Floorplanner::copytree(bool flag)
{
    Node* node;
    Node* bestNode;
    if(flag){
        _bestrotatelist.clear();
        _bestRoot = _blockArray[_blockName2Id[_root->getName()]]->getBestNode();
        for(int i = 0;i<_blockArray.size();i++){
            //recording the best rotate state for each block
            _bestrotatelist.push_back(_blockArray[i]->getRotate());
            node = _blockArray[i]->getNode();
            bestNode = _blockArray[i]->getBestNode();
            //clear the relationship for every best node
            bestNode->setLeft(NULL);
            bestNode->setRight(NULL);
            bestNode->setParent(NULL);
            Node* parentnode = node->getParent();
            Node* leftnode = node->getLeft();
            Node* rightnode = node->getRight();
            if(parentnode){
                bestNode->setParent(_blockArray[_blockName2Id[parentnode->getName()]]->getBestNode());
            }
            if(leftnode){
                bestNode->setLeft(_blockArray[_blockName2Id[leftnode->getName()]]->getBestNode());
            }
            if(rightnode){
                bestNode->setRight(_blockArray[_blockName2Id[rightnode->getName()]]->getBestNode());
            }
        }
    }
    else{
        for(int i = 0;i<_bestrotatelist.size();i++){
            _blockArray[i]->setRotate(_bestrotatelist[i]);
        }
        _root = _blockArray[_blockName2Id[_bestRoot->getName()]]->getNode();
        for(int i = 0;i<_blockArray.size();i++){
            node = _blockArray[i]->getNode();
            bestNode = _blockArray[i]->getBestNode();
            node->setLeft(NULL);
            node->setRight(NULL);
            node->setParent(NULL);
            Node* parentnode = bestNode->getParent();
            Node* leftnode = bestNode->getLeft();
            Node* rightnode = bestNode->getRight();
            if(parentnode){
                node->setParent(_blockArray[_blockName2Id[parentnode->getName()]]->getNode());
            }
            if(leftnode){
                node->setLeft(_blockArray[_blockName2Id[leftnode->getName()]]->getNode());
            }
            if(rightnode){
                node->setRight(_blockArray[_blockName2Id[rightnode->getName()]]->getNode());
            }
        }
    }
}

//directly changing the name blk pointer to swap the node
void Floorplanner::swapnode(Node* &node1, Node* &node2){
    string temp1 = node1->getName();
    string temp2 = node2->getName();
    Block* blk1 = _blockArray[_blockName2Id[node1->getName()]];
    Block* blk2 = _blockArray[_blockName2Id[node2->getName()]];
    node1->setName(temp2);
    node2->setName(temp1);
    blk1->setNode(node2);
    blk2->setNode(node1);

    return;
}

//delete a leaf node and insert it to external position
void Floorplanner::deleteandinsert(Node* &node1, Node* &node2)
{
    Node* node1p = node1->getParent();
    if(node1p->getLeft()){
        if(_blockName2Id[node1p->getLeft()->getName()]==_blockName2Id[node1->getName()]){
            //if the node is the left child of its parent clear the left child of its parent
            node1p->setLeft(NULL);
        }
        else{
            node1p->setRight(NULL);
        }
    }
    else{
        node1p->setRight(NULL);
    }
    //reset node1's parent to node2
    node1->setParent(node2);
    //if node2 has leftchild then insert node1 to its right position
    if(node2->getLeft()){
        node2->setRight(node1);
    }
    else if(node2->getRight()){
        node2->setLeft(node1);
    }
    else{
        //if node2 has no child then insert node1 to either position with half the chance
        int a = rand();
        if(a %2==1){
            node2->setRight(node1);
        }
        else{
            node2->setLeft(node1);
        }
    }    
}

void Floorplanner::swapsubtree(Node* &node){
    Node* left = node->getLeft();
    Node* right = node->getRight();
    node->setLeft(right);
    node->setRight(left);
}

//main function to plan
void Floorplanner::floorplan()
{
    //start the timer
    auto start = std::chrono::high_resolution_clock::now();
    
    //construct initial b tree complete binary tree
    int parent;
    _HPWL = 0;
    _chipHeight = 0;
    _chipWidth = 0;
    _chipArea = 0;
    _root = _blockArray[0]->getNode();
    
    for(int i = 1;i<_blockArray.size();i++){
        parent = (i-1)/2;
        if(i%2==1){
            _blockArray[i]->getNode()->setParent(_blockArray[parent]->getNode());
            _blockArray[parent]->getNode()->setLeft(_blockArray[i]->getNode());
        }
        else{
            _blockArray[i]->getNode()->setParent(_blockArray[parent]->getNode());
            _blockArray[parent]->getNode()->setRight(_blockArray[i]->getNode());
        }
    }
    //convert intial b tree to graph
    treetograph(_root);
    //copy the current tree to the best tree
    copytree(true);
    //calculate HPWL
    for(int i = 0;i<_netArray.size();i++){
        _HPWL += _netArray[i]->calcHPWL();
    }
    //calculate initial chip width and height
    for(int i = 0;i<_blockArray.size();i++){
        if(_blockArray[i]->getX2()>_chipWidth){
            _chipWidth = _blockArray[i]->getX2();
        }
        if(_blockArray[i]->getY2()>_chipHeight){
            _chipHeight = _blockArray[i]->getY2();
        }
    }
    //calculate initial chip area
    _chipArea = _chipHeight*_chipWidth;

    //for recording the iteration
    int iter = 0;

    //set initial _bestchipheight to infinity
    _bestChipHeight = INT_MAX;
    _bestChipWidth = INT_MAX;
    _bestChipArea = INT_MAX;
    _bestdiff = INT_MAX;
    _bestcost = INT_MAX;
    //if best chip postion fit the outline then break
    //trying to find legal initial position
    std::random_device ma;
    std::mt19937 mt(ma());
    
    while(_bestChipHeight>_ymax || _bestChipWidth>_xmax){
        iter++;
        //if iteration hit the limit than break
        // if(iter > 1000000){
        //     break;
        // }
        if(iter%15000==0){
            _HPWL = 0;
            _chipHeight = 0;
            _chipWidth = 0;
            _chipArea = 0;
            _bestChipHeight = INT_MAX;
            _bestChipWidth = INT_MAX;
            _bestChipArea = INT_MAX;
            _bestdiff = INT_MAX;
            _bestcost = INT_MAX;
            int temp;
            //randomize a whole new b tree
            vector<Block*> tempBlockArray = _blockArray;
            vector<Block*> newblockArray = {};

            while(tempBlockArray.size()>0){
                std::uniform_int_distribution<int> rg(0,tempBlockArray.size()-1);
                temp = rg(mt);
                newblockArray.push_back(tempBlockArray[temp]);
                tempBlockArray.erase(tempBlockArray.begin()+temp);
            }
            _root = newblockArray[0]->getNode();
            for(int i = 0;i<newblockArray.size();i++){
                newblockArray[i]->getNode()->setParent(NULL);
                newblockArray[i]->getNode()->setLeft(NULL);
                newblockArray[i]->getNode()->setRight(NULL);
            }
            for(int i = 1;i<newblockArray.size();i++){
                parent = (i-1)/2;
                if(i%2==1){
                    newblockArray[i]->getNode()->setParent(newblockArray[parent]->getNode());
                    newblockArray[parent]->getNode()->setLeft(newblockArray[i]->getNode());
                }
                else{
                    newblockArray[i]->getNode()->setParent(newblockArray[parent]->getNode());
                    newblockArray[parent]->getNode()->setRight(newblockArray[i]->getNode());
                }
            }
            //copy the current tree to the best tree
            copytree(true);
            
            //convert intial b tree to graph
            treetograph(_root);
        }
        else{
            //reset the chip parameter to 0 for recalculation
            //copy the best b tree to current tree for optimization
            copytree(false);
        }
        _chipHeight = 0;
        _chipWidth = 0;
        _chipArea = 0;
        _HPWL = 0;

        //randomly choose 5 node to delete and insert, swap or rotate
        std::uniform_int_distribution<int> rg(0,_blockArray.size()-1);
        int rn1 = rg(mt);
        int rn2 = rg(mt);
        int rn3 = rg(mt);
        int rn4 = rg(mt);
        int rn5 = rg(mt);
        int rn6 = rg(mt);
        //for delete and insert
        Node* node1 = _blockArray[rn1]->getNode();
        Node* node2 = _blockArray[rn2]->getNode();

        //for swap
        Node* node3 = _blockArray[rn3]->getNode();
        Node* node4 = _blockArray[rn4]->getNode();

        Node* node6 = _blockArray[rn6]->getNode();

        //randomly choose one operation
        std::uniform_int_distribution<int> rg3(0,3);
        int choice = rg3(mt);
        if(choice == 1){
            //node3 must not equal node4
            while(rn3 == rn4){
                rn4 = rg(mt);
                node4 = _blockArray[rn4]->getNode();
            }
            swapnode(node3,node4);
        }
        else if(choice == 2){
            //node 1 must be a leafnode
            while(node1->getLeft() || node1->getRight()){
                rn1 = rg(mt);
                node1 = _blockArray[rn1]->getNode();
            }
            //node2 must have null child and must not equal node1
            while((node2->getLeft() && node2->getRight())||rn1==rn2){
                rn2 = rg(mt);
                node2 = _blockArray[rn2]->getNode();
            }
            deleteandinsert(node1,node2);
        }
        else if(choice == 3){
            _blockArray[rn5]->rotate();
        }
        else{
            while(!node6->getLeft() && !node6->getRight()){
                rn6 = rg(mt);
                node6 = _blockArray[rn6]->getNode();
            }
            swapsubtree(node6);
        }

        //map tree to graph
        treetograph(_root);

        //calculate the parameter again
        for(int i = 0;i<_netArray.size();i++){
            _HPWL += _netArray[i]->calcHPWL();
        }
        for(int i = 0;i<_blockArray.size();i++){
            if(_blockArray[i]->getX2()>_chipWidth){
                _chipWidth = _blockArray[i]->getX2();
            }
            if(_blockArray[i]->getY2()>_chipHeight){
                _chipHeight = _blockArray[i]->getY2();
            }
        }
        _cost = _alpha*_chipHeight*_chipWidth+(1-_alpha)*_HPWL;
        _chipArea = _chipHeight*_chipWidth;
        int xdif = _chipWidth-_xmax;
        int ydif = _chipHeight-_ymax;
        // if(xdif<0){
        //     xdif = -xdif;
        // }
        // if(ydif<0){
        //     ydif = -ydif;
        // }
        int totaldif = xdif+ydif;
        if(_chipArea<_bestChipArea){
            //if the cost is better than the best aka the orginal tree
            //than copy the tree to the best tree
            copytree(true);

            //record the parameter
            _bestdiff = totaldif;
            _bestcost = _cost;
            _bestChipHeight = _chipHeight;
            _bestChipWidth = _chipWidth;
            _bestChipArea = _chipArea;
        }
        // cout << "iteration:  " << iter << endl;
        // cout.precision(10);
        // cout << "bestcost:   " << _bestcost << endl;
        // cout << "bestwidth:  "<< _bestChipWidth << " " << _xmax << endl;
        // cout << "bestheight: " << _bestChipHeight << " " << _ymax << endl;
    }

    //set the rotate state to best rotate state for each block
    for(int i = 0;i<_bestrotatelist.size();i++){
        _blockArray[i]->setRotate(_bestrotatelist[i]);
    }
    //map best b tree to graph
    treetograph(_bestRoot);

    //recalculate the parameters
    _HPWL = 0;
    _chipHeight = 0;
    _chipWidth = 0;
    _chipArea = 0;
    for(int i = 0;i<_netArray.size();i++){
        _HPWL += (double)(_netArray[i]->calcHPWL());
    }
    for(int i = 0;i<_blockArray.size();i++){
        if(_blockArray[i]->getX2()>_chipWidth){
            _chipWidth = _blockArray[i]->getX2();
        }
        if(_blockArray[i]->getY2()>_chipHeight){
            _chipHeight = _blockArray[i]->getY2();
        }
    }
    _chipArea = _chipHeight*_chipWidth;
    _cost = (_alpha*_chipHeight*_chipWidth);
    _cost += (1-_alpha)*_HPWL;
    _AreaNorm = _chipArea;
    _WireNorm = _HPWL;
    double _bestCostNorm = (_alpha*_chipArea/_AreaNorm) + (1-_alpha)*_HPWL/_WireNorm;

    //start Random Search

    //initial solution
    copytree(false);
    double T = 10000;
    int P = 1000;
    int limi = 10000000;
    if(iter>2000000){
        limi = 0;
    }
    iter=0;
    double DE;
    double cur_cost;
    double nei_cost;
    
    while(iter<limi){
        iter++;
        if(iter%15000==0){
            _HPWL = 0;
            _chipHeight = 0;
            _chipWidth = 0;
            _chipArea = 0;
            int temp;
            //randomize a whole new b tree
            vector<Block*> tempBlockArray = _blockArray;
            vector<Block*> newblockArray = {};

            while(tempBlockArray.size()>0){
                std::uniform_int_distribution<int> rg(0,tempBlockArray.size()-1);
                temp = rg(mt);
                newblockArray.push_back(tempBlockArray[temp]);
                tempBlockArray.erase(tempBlockArray.begin()+temp);
            }
            _root = newblockArray[0]->getNode();
            for(int i = 0;i<newblockArray.size();i++){
                newblockArray[i]->getNode()->setParent(NULL);
                newblockArray[i]->getNode()->setLeft(NULL);
                newblockArray[i]->getNode()->setRight(NULL);
            }
            for(int i = 1;i<newblockArray.size();i++){
                parent = (i-1)/2;
                if(i%2==1){
                    newblockArray[i]->getNode()->setParent(newblockArray[parent]->getNode());
                    newblockArray[parent]->getNode()->setLeft(newblockArray[i]->getNode());
                }
                else{
                    newblockArray[i]->getNode()->setParent(newblockArray[parent]->getNode());
                    newblockArray[parent]->getNode()->setRight(newblockArray[i]->getNode());
                }
            }
        }
        _chipHeight = 0;
        _chipWidth = 0;
        _chipArea = 0;
        _HPWL = 0;
        //map tree to graph
        treetograph(_root);
        //calculate the parameter again
        for(int i = 0;i<_netArray.size();i++){
            _HPWL += _netArray[i]->calcHPWL();
        }
        for(int i = 0;i<_blockArray.size();i++){
            if(_blockArray[i]->getX2()>_chipWidth){
                _chipWidth = _blockArray[i]->getX2();
            }
            if(_blockArray[i]->getY2()>_chipHeight){
                _chipHeight = _blockArray[i]->getY2();
            }
        }
        _cost = _alpha*_chipHeight*_chipWidth+(1-_alpha)*_HPWL;
        _chipArea = _chipHeight*_chipWidth;
        cur_cost = (_alpha*_chipArea/_AreaNorm) + (1-_alpha)*_HPWL/_WireNorm;

        //randomly choose 5 node to delete and insert, swap or rotate
        std::random_device ma;
        std::mt19937 mt(ma());
        std::uniform_int_distribution<int> rg(0,_blockArray.size()-1);
        int rn1 = rg(mt);
        int rn2 = rg(mt);
        int rn3 = rg(mt);
        int rn4 = rg(mt);
        int rn5 = rg(mt);
        int rn6 = rg(mt);
        //for delete and insert
        Node* node1 = _blockArray[rn1]->getNode();
        Node* node2 = _blockArray[rn2]->getNode();

        //for swap
        Node* node3 = _blockArray[rn3]->getNode();
        Node* node4 = _blockArray[rn4]->getNode();

        Node* node6 = _blockArray[rn6]->getNode();

        Node* tempParent;
        //randomly choose one operation
        std::uniform_int_distribution<int> rg3(0,3);
        int choice = rg3(mt);
        if(choice == 1){
            //node3 must not equal node4
            while(rn3 == rn4){
                rn4 = rg(mt);
                node4 = _blockArray[rn4]->getNode();
            }
            swapnode(node3,node4);
        }
        else if(choice == 2){
            //node 1 must be a leafnode
            while(node1->getLeft() || node1->getRight()){
                rn1 = rg(mt);
                node1 = _blockArray[rn1]->getNode();
            }
            //node2 must have null child and must not equal node1
            while((node2->getLeft() && node2->getRight())||rn1==rn2){
                rn2 = rg(mt);
                node2 = _blockArray[rn2]->getNode();
            }
            tempParent = node1->getParent();
            deleteandinsert(node1,node2);
        }
        else if(choice == 3){
            _blockArray[rn5]->rotate();
        }
        else{
            while(!node6->getLeft() && !node6->getRight()){
                rn6 = rg(mt);
                node6 = _blockArray[rn6]->getNode();
            }
            swapsubtree(node6);
        }
        _chipHeight = 0;
        _chipWidth = 0;
        _chipArea = 0;
        _HPWL = 0;
        //map tree to graph
        treetograph(_root);
        //calculate the parameter again
        for(int i = 0;i<_netArray.size();i++){
            _HPWL += _netArray[i]->calcHPWL();
        }
        for(int i = 0;i<_blockArray.size();i++){
            if(_blockArray[i]->getX2()>_chipWidth){
                _chipWidth = _blockArray[i]->getX2();
            }
            if(_blockArray[i]->getY2()>_chipHeight){
                _chipHeight = _blockArray[i]->getY2();
            }
        }
        _cost = _alpha*_chipHeight*_chipWidth+(1-_alpha)*_HPWL;
        _chipArea = _chipHeight*_chipWidth;
        nei_cost = (_alpha*_chipArea/_AreaNorm) + (1-_alpha)*_HPWL/_WireNorm;
        DE = nei_cost - cur_cost;

        //uphill
        if(DE>0){
            if(choice == 1){
                swapnode(node3,node4);
            }
            else if(choice == 2){
                deleteandinsert(node1,tempParent);
            }
            else if(choice == 3){
                _blockArray[rn5]->rotate();
            }
            else{
                swapsubtree(node6);
            }
        }

        _chipHeight = 0;
        _chipWidth = 0;
        _chipArea = 0;
        _HPWL = 0;
        treetograph(_root);
        //calculate the parameter again
        for(int i = 0;i<_netArray.size();i++){
            _HPWL += _netArray[i]->calcHPWL();
        }
        for(int i = 0;i<_blockArray.size();i++){
            if(_blockArray[i]->getX2()>_chipWidth){
                _chipWidth = _blockArray[i]->getX2();
            }
            if(_blockArray[i]->getY2()>_chipHeight){
                _chipHeight = _blockArray[i]->getY2();
            }
        }
        _cost = _alpha*_chipHeight*_chipWidth+(1-_alpha)*_HPWL;
        _chipArea = _chipHeight*_chipWidth;

        cur_cost = (_alpha*_chipArea/_AreaNorm) + (1-_alpha)*_HPWL/_WireNorm;
        if(cur_cost<_bestCostNorm && _chipHeight<=_ymax && _chipWidth <= _xmax){
            //if the cost is better than the best aka the orginal tree
            //than copy the tree to the best tree
            copytree(true);

            //record the parameter
            _bestCostNorm = cur_cost;
            _bestcost = _cost;
            _bestChipHeight = _chipHeight;
            _bestChipWidth = _chipWidth;
            _bestChipArea = _chipArea;
        }
        // cout << "iteration:  " << iter << endl;
        // cout.precision(10);
        // cout << "normcost:   " << _bestCostNorm << endl;
        // cout << "bestcost:   " << _bestcost << endl;
        // cout << "bestwidth:  "<< _bestChipWidth << " " << _xmax << endl;
        // cout << "bestheight: " << _bestChipHeight << " " << _ymax << endl;
    }
    for(int i = 0;i<_bestrotatelist.size();i++){
        _blockArray[i]->setRotate(_bestrotatelist[i]);
    }
    treetograph(_bestRoot);
    _HPWL = 0;
    _chipHeight = 0;
    _chipWidth = 0;
    _chipArea = 0;
    for(int i = 0;i<_netArray.size();i++){
        _HPWL += (double)(_netArray[i]->calcHPWL());
    }
    for(int i = 0;i<_blockArray.size();i++){
        if(_blockArray[i]->getX2()>_chipWidth){
            _chipWidth = _blockArray[i]->getX2();
        }
        if(_blockArray[i]->getY2()>_chipHeight){
            _chipHeight = _blockArray[i]->getY2();
        }
    }
    _chipArea = _chipHeight*_chipWidth;
    _cost = (_alpha*_chipHeight*_chipWidth);
    _cost += (1-_alpha)*_HPWL;

    // //start SA
    // _WireNorm = _HPWL;
    // _AreaNorm = _chipArea;
    // _bestCostFunc = getcostFunc();
    // //initial solution
    // copytree(false);
    // double T = 10000000;
    // int P = 1000;
    // int limi = 10000000;
    // if(iter>10000000){
    //     limi =  0;
    // }
    // iter=0;
    // double DE;
    // double cur_cost;
    // double nei_cost;
    // while(T>0.01 && iter<limi){
    //     for(int i = 0;i<P;i++){
    //         iter++;
    //         cur_cost = getcostFunc();

    //         //randomly choose 5 node to delete and insert, swap or rotate
    //         std::random_device ma;
    //         std::mt19937 mt(ma());
    //         std::uniform_int_distribution<int> rg(0,_blockArray.size()-1);
    //         int rn1 = rg(mt);
    //         int rn2 = rg(mt);
    //         int rn3 = rg(mt);
    //         int rn4 = rg(mt);
    //         int rn5 = rg(mt);
    //         int rn6 = rg(mt);
    //         //for delete and insert
    //         Node* node1 = _blockArray[rn1]->getNode();
    //         Node* node2 = _blockArray[rn2]->getNode();

    //         //for swap
    //         Node* node3 = _blockArray[rn3]->getNode();
    //         Node* node4 = _blockArray[rn4]->getNode();

    //         Node* node6 = _blockArray[rn6]->getNode();

    //         Node* tempParent;
    //         //randomly choose one operation
    //         std::uniform_int_distribution<int> rg3(0,3);
    //         int choice = rg3(mt);
    //         if(choice == 1){
    //             //node3 must not equal node4
    //             while(rn3 == rn4){
    //                 rn4 = rg(mt);
    //                 node4 = _blockArray[rn4]->getNode();
    //             }
    //             swapnode(node3,node4);
    //         }
    //         else if(choice == 2){
    //             //node 1 must be a leafnode
    //             while(node1->getLeft() || node1->getRight()){
    //                 rn1 = rg(mt);
    //                 node1 = _blockArray[rn1]->getNode();
    //             }
    //             //node2 must have null child and must not equal node1
    //             while((node2->getLeft() && node2->getRight())||rn1==rn2){
    //                 rn2 = rg(mt);
    //                 node2 = _blockArray[rn2]->getNode();
    //             }
    //             tempParent = node1->getParent();
    //             deleteandinsert(node1,node2);
    //         }
    //         else if(choice == 3){
    //             _blockArray[rn5]->rotate();
    //         }
    //         else{
    //             while(!node6->getLeft() && !node6->getRight()){
    //                 rn6 = rg(mt);
    //                 node6 = _blockArray[rn6]->getNode();
    //             }
    //             swapsubtree(node6);
    //         }
    //         _chipHeight = 0;
    //         _chipWidth = 0;
    //         _chipArea = 0;
    //         _HPWL = 0;
    //         //map tree to graph
    //         treetograph(_root);
    //         //calculate the parameter again
    //         for(int i = 0;i<_netArray.size();i++){
    //             _HPWL += _netArray[i]->calcHPWL();
    //         }
    //         for(int i = 0;i<_blockArray.size();i++){
    //             if(_blockArray[i]->getX2()>_chipWidth){
    //                 _chipWidth = _blockArray[i]->getX2();
    //             }
    //             if(_blockArray[i]->getY2()>_chipHeight){
    //                 _chipHeight = _blockArray[i]->getY2();
    //             }
    //         }
    //         _cost = _alpha*_chipHeight*_chipWidth+(1-_alpha)*_HPWL;
    //         _chipArea = _chipHeight*_chipWidth;
    //         nei_cost = getcostFunc();
    //         DE = nei_cost - cur_cost;

    //         std::uniform_real_distribution<double> pb(0,1);
    //         int prob = pb(mt);
    //         //cout <<"de: "<<DE<<endl;
    //         //cout <<"prob: "<< exp(-DE/T)<<endl;
    //         //cout <<"T: " << T << endl;
    //         //uphill
    //         if(DE>0){
    //             if(prob>(double)exp(-DE/T)){
    //                 if(choice == 1){
    //                     swapnode(node3,node4);
    //                 }
    //                 else if(choice == 2){
    //                     deleteandinsert(node1,tempParent);
    //                 }
    //                 else if(choice == 3){
    //                     _blockArray[rn5]->rotate();
    //                 }
    //                 else{
    //                     swapsubtree(node6);
    //                 }
    //             }
    //         }

    //         _chipHeight = 0;
    //         _chipWidth = 0;
    //         _chipArea = 0;
    //         _HPWL = 0;
    //         treetograph(_root);
    //         //calculate the parameter again
    //         for(int i = 0;i<_netArray.size();i++){
    //             _HPWL += _netArray[i]->calcHPWL();
    //         }
    //         for(int i = 0;i<_blockArray.size();i++){
    //             if(_blockArray[i]->getX2()>_chipWidth){
    //                 _chipWidth = _blockArray[i]->getX2();
    //             }
    //             if(_blockArray[i]->getY2()>_chipHeight){
    //                 _chipHeight = _blockArray[i]->getY2();
    //             }
    //         }
    //         _cost = _alpha*_chipHeight*_chipWidth+(1-_alpha)*_HPWL;
    //         _chipArea = _chipHeight*_chipWidth;

    //         cur_cost = getcostFunc();
    //         if(cur_cost<_bestCostFunc && _chipWidth<=_xmax && _chipHeight<=_ymax){
    //             //if the cost is better than the best aka the orginal tree
    //             //than copy the tree to the best tree
    //             copytree(true);

    //             //record the parameter
    //             _bestCostFunc = cur_cost;
    //             _bestcost = _cost;
    //             _bestChipHeight = _chipHeight;
    //             _bestChipWidth = _chipWidth;
    //             _bestChipArea = _chipArea;
    //         }
    //         cout << "iteration:  " << iter << endl;
    //         cout.precision(10);
    //         cout << "bestcost:   " << _bestcost << endl;
    //         cout << "bestwidth:  "<< _bestChipWidth << " " << _xmax << endl;
    //         cout << "bestheight: " << _bestChipHeight << " " << _ymax << endl;
    //     }
    //     T = 0.95*T;
    // }
    // for(int i = 0;i<_bestrotatelist.size();i++){
    //     _blockArray[i]->setRotate(_bestrotatelist[i]);
    // }
    // treetograph(_bestRoot);
    // _HPWL = 0;
    // _chipHeight = 0;
    // _chipWidth = 0;
    // _chipArea = 0;
    // for(int i = 0;i<_netArray.size();i++){
    //     _HPWL += (double)(_netArray[i]->calcHPWL());
    // }
    // for(int i = 0;i<_blockArray.size();i++){
    //     if(_blockArray[i]->getX2()>_chipWidth){
    //         _chipWidth = _blockArray[i]->getX2();
    //     }
    //     if(_blockArray[i]->getY2()>_chipHeight){
    //         _chipHeight = _blockArray[i]->getY2();
    //     }
    // }
    // _chipArea = _chipHeight*_chipWidth;
    // _cost = (_alpha*_chipHeight*_chipWidth);
    // _cost += (1-_alpha)*_HPWL;

    //end timer
    auto end = std::chrono::high_resolution_clock::now();
    //calculate running time
    _runningTime = (double)std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()/1000;
}

void Floorplanner::printSummary() const
{
    cout.precision(10);
    cout << "Cost: " << _cost << endl;
    cout.precision(10);
    cout << "HPWL: " << _HPWL << endl;
    cout << "chipwidth: " << _chipWidth << endl;
    cout << "chipheight: " << _chipHeight << endl;
    cout.precision(10);
    cout << "chiparea: " << _chipArea << endl;
    for (size_t i = 0, end_i = _blockArray.size(); i < end_i; ++i) {
        cout << _blockArray[i]->getName() << ": ";
        cout << _blockArray[i]->getX1() << " " << _blockArray[i]->getY1() << " " << _blockArray[i]->getX2() << " " << _blockArray[i]->getY2();
        //cout << setw(8) << _blockArray[i]->getWidth() << " " << _blockArray[i]->getHeight();
        cout << endl;
    }
    return;
}

void Floorplanner::reportBlk() const
{
    cout << "Number of blocks: " << _blockNum << endl;
    for (size_t i = 0, end_i = _blockArray.size(); i < end_i; ++i) {
        cout << setw(8) << _blockArray[i]->getName() << ": ";
        cout << setw(8) << _blockArray[i]->getX1() << " " << _blockArray[i]->getY1() << " " << _blockArray[i]->getX2() << " " << _blockArray[i]->getY2();
        //cout << setw(8) << _blockArray[i]->getWidth() << " " << _blockArray[i]->getHeight();
        cout << endl;
    }
    return;
}

void Floorplanner::reportTerm() const
{
    cout << "Number of terminals: " << _terminalNum << endl;
    for (size_t i = 0, end_i = _terminalArray.size(); i < end_i; ++i) {
        cout << setw(8) << _terminalArray[i]->getName() << ": ";
        cout << setw(8) << _terminalArray[i]->getX1() << " " << _terminalArray[i]->getY1();
        cout << endl;
    }
    return;
}

void Floorplanner::reportNet() const
{
    cout << "Number of nets: " << _netNum << endl;
    for (size_t i = 0, end_i = _netArray.size(); i < end_i; ++i) {
        cout << setw(8) << i << ": ";
        vector<Terminal*> termList = _netArray[i]->getTermList();
        for (size_t j = 0, end_j = termList.size(); j < end_j; ++j) {
            cout << setw(8) << termList[j]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Floorplanner::writeResult(fstream& outFile)
{
    stringstream buff;
    buff.precision(10);
    buff << _cost;
    outFile << buff.str() << '\n';
    buff.str("");
    buff.precision(10);
    buff << _HPWL;
    outFile << buff.str() << '\n';

    buff.str("");
    buff.precision(10);
    buff << _chipArea;
    outFile << buff.str() << '\n';
    buff.str("");
    buff << _chipWidth;
    outFile << buff.str() << ' ';
    buff.str("");
    buff << _chipHeight;
    outFile << buff.str() << '\n';
    buff.str("");
    buff << _runningTime;
    outFile << buff.str() << '\n';
    for (size_t i = 0, end_i = _blockArray.size(); i < end_i; ++i) {
        outFile << _blockArray[i]->getName() << " ";
        outFile << _blockArray[i]->getX1() << " " << _blockArray[i]->getY1() << " " << _blockArray[i]->getX2() << " " << _blockArray[i]->getY2();
        //cout << setw(8) << _blockArray[i]->getWidth() << " " << _blockArray[i]->getHeight();
        outFile << endl;
    }
    outFile << "\n";
    return;
}

void Floorplanner::clear()
{
    for (size_t i = 0, end = _blockArray.size(); i < end; ++i) {
        delete _blockArray[i];
    }
    for (size_t i = 0, end = _terminalArray.size(); i < end; ++i) {
        delete _terminalArray[i];
    }
    for (size_t i = 0, end = _netArray.size(); i < end; ++i) {
        delete _netArray[i];
    }
    return;
}

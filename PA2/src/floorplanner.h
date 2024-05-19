#ifndef FLOORPLANNER_H
#define FLOORPLANNER_H

#include <fstream>
#include <vector>
#include <map>
#include <climits>
#include <ctime>
#include "module.h"
using namespace std;

class Floorplanner
{
public:
    // constructor and destructor
    Floorplanner(fstream& blkFile, fstream& netFile, double alpha) :
        _blockNum(0),_netNum(0),_terminalNum(0),_xmax(0),_ymax(0),_bestcost(INT_MAX){
            parseBlk(blkFile);
            parseNet(netFile);
            _alpha = alpha;
        }
    ~Floorplanner() {
        clear();
    }

    // modify method
    void parseBlk(fstream& blkFile);
    void parseNet(fstream& netFile);
    void recur(Node* root);
    void treetograph(Node* root);
    void copyrecur(Node* &node, Node* &bestNode);
    void copytree(bool flag);
    void swapnode(Node* &node1, Node* &node2);
    void deleteandinsert(Node* &node1, Node* &node2);
    void swapsubtree(Node* &node);
    void floorplan();

    // member functions about reporting
    void printSummary() const;
    void reportBlk() const;
    void reportTerm() const;
    void reportNet() const;
    void writeResult(fstream& outFile);

    double getcostFunc(){
        return _alpha*(_chipArea/_AreaNorm)+(1-_alpha)*(_HPWL/_WireNorm);
    };

private:
    double _alpha; //alpha
    double _chipArea;
    int _blockNum;
    int _netNum;
    int _terminalNum;
    double _xmax; //outline
    double _ymax;
    double _HPWL;
    double _chipWidth;
    double _chipHeight;
    double _cost;
    double _bestcost;
    double _bestChipWidth;
    double _bestChipHeight;
    double _bestChipArea;
    double _bestdiff;
    double _runningTime;
    double _AreaNorm;
    double _WireNorm;
    double _bestCostFunc;
    vector<bool>        _bestrotatelist; //rotate condition for best b tree
    vector<Net*>        _netArray;      // net array of the circuit
    vector<Block*>      _blockArray;     // block array of the circuit
    vector<Terminal*>   _terminalArray;
    map<string, int>    _blockName2Id;    // mapping from block name to id
    map<string, int>    _terminalName2Id;   // mapping from terminal name to id
    Node*               _root; //for storing the current b tree
    Node*               _bestRoot; //for storing the best b tree
    vector<vector<double>> doublyLinkedlist; //for finding coordinate y
    // Clean up
    void clear();
};

#endif  
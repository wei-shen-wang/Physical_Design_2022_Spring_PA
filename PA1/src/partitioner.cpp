#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <map>
#include "cell.h"
#include "net.h"
#include "partitioner.h"
using namespace std;


void Partitioner::parseInput(fstream& inFile)
{
    string str;
    // Set balance factor
    inFile >> str;
    _bFactor = stod(str);

    // Set up whole circuit
    while (inFile >> str) {
        if (str == "NET") {
            string netName, cellName, tmpCellName = "";
            inFile >> netName;
            int netId = _netNum;
            _netArray.push_back(new Net(netName));
            _netName2Id[netName] = netId;
            while (inFile >> cellName) {
                if (cellName == ";") {
                    tmpCellName = "";
                    break;
                }
                else {
                    // a newly seen cell
                    if (_cellName2Id.count(cellName) == 0) {
                        int cellId = _cellNum;
                        _cellArray.push_back(new Cell(cellName, 0, cellId));
                        _cellName2Id[cellName] = cellId;
                        _cellArray[cellId]->addNet(netId);
                        _cellArray[cellId]->incPinNum();
                        _netArray[netId]->addCell(cellId);
                        ++_cellNum;
                        tmpCellName = cellName;
                    }
                    // an existed cell
                    else {
                        if (cellName != tmpCellName) {
                            assert(_cellName2Id.count(cellName) == 1);
                            int cellId = _cellName2Id[cellName];
                            _cellArray[cellId]->addNet(netId);
                            _cellArray[cellId]->incPinNum();
                            _netArray[netId]->addCell(cellId);
                            tmpCellName = cellName;
                        }
                    }
                }
            }
            ++_netNum;
        }
    }
    return;
}

void Partitioner::Update_Gain(Node* &head, Node* &tail, Cell* &tempcell, Net* &tempnet, Cell* &basecell){
    int tempgain;
    _moveNum++;
    Node* cur;
    //cell to be moved
    basecell = _cellArray[_maxGainCell->getId()];
    //if unbalanced then choose cell from other part with max gain
    if(_partSize[0]<=(_cellArray.size())*(1-_bFactor)/2){
        _flag = true;
        for(int i = _maxPinNum;i>=-_maxPinNum;i--){
            if(_bList[1][i]->getNext()->getId() != -2){
                basecell = _cellArray[_bList[1][i]->getNext()->getId()];
                _flag = false;
                break;
            }
        }
        
    }
    if(_partSize[0]>=(_cellArray.size())*(1+_bFactor)/2){
        _flag = true;
        for(int i = _maxPinNum;i>=-_maxPinNum;i--){
            if(_bList[0][i]->getNext()->getId() != -2){
                basecell = _cellArray[_bList[0][i]->getNext()->getId()];
                _flag = false;
                break;
            }
        }
        
    }
    //cout << basecell->getGain()<<endl;
    //add the gain
    _accGain += basecell->getGain();
    //add the move to move stack
    _moveStack.push_back(_cellName2Id[basecell->getName()]);
    //lock the cell to be moved
    basecell->lock();
    //remove the cell to be moved from the bucketlist because it will not be used anymore
    removeNode(basecell);
    vector<int> tempnetList = basecell->getNetList();
    vector<int> tempcellList;
    if(basecell->getPart()){
        //f=b,t=a
        _unlockNum[1]--;
        for(int i = 0;i<tempnetList.size();i++){
            tempnet = _netArray[tempnetList[i]];
            tempcellList = tempnet->getCellList();
            if(tempnet->getPartCount(0)==0){
                //ta=0
                for(int j = 0;j<tempcellList.size();j++){
                    tempcell = _cellArray[tempcellList[j]];
                    //increment all free cells on the net
                    if(!tempcell->getLock()){
                        removeNode(tempcell);
                        tempcell->incGain();
                        insertNode(tempcell);
                    }
                }
            }
            else if(tempnet->getPartCount(0)==1){
                //ta=1
                for(int j = 0;j<tempcellList.size();j++){
                    tempcell = _cellArray[tempcellList[j]];
                    //dec the only ta cell on the net
                    if(!tempcell->getLock() && !tempcell->getPart()){
                        removeNode(tempcell);
                        tempcell->decGain();
                        insertNode(tempcell);
                    }
                }
            }
            //fb--,ta++
            tempnet->decPartCount(1);
            tempnet->incPartCount(0);
            if(tempnet->getPartCount(1)==0){
                //fb=0
                for(int j = 0;j<tempcellList.size();j++){
                    tempcell = _cellArray[tempcellList[j]];
                    //dec all free cell on the net
                    if(!tempcell->getLock()){
                        removeNode(tempcell);
                        tempcell->decGain();
                        insertNode(tempcell);
                    }
                }
            }
            else if(tempnet->getPartCount(1)==1){
                //fb=1
                for(int j = 0;j<tempcellList.size();j++){
                    tempcell = _cellArray[tempcellList[j]];
                    //inc the only fb cell on the net
                    if(!tempcell->getLock() && tempcell->getPart()){
                        removeNode(tempcell);
                        tempcell->incGain();
                        insertNode(tempcell);
                    }
                }
            }
        }
        //move the cell
        basecell->move();
        _partSize[1]--;
        _partSize[0]++;
    }
    else{
        //f=a,t=b
        _unlockNum[0]--;
        for(int i = 0;i<tempnetList.size();i++){
            tempnet = _netArray[tempnetList[i]];
            tempcellList = tempnet->getCellList();
            if(tempnet->getPartCount(1)==0){
                //tb=0
                for(int j = 0;j<tempcellList.size();j++){
                    tempcell = _cellArray[tempcellList[j]];
                    //inc all free cells on the net
                    if(!tempcell->getLock()){
                        removeNode(tempcell);
                        tempcell->incGain();
                        insertNode(tempcell);
                    }
                }
            }
            else if(tempnet->getPartCount(1)==1){
                //tb=1
                for(int j = 0;j<tempcellList.size();j++){
                    tempcell = _cellArray[tempcellList[j]];
                    //dec the only tb cell on the net
                    if(!tempcell->getLock() && tempcell->getPart()){
                        removeNode(tempcell);
                        tempcell->decGain();
                        insertNode(tempcell);
                    }
                }
            }
            //fa--,tb++
            tempnet->decPartCount(0);
            tempnet->incPartCount(1);
            if(tempnet->getPartCount(0)==0){
                //fa=0
                for(int j = 0;j<tempcellList.size();j++){
                    tempcell = _cellArray[tempcellList[j]];
                    //dec all free cells on the net
                    if(!tempcell->getLock()){
                        removeNode(tempcell);
                        tempcell->decGain();
                        insertNode(tempcell);
                    }
                }
            }
            else if(tempnet->getPartCount(0)==1){
                //fa=1
                for(int j = 0;j<tempcellList.size();j++){
                    tempcell = _cellArray[tempcellList[j]];
                    //inc the only fa free cell on the net
                    if(!tempcell->getLock() && !tempcell->getPart()){
                        removeNode(tempcell);
                        tempcell->incGain();
                        insertNode(tempcell);
                    }
                }
            }
        }
        basecell->move();
        _partSize[0]--;
        _partSize[1]++;
    }
    if(_accGain>_maxAccGain){
        _maxAccGain = _accGain;
        _bestMoveNum = _moveNum;
        //cout<<"bestmove: "<<_bestMoveNum<<endl;
        //cout<<"maxgain: "<<_maxAccGain<<endl;
    }
}
//insert node to its place of bucketlist
void Partitioner::insertNode(Cell* &tempcell)
{
    Node* head;
    Node* cur;
    Node* tail;
    int tempgain = tempcell->getGain();
    cur = tempcell->getNode();
    if(tempcell->getPart()){
        head = _bList[1][tempgain];
    }
    else{
        head = _bList[0][tempgain];
    }
    tail = head->getNext();
    head->setNext(cur);
    cur->setPrev(head);
    cur->setNext(tail);
    tail->setPrev(cur);
}

//remove the node from the bucket list
void Partitioner::removeNode(Cell* &tempcell)
{
    Node* cur = tempcell->getNode();
    Node* head = cur->getPrev();
    Node* tail = cur->getNext();
    head->setNext(tail);
    tail->setPrev(head);
    cur->setNext(NULL);
    cur->setPrev(NULL);
}

//find the cell with max gain through the bucket list
void Partitioner::updateMaxGainCell()
{
    Node* node0 = NULL;
    Node* node1 = NULL;
    Cell* cell0;
    Cell* cell1;
    for(int i = _maxPinNum;i>=-_maxPinNum;i--){
        if(_bList[0][i]->getNext()->getId() != -2){
            node0 = _bList[0][i]->getNext();
            break;
        }
    }
    for(int i = _maxPinNum;i>=-_maxPinNum;i--){
        if(_bList[1][i]->getNext()->getId() != -2){
            node1 = _bList[1][i]->getNext();
            break;
        }
    }
    //if node0 is not assigned to any node
    if(!node0){
        _maxGainCell = node1;
        return;
    }
    if(!node1){
        _maxGainCell = node0;
        return;
    }
    //cout<<node0->getId()<<" "<<node1->getId()<<endl;
    cell0 = _cellArray[node0->getId()];
    cell1 = _cellArray[node1->getId()];

    if(cell0->getGain()>cell1->getGain()){
        _maxGainCell = node0;
    }
    else if(cell1->getGain()>cell0->getGain()){
        _maxGainCell = node1;
    }
    else{
        //if the gain is equal then balance
        if(_partSize[0]>_partSize[1]){
            _maxGainCell = node0;
        }
        else{
            _maxGainCell = node1;
        }
    }
    //cout<<_maxGainCell->getId()<<endl;
}

//reconstruct the bucketlist
void Partitioner::resetBucketList()
{
    int tempgain;
    Node* tail;
    Node* head;
    Node* cur;

    //delete the head and tail of each bucket
    for(int i = -_maxPinNum;i<=_maxPinNum;i++){
        if(_bList[0][i]){
            delete _bList[0][i]->getPrev();
            delete _bList[1][i]->getPrev();
            delete _bList[0][i];
            delete _bList[1][i];
        }

    }
    _bList[0].clear();
    _bList[1].clear();

    //initialize head and tail for each bucket
    for(int i = -_maxPinNum;i<=_maxPinNum;i++){
        _bList[0][i] = new Node(-1);
        _bList[1][i] = new Node(-1);
        _bList[0][i]->setNext(new Node(-2));
        _bList[0][i]->setPrev(_bList[0][i]->getNext());
        _bList[0][i]->getNext()->setPrev(_bList[0][i]);
        _bList[1][i]->setNext(new Node(-2));
        _bList[1][i]->setPrev(_bList[1][i]->getNext());
        _bList[1][i]->getNext()->setPrev(_bList[1][i]);
    }

    //insert every node into the bucket list
    for(int i = 0;i<_cellArray.size();i++){
        insertNode(_cellArray[i]);
    }
}

void Partitioner::partition()
{
    Cell* tempcell;
    Cell* basecell;
    Net* tempnet;
    int tempgain;
    Node* tail;
    Node* head;
    Node* cur;
    int upper = (_cellArray.size())*(1+_bFactor)/2;
    int lower = (_cellArray.size())*(1-_bFactor)/2;
    int limi = (_cellArray.size())*1/2;
    //give the graph an initial cut
    _moveStack.clear();
    _cutSize = 0;
    _partSize[0] = 0;
    _partSize[1] = 0;
    _bestMoveNum = 0;
    _moveNum = 0;
    _unlockNum[0] = 0;
    _unlockNum[1] = 0;
    _iterNum = 0;
    // for(int i = 0;i<_cellArray.size();i++){
    //     //find the max pin number for construction of the bucket list
    //     _cellArray[i]->setGain(0);
    //     if(_cellArray[i]->getNetList().size()>_maxPinNum){
    //         _maxPinNum = _cellArray[i]->getNetList().size();
    //     }
    // }
    // for(int i = 0;i<partA.size();i++){
    //     _cellArray[partA[i]]->setPart(false);
    // }
    // for(int i = 0;i<partB.size();i++){
    //     _cellArray[partB[i]]->setPart(true);
    // }
    
    for(int i = 0;i<_cellArray.size();i++){
        //find the max pin number for construction of the bucket list
        _cellArray[i]->setGain(0);
        if(_cellArray[i]->getNetList().size()>_maxPinNum){
            _maxPinNum = _cellArray[i]->getNetList().size();
        }
        //bucket
        if(i<=limi){
            _cellArray[i]->setPart(false);
            _partSize[0]++;
            _unlockNum[0]++;
        }
        else{
            _cellArray[i]->setPart(true);
            _partSize[1]++;
            _unlockNum[1]++;
        }
    }

    //calculate the part count for each net and total cutsize and each cell gain
    for(int i = 0;i<_netArray.size();i++){
        vector<int> cellList = _netArray[i]->getCellList();
        _netArray[i]->setPartCount(1,0);
        _netArray[i]->setPartCount(0,0);
        for(int j = 0;j<cellList.size();j++){
            tempcell = _cellArray[cellList[j]];
            if(tempcell->getPart()){
                _netArray[i]->incPartCount(1);
            }
            else{
                _netArray[i]->incPartCount(0);
            }
        }
        if(_netArray[i]->getPartCount(1)>0 && _netArray[i]->getPartCount(0)>0){
            _cutSize++;
        }
        for(int j = 0; j <cellList.size();j++){
            tempcell = _cellArray[cellList[j]];
            if(tempcell->getPart()){
                if(_netArray[i]->getPartCount(1)==1){
                    tempcell->incGain();
                }
                if(_netArray[i]->getPartCount(0)==0){
                    tempcell->decGain();
                }
            }
            else{
                if(_netArray[i]->getPartCount(0)==1){
                    tempcell->incGain();
                }
                if(_netArray[i]->getPartCount(1)==0){
                    tempcell->decGain();
                }
            }
        }
    }
    //cout <<"initial cut size:"<< _cutSize << endl;

    resetBucketList();
    updateMaxGainCell();

    while(true){
        _iterNum++;
        //cout<<"iteration: "<<_iterNum<<" "<<_unlockNum[0]<<" "<<_unlockNum[1]<<endl;
        while(_unlockNum[0]>0 || _unlockNum[1]>0){
            //cout<<_unlockNum[0]<<" "<<_unlockNum[1]<<endl;
            Update_Gain(head,tail,tempcell,tempnet,basecell);
            updateMaxGainCell();
        }
        //cout<<"iteration: "<<_iterNum<<endl;
        //cout <<"gain: "<< _maxAccGain<<endl;
        //cout <<"bestmove: "<<_bestMoveNum<<endl;
        //cout<<"cutsize: "<<_cutSize<<endl;
        //cout<<"A: "<<_partSize[0]<<"  B: "<<_partSize[1]<<endl<<endl;
        for(int i = _moveStack.size()-1;i>=_bestMoveNum;i--){
            _cellArray[_moveStack[i]]->move();
        }

        
        _moveStack.clear();
        _cutSize = 0;
        _partSize[0] = 0;
        _partSize[1] = 0;
        _bestMoveNum = 0;
        _moveNum = 0;
        _unlockNum[0] = 0;
        _unlockNum[1] = 0;
        for(int i = 0;i<_cellArray.size();i++){
            _cellArray[i]->setGain(0);
            _cellArray[i]->unlock();
            if(_cellArray[i]->getPart()){
                _partSize[1]++;
                _unlockNum[1]++;
            }
            else{
                _partSize[0]++;
                _unlockNum[0]++;
            }
        }
        //calculate the part count for each net and total cutsize
        for(int i = 0;i<_netArray.size();i++){
            _netArray[i]->setPartCount(1,0);
            _netArray[i]->setPartCount(0,0);
            vector<int> cellList = _netArray[i]->getCellList();
            for(int j = 0;j<cellList.size();j++){                
                tempcell = _cellArray[cellList[j]];
                if(tempcell->getPart()){
                    _netArray[i]->incPartCount(1);
                }
                else{
                    _netArray[i]->incPartCount(0);
                }
            }
            //cout << _netArray[i]->getName()<< " " << _netArray[i]->getPartCount(0) << " " << _netArray[i]->getPartCount(1)<<endl;
            if(_netArray[i]->getPartCount(1)>0 && _netArray[i]->getPartCount(0)>0){
                _cutSize++;
            }
            for(int j = 0; j <cellList.size();j++){
                tempcell = _cellArray[cellList[j]];
                if(tempcell->getPart()){
                    if(_netArray[i]->getPartCount(1)==1){
                        tempcell->incGain();
                    }
                    if(_netArray[i]->getPartCount(0)==0){
                        tempcell->decGain();
                    }
                }
                else{
                    if(_netArray[i]->getPartCount(0)==1){
                        tempcell->incGain();
                    }
                    if(_netArray[i]->getPartCount(1)==0){
                        tempcell->decGain();
                    }
                }
            }
        }
        resetBucketList();
        updateMaxGainCell();
        
        //cout<<" cutsize: "<<_cutSize<<endl;
        if(_maxAccGain<=0){
            _maxAccGain = 0;
            _accGain = 0;
            break;
        }
        _maxAccGain = 0;
        _accGain = 0;
    }
}

void Partitioner::printSummary() const
{
    cout << endl;
    cout << "==================== Summary ====================" << endl;
    cout << " Cutsize: " << _cutSize << endl;
    cout << " Total cell number: " << _cellNum << endl;
    cout << " Total net number:  " << _netNum << endl;
    cout << " Cell Number of partition A: " << _partSize[0] << endl;
    cout << " Cell Number of partition B: " << _partSize[1] << endl;
    cout << "=================================================" << endl;
    cout << endl;
    return;
}

void Partitioner::reportNet() const
{
    cout << "Number of nets: " << _netNum << endl;
    for (size_t i = 0, end_i = _netArray.size(); i < end_i; ++i) {
        cout << setw(8) << _netArray[i]->getName() << ": ";
        vector<int> cellList = _netArray[i]->getCellList();
        for (size_t j = 0, end_j = cellList.size(); j < end_j; ++j) {
            cout << setw(8) << _cellArray[cellList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::reportCell() const
{
    cout << "Number of cells: " << _cellNum << endl;
    for (size_t i = 0, end_i = _cellArray.size(); i < end_i; ++i) {
        cout << setw(8) << _cellArray[i]->getName() << ": ";
        vector<int> netList = _cellArray[i]->getNetList();
        for (size_t j = 0, end_j = netList.size(); j < end_j; ++j) {
            cout << setw(8) << _netArray[netList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::writeResult(fstream& outFile)
{
    stringstream buff;
    buff << _cutSize;
    outFile << "Cutsize = " << buff.str() << '\n';
    buff.str("");
    buff << _partSize[0];
    outFile << "G1 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 0) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    buff.str("");
    buff << _partSize[1];
    outFile << "G2 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 1) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    return;
}

void Partitioner::clear()
{
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        delete _cellArray[i];
    }
    for (size_t i = 0, end = _netArray.size(); i < end; ++i) {
        delete _netArray[i];
    }
    return;
}

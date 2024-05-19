#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <climits>
#include "partitioner.h"
using namespace std;

int main(int argc, char** argv)
{
    fstream input, output;

    if (argc == 3) {
        input.open(argv[1], ios::in);
        output.open(argv[2], ios::out);
        if (!input) {
            cerr << "Cannot open the input file \"" << argv[1]
                 << "\". The program will be terminated..." << endl;
            exit(1);
        }
        if (!output) {
            cerr << "Cannot open the output file \"" << argv[2]
                 << "\". The program will be terminated..." << endl;
            exit(1);
        }
    }
    else {
        cerr << "Usage: ./fm <input file> <output file>" << endl;
        exit(1);
    }
    std::random_device ma;
    std::mt19937 mt(ma());
    std::uniform_real_distribution<double> rg(0,1);
    double rn = rg(mt);
    int bestCutSize = INT_MAX;
    int bestIndex;
    Partitioner* partitioner = new Partitioner(input);
    // int Cellnum = partitioner->getCellNum();
    // int limi = Cellnum/2;
    // vector<int> partA;
    // vector<int> partB;
    // vector<vector<vector<int>>> all;
    // for(int j = 0;j<10;j++){
    //     partA.clear();
    //     partB.clear();
    //     for(int i = 0;i<Cellnum;i++){
    //         rn = rg(mt);
    //         if(partA.size()==limi){
    //             partB.push_back(i);
    //         }
    //         else if(partB.size()==Cellnum-limi){
    //             partA.push_back(i);
    //         }
    //         else{
    //             if(rn>0.5){
    //                 partB.push_back(i);
    //             }
    //             else{
    //                 partA.push_back(i);
    //             }
    //         }
    //     }
    //     all.push_back({partA,partB});
    //     partitioner->partition(partA,partB);
    //     if(partitioner->getCutSize()<bestCutSize){
    //         bestCutSize = partitioner->getCutSize();
    //         bestIndex = j;
    //     }
    // }
    partitioner->partition();
    
    //partitioner->reportNet();
    //partitioner->reportCell();
    partitioner->printSummary();
    partitioner->writeResult(output);
    delete partitioner;
    return 0;
}

#include "GlobalPlacer.h"
#include "ExampleFunction.h"
#include "NumericalOptimizer.h"
#include <random>
#include <cmath>
#include <ctime>
using namespace std;
GlobalPlacer::GlobalPlacer(Placement &placement)
	:_placement(placement)
{

}


void GlobalPlacer::place()
{
	///////////////////////////////////////////////////////////////////
	// The following example is only for analytical methods.
	// if you use other methods, you can skip and delete it directly.
	//////////////////////////////////////////////////////////////////

	

    unsigned modulenum = _placement.numModules();
    vector<double> x(modulenum*2); // solution vector, size: num_blocks*2 
                         // each 2 variables represent the X and Y dimensions of a block

    //initialize the vector
    double placew = _placement.boundryRight()-_placement.boundryLeft();
    double placeh = _placement.boundryTop()-_placement.boundryBottom();
    double midw = (_placement.boundryRight()+_placement.boundryLeft())/2;
    double midh = (_placement.boundryTop()+_placement.boundryBottom())/2;
    int num = sqrt(modulenum);
    double space = sqrt(modulenum/(placew*placeh));
    int xnum = space*placew;
    int ynum = space*placeh;

    //random
    vector<int> moduleorder(modulenum);
    vector<int> best;
    vector<int> current;
    vector<int> next;
    int wirelength;
    int bestwirelength;
    int nextwirelength;
    int iteration = 0;
    for(int i = 0;i<moduleorder.size();i++){
        moduleorder[i] = i;
    }
    for(int i = 0;i<modulenum;i++){
        int rnum = moduleorder[i];
        int k = i%xnum;
        int j = i/ynum;
        x[2*rnum] = _placement.boundryLeft()+k*placew/xnum;
        x[2*rnum+1] = _placement.boundryBottom()+j*_placement.module(rnum).height();
    }
    for(int i = 0;i<modulenum*2;i+=2){
        _placement.module(i/2).setPosition(x[i],x[i+1]);
    }
    bestwirelength = _placement.computeHpwl();
    wirelength = _placement.computeHpwl();
    std::random_device ma;
    std::mt19937 mt(ma());
    std::uniform_int_distribution<int> rg(0,modulenum-1);
    time_t start = time(NULL); 
    int choice, choice2, temp, rnum, k, j;
    double a, b;
    while((time(NULL)-start)<6600){
        iteration++;
        choice = rg(mt);
        choice2 = rg(mt);
        temp = moduleorder[choice];
        moduleorder[choice] = moduleorder[choice2];
        moduleorder[choice2] = temp;
        for(int i = 0;i<modulenum;i++){
            rnum = moduleorder[i];
            k = i%xnum;
            j = i/ynum;
            a = _placement.boundryLeft()+k*placew/xnum;
            b = _placement.boundryBottom()+j*_placement.module(rnum).height();
            _placement.module(rnum).setPosition(a,b);
        }
        wirelength = _placement.computeHpwl();

        if(wirelength<bestwirelength){
            bestwirelength = wirelength;
            //cout <<"T: "<<T<<endl;
            //cout <<"time: "<<time(NULL)-start<<endl;
            //cout <<"iteration: "<< iteration << " best: " << bestwirelength <<endl;
        }
        else{
            temp = moduleorder[choice];
            moduleorder[choice] = moduleorder[choice2];
            moduleorder[choice2] = temp;
        }
    }
    for(int i = 0;i<modulenum;i++){
        rnum = moduleorder[i];
        k = i%xnum;
        j = i/ynum;
        x[2*rnum] = _placement.boundryLeft()+k*placew/xnum;
        x[2*rnum+1] = _placement.boundryBottom()+j*_placement.module(rnum).height();
    }
    for(int i = 0;i<modulenum*2;i+=2){
        _placement.module(i/2).setPosition(x[i],x[i+1]);
    }
    //random

    // //main version

    // //start from middle
    // for(int i = 0;i<modulenum;i++){
    //     x[2*i] = midw;
    //     x[2*i+1] = midh;
    // }

    // //main version
    // ExampleFunction ef(x,_placement); // require to define the object function and gradient function
    // NumericalOptimizer no(ef);
    // no.setX(x); // set initial solution
    // no.setNumIteration(25); // user-specified parameter
    // no.setStepSizeBound(1000000); // user-specified parameter
    // //20 1000000 quadratic
    // //25 1000000 wa
    // no.solve(); // Conjugate Gradient solver

    // cout << "Current solution:" << endl;
    // for (unsigned i = 0; i < no.dimension(); i++) {
    //     cout << "x[" << i << "] = " << no.x(i) << endl;
    // }
    // int w = placew/num;
    // int h = placeh/num;
    // for(int i = 0;i<modulenum;i++){
    //     int k = x[2*i]/w;
    //     int j = x[2*i+1]/h;
    //     x[2*i] = _placement.boundryLeft()+k*placew/num;
    //     x[2*i+1] = _placement.boundryBottom()+j*placeh/num;
    // }
    // for(int i = 0;i<modulenum*2;i+=2){
    //     _placement.module(i/2).setPosition(no.x(i),no.x(i+1));
    // }
    // cout << "Objective: " << no.objective() << endl;
	////////////////////////////////////////////////////////////////


}


void GlobalPlacer::plotPlacementResult( const string outfilename, bool isPrompt )
{
    ofstream outfile( outfilename.c_str() , ios::out );
    outfile << " " << endl;
    outfile << "set title \"wirelength = " << _placement.computeHpwl() << "\"" << endl;
    outfile << "set size ratio 1" << endl;
    outfile << "set nokey" << endl << endl;
    outfile << "plot[:][:] '-' w l lt 3 lw 2, '-' w l lt 1" << endl << endl;
    outfile << "# bounding box" << endl;
    plotBoxPLT( outfile, _placement.boundryLeft(), _placement.boundryBottom(), _placement.boundryRight(), _placement.boundryTop() );
    outfile << "EOF" << endl;
    outfile << "# modules" << endl << "0.00, 0.00" << endl << endl;
    for( size_t i = 0; i < _placement.numModules(); ++i ){
        Module &module = _placement.module(i);
        plotBoxPLT( outfile, module.x(), module.y(), module.x() + module.width(), module.y() + module.height() );
    }
    outfile << "EOF" << endl;
    outfile << "pause -1 'Press any key to close.'" << endl;
    outfile.close();

    if( isPrompt ){
        char cmd[ 200 ];
        sprintf( cmd, "gnuplot %s", outfilename.c_str() );
        if( !system( cmd ) ) { cout << "Fail to execute: \"" << cmd << "\"." << endl; }
    }
}

void GlobalPlacer::plotBoxPLT( ofstream& stream, double x1, double y1, double x2, double y2 )
{
    stream << x1 << ", " << y1 << endl << x2 << ", " << y1 << endl
           << x2 << ", " << y2 << endl << x1 << ", " << y2 << endl
           << x1 << ", " << y1 << endl << endl;
}

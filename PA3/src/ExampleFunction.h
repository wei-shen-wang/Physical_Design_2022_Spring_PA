#define _GLIBCXX_USE_CXX11_ABI 0
#ifndef EXAMPLEFUNCTION_H
#define EXAMPLEFUNCTION_H

#include "NumericalOptimizerInterface.h"
#include "Placement.h"
class ExampleFunction : public NumericalOptimizerInterface
{
public:
    ExampleFunction(const vector<double> &x,Placement &_placement);

    void evaluateFG(const vector<double> &x, double &f, vector<double> &g);
    void evaluateF(const vector<double> &x, double &f);
    unsigned d;
    Placement p;
    unsigned dimension();

};
#endif // EXAMPLEFUNCTION_H

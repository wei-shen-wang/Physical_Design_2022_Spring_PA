#include "ExampleFunction.h"
#include <iostream>
#include <cmath>
// minimize 3*x^2 + 2*x*y + 2*y^2 + 7
using namespace std;
double lambda = 10;
double alpha = 300;
double beta = 1;
double mb = 0;
int iteration = 0;
//10000 300 quadratic
//100000-10000 1000-300 quadratic
ExampleFunction::ExampleFunction(const vector<double> &x,Placement &_placement)
{
    p = _placement;
    d = x.size();
}

void ExampleFunction::evaluateFG(const vector<double> &x, double &f, vector<double> &g)
{
    f = 0;
    for(int i = 0;i<g.size()/2;i++){
        g[2*i] = 0;
        g[2*i+1] = 0;
    }
    double ka = 0;
    double bp = 0;
    for(int i = 0;i<x.size()/2;i++){
        if(x[2*i]>p.boundryRight()){
            bp+=(x[2*i]-p.boundryRight())*(x[2*i]-p.boundryRight());
        }
        if(x[2*i]<p.boundryLeft()){
            bp+=(x[2*i]-p.boundryLeft())*(x[2*i]-p.boundryLeft());
        }
        if(x[2*i+1]>p.boundryTop()){
            bp+=(x[2*i+1]-p.boundryTop())*(x[2*i+1]-p.boundryTop());
        }
        if(x[2*i+1]<p.boundryBottom()){
            bp+=(x[2*i+1]-p.boundryBottom())*(x[2*i+1]-p.boundryBottom());
        }
    }
    cout << "penalty:    " << bp << endl;
    //density
    // for(int i = 0;i<x.size()/2;i++){
    //     for(int j = 0;j<x.size()/2;j++){
    //         if(i == j){
    //             continue;
    //         }
    //         if(x[2*i]==x[2*j] && x[2*i+1] == x[2*j+1] && i>j){
    //             continue;
    //         }
    //         else if((x[2*j]+p.module(j).width())==(x[2*i]+p.module(i).width()) && x[2*i+1]==x[2*j+1] && i>j){
    //             continue;
    //         }
    //         if(x[2*i]<=x[2*j] && x[2*i+1]<=x[2*j+1]){
    //             if((x[2*i]+p.module(i).width()-x[2*j])>0 && (x[2*i+1]+p.module(i).height()-x[2*j+1])>0){
    //                 f+=(x[2*i]+p.module(i).width()-x[2*j])*(x[2*i+1]+p.module(i).height()-x[2*j+1]);
    //             }
    //         }
    //         else if((x[2*j]+p.module(j).width())<=(x[2*i]+p.module(i).width()) && x[2*i+1]<=x[2*j+1]){
    //             if(x[2*i]<(x[2*j]+p.module(j).width()) && x[2*j+1]<(x[2*i+1]+p.module(i).height())){
    //                 f+=(x[2*j]+p.module(j).width()-x[2*i])*(x[2*i+1]+p.module(i).height()-x[2*j+1]);
    //             }
    //         }
            
    //     }
    // }

    //bin density
    double wb,wv,ax,bx,ay,by,xi,xb,yi,yb,hb,hv,dx,dy,ox,oy,bin,gox,goy;
    wb=wv=ax=bx=ay=by=xi=xb=yi=yb=hb=hv=dx=dy=ox=oy=bin=gox=goy=0;
    double placew = p.boundryRight()-p.boundryLeft();
    double placeh = p.boundryTop()-p.boundryBottom();
    int modulenum = p.numModules();
    int num = sqrt(modulenum);
    double space = sqrt(modulenum/(placew*placeh));
    int xnum = space*placew;
    int ynum = space*placeh;
    wb = placew/xnum;
    hb = placeh/ynum;
    for(int i = 0;i<xnum;i++){
        for(int j = 0;j<ynum;j++){
            xb = p.boundryLeft()+i*placew/xnum;
            yb = p.boundryBottom()+j*placeh/ynum;
            bin = 0;
            for(int k = 0;k<modulenum;k++){
                xi = x[2*k];
                yi = x[2*k+1];
                wv = p.module(k).width();
                hv = p.module(k).height();
                dx = xi-xb;
                dy = yi-yb;
                ax = 4/((wv+2*wb)*(wv+4*wb));
                bx = 2/(wb*(wv+4*wb));
                ay = 4/((hv+2*hb)*(hv+4*hb));
                by = 2/(hb*(hv+4*hb));
                ox = 0;
                oy = 0;
                if(dx<0){
                    dx = -dx;
                }
                if(dy<0){
                    dy = -dy;
                }
                if(0<=dx && dx<=(wb+wv/2)){
                    ox = 1-ax*dx*dx;
                }
                else if((wb+wv/2)<dx && dx<(wv/2+2*wb)){
                    ox = bx*(dx-wv/2-2*wb)*(dx-wv/2-2*wb);
                }
                if(0<=dy && dy<=(hb+hv/2)){
                    oy = 1-ay*dy*dy;
                }
                else if((hb+hv/2)<dy && dy<(hv/2+2*hb)){
                    oy = by*(dy-hv/2-2*hb)*(dy-hv/2-2*hb);
                }
                bin += ox*oy;
            }
            f+=(bin-mb)*(bin-mb);

            for(int k = 0;k<modulenum;k++){
                xi = x[2*k];
                yi = x[2*k+1];
                wv = p.module(k).width();
                hv = p.module(k).height();
                dx = xi-xb;
                dy = yi-yb;
                ax = 4/((wv+2*wb)*(wv+4*wb));
                bx = 2/(wb*(wv+4*wb));
                ay = 4/((hv+2*hb)*(hv+4*hb));
                by = 2/(hb*(hv+4*hb));
                ox = 0;
                oy = 0;
                gox = 0;
                goy = 0;
                if(dx<0){
                    dx = -dx;
                }
                if(dy<0){
                    dy = -dy;
                }
                if(0<=dx && dx<=(wb+wv/2)){
                    ox = 1-ax*dx*dx;
                }
                else if((wb+wv/2)<dx && dx<(wv/2+2*wb)){
                    ox = bx*(dx-wv/2-2*wb)*(dx-wv/2-2*wb);
                }
                if(0<=dy && dy<=(hb+hv/2)){
                    oy = 1-ay*dy*dy;
                }
                else if((hb+hv/2)<dy && dy<(hv/2+2*hb)){
                    oy = by*(dy-hv/2-2*hb)*(dy-hv/2-2*hb);
                }
                if(0<=dx && dx<=(wb+wv/2)){
                    gox = -2*ax*(xi-xb);
                }
                else if((wb+wv/2)<dx && dx<(wv/2+2*wb)){
                    if(xi>xb){
                        gox = 2*bx*(xi-xb-wv/2-2*wb);
                    }
                    else{
                        gox = 2*bx*(xi-xb+wv/2+2*wb);
                    }
                }
                if(0<=dy && dy<=(hb+hv/2)){
                    goy = -2*ay*(yi-yb);
                }
                else if((hb+hv/2)<dy && dy<(hv/2+2*hb)){
                    if(yi>yb){
                        goy = 2*by*(yi-yb-hv/2-2*hb);
                    }
                    else{
                        goy = 2*by*(yi-yb+hv/2+2*hb);
                    }
                }
                g[2*k] += 2*lambda*(bin-mb)*oy*gox;
                g[2*k+1] += 2*lambda*(bin-mb)*ox*goy;
            }
        }
    }


    cout.precision(14);
    cout << "overlap:    " << f << endl;
    double xmaxmom = 0;
    double xmaxson = 0;
    double xminmom = 0;
    double xminson = 0;
    double ymaxmom = 0;
    double ymaxson = 0;
    double yminmom = 0;
    double yminson = 0;
    for(int i = 0;i<p.numNets();i++){
        xmaxmom = 0;
        xmaxson = 0;
        xminmom = 0;
        xminson = 0;
        ymaxmom = 0;
        ymaxson = 0;
        yminmom = 0;
        yminson = 0;
        for(int j = 0;j<p.net(i).numPins();j++){
            //for(int k = j+1;k<p.net(i).numPins();k++){
                int m1 = p.net(i).pin(j).moduleId();
                //int m2 = p.net(i).pin(k).moduleId();
                //wa
                if(!isnan(x[2*m1]*exp(x[2*m1]/100))){
                    xmaxson += x[2*m1]*exp(x[2*m1]/100);
                }
                if(!isnan(exp(x[2*m1]/100))){
                    xmaxmom += exp(x[2*m1]/100);
                }
                if(!isnan(x[2*m1]*exp(-x[2*m1]/100))){
                    xminson += x[2*m1]*exp(-x[2*m1]/100);
                }
                if(!isnan(exp(-x[2*m1]/100))){
                    xminmom += exp(-x[2*m1]/100);
                }
                if(!isnan(x[2*m1+1]*exp(x[2*m1+1]/100))){
                    ymaxson += x[2*m1+1]*exp(x[2*m1+1]/100);
                }
                if(!isnan(exp(x[2*m1+1]/100))){
                    ymaxmom += exp(x[2*m1+1]/100);
                }
                if(!isnan(x[2*m1+1]*exp(-x[2*m1+1]/100))){
                    yminson += x[2*m1+1]*exp(-x[2*m1+1]/100);
                }
                if(!isnan(exp(-x[2*m1+1]/100))){
                    yminmom += exp(-x[2*m1+1]/100);
                }
                
                //quadratic
                //ka+=0.5*(x[2*m1]-x[2*m2])*(x[2*m1]-x[2*m2]);
                //ka+=0.5*(x[2*m1+1]-x[2*m2+1])*(x[2*m1+1]-x[2*m2+1]);
            //}
        }
        //wa
        for(int j = 0;j<p.net(i).numPins();j++){
            int m1 = p.net(i).pin(j).moduleId();
            double xmax = ((exp(x[2*m1]/100)+(x[2*m1]/100)*exp(x[2*m1]/100))*xmaxmom-(exp(x[2*m1]/100)/100)*xmaxson)/((xmaxmom)*(xmaxmom));
            double xmin = ((exp(-x[2*m1]/100)-(x[2*m1]/100)*exp(-x[2*m1]/100))*xminmom+(exp(-x[2*m1]/100)/100)*xminson)/((xminmom)*(xminmom));
            double ymax = ((exp(x[2*m1+1]/100)+(x[2*m1+1]/100)*exp(x[2*m1+1]/100))*ymaxmom-(exp(x[2*m1+1]/100)/100)*ymaxson)/((ymaxmom)*(ymaxmom));
            double ymin = ((exp(-x[2*m1+1]/100)-(x[2*m1+1]/100)*exp(-x[2*m1+1]/100))*yminmom+(exp(-x[2*m1+1]/100)/100)*yminson)/((yminmom)*(yminmom));
            if(isnan(xmax) || isinf(xmax)){
                xmax = 0;
            }
            if(isnan(ymax) || isinf(ymax)){
                ymax = 0;
            }
            if(isnan(xmin) || isinf(xmin)){
                xmin = 0;
            }
            if(isnan(ymin) || isinf(ymin)){
                ymin = 0;
            }
            g[2*m1] += beta*xmax;
            g[2*m1] -= beta*xmin;
            g[2*m1+1] += beta*ymax;
            g[2*m1+1] -= beta*ymin;
            //cout << g[2*m1]<<endl;
        }
        ka+=(xmaxson/xmaxmom)-xminson/xminmom+ymaxson/ymaxmom-yminson/yminmom;
    }
    //alpha = ka/(bp*1000);
    //lambda = ka/(f*10);
    f = f*lambda + ka*beta + bp*alpha;
    cout.precision(14);
    cout << "wirelength: " << ka <<endl;

    for(int i = 0;i<x.size()/2;i++){
        if(x[2*i]>p.boundryRight()){
            g[2*i]+=2*alpha*(x[2*i]-p.boundryRight());
        }
        if(x[2*i]<p.boundryLeft()){
            g[2*i]+=2*alpha*(x[2*i]-p.boundryLeft());
        }
        if(x[2*i+1]>p.boundryTop()){
            g[2*i+1]+=2*alpha*(x[2*i+1]-p.boundryTop());
        }
        if(x[2*i+1]<p.boundryBottom()){
            g[2*i+1]+=2*alpha*(x[2*i+1]-p.boundryBottom());
        }
    }
    //density
    // for(int i = 0;i<g.size()/2;i++){
    //     for(int j = 0;j<x.size()/2;j++){
    //         if(i == j){
    //             continue;
    //         }
    //         if(x[2*i]==x[2*j] && x[2*i+1] == x[2*j+1] && i>j){
    //             continue;
    //         }
    //         else if((x[2*j]+p.module(j).width())==(x[2*i]+p.module(i).width()) && x[2*i+1]==x[2*j+1] && i>j){
    //             continue;
    //         }
    //         if(x[2*i]<=x[2*j] && x[2*i+1]<=x[2*j+1]){
    //             if((x[2*i]+p.module(i).width()-x[2*j])>0 && (x[2*i+1]+p.module(i).height()-x[2*j+1])>0){
    //                 g[2*i]+=lambda*(x[2*i+1]+p.module(i).height()-x[2*j+1]);
    //                 g[2*i+1]+=lambda*(x[2*i]+p.module(i).width()-x[2*j]);
    //                 g[2*j]+=-lambda*(x[2*i+1]+p.module(i).height()-x[2*j+1]);
    //                 g[2*j+1]+=-lambda*(x[2*i]+p.module(i).width()-x[2*j]);
    //             }
    //         }
    //         else if((x[2*j]+p.module(j).width())<=(x[2*i]+p.module(i).width()) && x[2*i+1]<=x[2*j+1]){
    //             if(x[2*i]<(x[2*j]+p.module(j).width()) && x[2*j+1]<(x[2*i+1]+p.module(i).height())){
    //                 g[2*i]+=-lambda*(x[2*i+1]+p.module(i).height()-x[2*j+1]);
    //                 g[2*i+1]+=lambda*(x[2*j]+p.module(j).width()-x[2*i]);
    //                 g[2*j]+=lambda*(x[2*i+1]+p.module(i).height()-x[2*j+1]);
    //                 g[2*j+1]+=-lambda*(x[2*j]+p.module(j).width()-x[2*i]);
    //             }
    //         }
    //     }
    // }
    // // quadratic
    // for(int i = 0;i<p.numNets();i++){
    //     for(int j = 0;j<p.net(i).numPins();j++){
    //         for(int k = j+1;k<p.net(i).numPins();k++){
    //             int m1 = p.net(i).pin(j).moduleId();
    //             int m2 = p.net(i).pin(k).moduleId();
    //             double for2i = x[2*m1]-x[2*m2];
    //             double for2i1 = x[2*m1+1]-x[2*m2+1];
    //             g[2*m1]+=(beta*for2i);
    //             g[2*m1+1]+=(beta*for2i1);
    //             g[2*m2]+=(-beta*for2i);
    //             g[2*m2+1]+=(-beta*for2i1);
    //             //cout<<g[2*m1]<<endl;
    //         }
            
    //     }
    // }
    iteration++;
    // if(iteration%5==0){
    //     lambda = lambda*10;
    // }
    cout << "iteration: " << iteration<<endl;
}

void ExampleFunction::evaluateF(const vector<double> &x, double &f)
{
    //quadratic
    f = 0;
    double ka = 0;
    double bp = 0;
    for(int i = 0;i<x.size()/2;i++){
        if(x[2*i]>p.boundryRight()){
            bp+=(x[2*i]-p.boundryRight())*(x[2*i]-p.boundryRight());
        }
        if(x[2*i]<p.boundryLeft()){
            bp+=(x[2*i]-p.boundryLeft())*(x[2*i]-p.boundryLeft());
        }
        if(x[2*i+1]>p.boundryTop()){
            bp+=(x[2*i+1]-p.boundryTop())*(x[2*i+1]-p.boundryTop());
        }
        if(x[2*i+1]<p.boundryBottom()){
            bp+=(x[2*i+1]-p.boundryBottom())*(x[2*i+1]-p.boundryBottom());
        }
    }

    //bin density
    double wb,wv,ax,bx,ay,by,xi,xb,yi,yb,hb,hv,dx,dy,ox,oy,bin;
    wb=wv=ax=bx=ay=by=xi=xb=yi=yb=hb=hv=dx=dy=ox=oy=bin=0;
    double placew = p.boundryRight()-p.boundryLeft();
    double placeh = p.boundryTop()-p.boundryBottom();
    int modulenum = p.numModules();
    int num = sqrt(modulenum);
    double space = sqrt(modulenum/(placew*placeh));
    int xnum = space*placew;
    int ynum = space*placeh;
    wb = placew/xnum;
    hb = placeh/ynum;
    for(int i = 0;i<xnum;i++){
        for(int j = 0;j<ynum;j++){
            xb = p.boundryLeft()+i*placew/xnum;
            yb = p.boundryBottom()+j*placeh/ynum;
            bin = 0;
            for(int k = 0;k<modulenum;k++){
                xi = x[2*k];
                yi = x[2*k+1];
                wv = p.module(k).width();
                hv = p.module(k).height();
                dx = xi-xb;
                dy = yi-yb;
                ax = 4/((wv+2*wb)*(wv+4*wb));
                bx = 2/(wb*(wv+4*wb));
                ay = 4/((hv+2*hb)*(hv+4*hb));
                by = 2/(hb*(hv+4*hb));
                ox = 0;
                oy = 0;
                if(dx<0){
                    dx = -dx;
                }
                if(dy<0){
                    dy = -dy;
                }
                if(0<=dx && dx<=(wb+wv/2)){
                    ox = 1-ax*dx*dx;
                }
                else if((wb+wv/2)<dx && dx<(wv/2+2*wb)){
                    ox = bx*(dx-wv/2-2*wb)*(dx-wv/2-2*wb);
                }
                if(0<=dy && dy<=(hb+hv/2)){
                    oy = 1-ay*dy*dy;
                }
                else if((hb+hv/2)<dy && dy<(hv/2+2*hb)){
                    oy = by*(dy-hv/2-2*hb)*(dy-hv/2-2*hb);
                }
                bin += ox*oy;
            }
            f+=(bin-mb)*(bin-mb);
        }
    }
    //density
    // for(int i = 0;i<x.size()/2;i++){
    //     for(int j = 0;j<x.size()/2;j++){
    //         if(i == j){
    //             continue;
    //         }
    //         if(x[2*i]==x[2*j] && x[2*i+1] == x[2*j+1] && i>j){
    //             continue;
    //         }
    //         else if((x[2*j]+p.module(j).width())==(x[2*i]+p.module(i).width()) && x[2*i+1]==x[2*j+1] && i>j){
    //             continue;
    //         }
    //         if(x[2*i]<=x[2*j] && x[2*i+1]<=x[2*j+1]){
    //             if((x[2*i]+p.module(i).width()-x[2*j])>0 && (x[2*i+1]+p.module(i).height()-x[2*j+1])>0){
    //                 f+=(x[2*i]+p.module(i).width()-x[2*j])*(x[2*i+1]+p.module(i).height()-x[2*j+1]);
    //             }
    //         }
    //         else if((x[2*j]+p.module(j).width())<=(x[2*i]+p.module(i).width()) && x[2*i+1]<=x[2*j+1]){
    //             if(x[2*i]<(x[2*j]+p.module(j).width()) && x[2*j+1]<(x[2*i+1]+p.module(i).height())){
    //                 f+=(x[2*j]+p.module(j).width()-x[2*i])*(x[2*i+1]+p.module(i).height()-x[2*j+1]);
    //             }
    //         }

    //     }
    // }
    double xmaxmom = 0;
    double xmaxson = 0;
    double xminmom = 0;
    double xminson = 0;
    double ymaxmom = 0;
    double ymaxson = 0;
    double yminmom = 0;
    double yminson = 0;
    for(int i = 0;i<p.numNets();i++){
        xmaxmom = 0;
        xmaxson = 0;
        xminmom = 0;
        xminson = 0;
        ymaxmom = 0;
        ymaxson = 0;
        yminmom = 0;
        yminson = 0;
        for(int j = 0;j<p.net(i).numPins();j++){
            //for(int k = j+1;k<p.net(i).numPins();k++){
                int m1 = p.net(i).pin(j).moduleId();
                //int m2 = p.net(i).pin(k).moduleId();
                //wa
                if(!isnan(x[2*m1]*exp(x[2*m1]/100))){
                    xmaxson += x[2*m1]*exp(x[2*m1]/100);
                }
                if(!isnan(exp(x[2*m1]/100))){
                    xmaxmom += exp(x[2*m1]/100);
                }
                if(!isnan(x[2*m1]*exp(-x[2*m1]/100))){
                    xminson += x[2*m1]*exp(-x[2*m1]/100);
                }
                if(!isnan(exp(-x[2*m1]/100))){
                    xminmom += exp(-x[2*m1]/100);
                }
                if(!isnan(x[2*m1+1]*exp(x[2*m1+1]/100))){
                    ymaxson += x[2*m1+1]*exp(x[2*m1+1]/100);
                }
                if(!isnan(exp(x[2*m1+1]/100))){
                    ymaxmom += exp(x[2*m1+1]/100);
                }
                if(!isnan(x[2*m1+1]*exp(-x[2*m1+1]/100))){
                    yminson += x[2*m1+1]*exp(-x[2*m1+1]/100);
                }
                if(!isnan(exp(-x[2*m1+1]/100))){
                    yminmom += exp(-x[2*m1+1]/100);
                }
                // //quadratic
                // ka+=0.5*(x[2*m1]-x[2*m2])*(x[2*m1]-x[2*m2]);
                // ka+=0.5*(x[2*m1+1]-x[2*m2+1])*(x[2*m1+1]-x[2*m2+1]);
            //}
        }
        //wa
        ka+=(xmaxson/xmaxmom)-xminson/xminmom+ymaxson/ymaxmom-yminson/yminmom;
    }
    //alpha = ka/(bp*1000);
    //cout << alpha << " " << lambda << endl;
    //lambda = ka/(f*10);
    cout << ka << endl;
    f = f*lambda + (ka*beta) + bp*alpha;
}

unsigned ExampleFunction::dimension()
{
    return d; // num_blocks*2 
    // each two dimension represent the X and Y dimensions of each block
}

#ifndef MODULE_H
#define MODULE_H

#include <vector>
#include <string>
using namespace std;

class Node
{
    friend class Block;

public:
    // Constructor and destructor
    Node(string& name) :
        _name(name), _parent(NULL), _left(NULL), _right(NULL) { }
    ~Node() { }

    // Basic access methods
    string getName() const  { return _name; }
    Node* getLeft() const   { return _left; }
    Node* getRight() const  { return _right; }
    Node* getParent() const  { return _parent; }

    // Set functions
    void setName(string& name){ _name = name; }
    void setLeft(Node* left)  { _left = left; }
    void setRight(Node* right){ _right = right; }
    void setParent(Node* parent){ _parent = parent; }

private:
    string      _name;    // id of the node (indicating the cell)
    Node*       _left;  // pointer to the previous node
    Node*       _right;  // pointer to the next node
    Node*       _parent;
};

class Terminal
{
public:
    // constructor and destructor
    Terminal(string& name, double x,double y) :
        _name(name), _x1(x), _y1(y), _x2(x), _y2(y) { }
    ~Terminal()  { }

    // basic access methods
    const string getName()  { return _name; }
    const double getX1()    { return _x1; }
    const double getX2()    { return _x2; }
    const double getY1()    { return _y1; }
    const double getY2()    { return _y2; }

    // set functions
    void setName(string& name) { _name = name; }
    void setPos(double x1, double y1, double x2, double y2) {
        _x1 = x1;   _y1 = y1;
        _x2 = x2;   _y2 = y2;
    }

protected:
    string      _name;      // module name
    double      _x1;        // min x coordinate of the terminal
    double      _y1;        // min y coordinate of the terminal
    double      _x2;        // max x coordinate of the terminal
    double      _y2;        // max y coordinate of the terminal
};


class Block : public Terminal
{
public:
    // constructor and destructor
    Block(string& name, double w, double h) :
        Terminal(name, 0, 0), _w(w), _h(h), _rotate(false) { 
            _node = new Node(name);
            _bestNode = new Node(name);
        }
    ~Block() { }

    // basic access methods
    const double getWidth()  { return _rotate? _h: _w; }
    const double getHeight() { return _rotate? _w: _h; }
    const int getArea()  { return _h * _w; }
    static int getMaxX() { return _maxX; }
    static int getMaxY() { return _maxY; }
    Node* getNode() const   { return _node; }
    Node* getBestNode() const { return _bestNode; }
    bool getRotate() const {return _rotate; }
    // set functions
    void setWidth(double w)         { _w = w; }
    void setHeight(double h)        { _h = h; }
    void setNode(Node* node)     { _node = node; }
    void setBestNode(Node* node) { _bestNode = node; }
    void setRotate(bool rotate)  { _rotate = rotate; }
    void rotate()                { _rotate = !_rotate; }
    static void setMaxX(int x)   { _maxX = x; }
    static void setMaxY(int y)   { _maxY = y; }


private:
    double          _w;         // width of the block
    double          _h;         // height of the block
    static int   _maxX;      // maximum x coordinate for all blocks
    static int   _maxY;      // maximum y coordinate for all blocks
    bool         _rotate;
    Node*        _node;
    Node*        _bestNode;
};


class Net
{
public:
    // constructor and destructor
    Net()   { }
    ~Net()  { }

    // basic access methods
    const vector<Terminal*> getTermList()   { return _termList; }

    // modify methods
    void addTerm(Terminal* term) { _termList.push_back(term); }

    // other member functions
    double calcHPWL(){
        double ans = 0;
        double minx = INT_MAX;
        double maxx = 0;
        double miny = INT_MAX;
        double maxy = 0;
        double avgx;
        double avgy;
        for(int i = 0;i<_termList.size();i++){
            avgx = (double)((double)(_termList[i]->getX1())+(double)(_termList[i]->getX2()))/2;
            avgy = (double)((double)(_termList[i]->getY1())+(double)(_termList[i]->getY2()))/2;
            if(avgx<minx){
                minx = avgx;
            }
            if(avgy<miny){
                miny = avgy;
            }
            if(avgx>maxx){
                maxx = avgx;
            }
            if(avgy>maxy){
                maxy = avgy;
            }
        }
        ans = (double)(maxx - minx + maxy - miny);
        return ans;
    };

private:
    vector<Terminal*>   _termList;  // list of terminals the net is connected to
};

#endif  // MODULE_H

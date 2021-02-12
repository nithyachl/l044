#include <iostream>
#include <vector>
#include <utility>
#include <fstream>
#include <iterator>
#include <iomanip>
#include <ios>
#include <cmath>
#include <list>
#include <cstdio>
#include <ctime>
#include <math.h>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <tuple>
#include <chrono>
#include <unordered_map>
#include <queue>

using namespace std;


class Point {

public:
    double xval;
    double yval;
    int cluster;
    double minDist;
    vector<int> possibleCentroids;




    Point()
    {
        cluster = -1;
        minDist = __DBL_MAX__;
    }
//    Point(const Point &other)
//    {
//        this->xval = other.xval;
//        this->yval = other.yval;
//        this->cluster = other.cluster;
//        this->minDist = other.minDist;
//    }

    double distance(Point p) {
        return sqrt((p.xval - xval) * (p.xval - xval) + (p.yval - yval) * (p.yval - yval));
    }

};

void drawCircle(double r, double xcenter, double ycenter, int color);

const int k = 2;
int n = 0;
int in = 0;

struct Node
{
    int point[200];
    Node *left, *right;
    int xmin ;
    int xmax ;
    int ymin ;
    int ymax ;
    Point p;


};


Node *newNode(int arr[], int xmin, int xmax, int ymin, int ymax, int i1, int i2)
{
    struct Node* temp = new Node;

    for (int i=0; i<k; i++)
        temp->point[i] = arr[i];

    temp->left = temp->right = nullptr;
    temp->xmin = xmin;
    temp->xmax = xmax;

    temp->ymin = ymin;
    temp->ymax = ymax;
    temp->p.xval = i1;
    temp->p.yval = i2;

    return temp;
}

void part4();

void traverse(Node *root, int level, int i);

void drawExtendedHorizontal(int i, int i1, int i2);

void drawExtendedVertical(int i, int i1, int i2);

void generatePoints();

void readPoints();
std::vector<Point> total;



void kMeansClustering(vector<Point> *pVector);

int getHeight(Node *pNode);

Node *insertRec(Node *root, int point[], unsigned int depth, int xmin, int xmax, int ymin, int ymax, int i, int i1)
{

    if (root == nullptr)
        return newNode(point, xmin, xmax, ymin, ymax, i, i1);


    unsigned cd = depth % 2;
    cout << "depth: " << depth << endl;

    if (point[cd] < (root->point[cd])) {
        cout << "root point left: " << root->point[cd] << endl;
        if(depth%2 == 0)
            root->left = insertRec(root->left, point, depth + 1, xmin, root->point[cd], ymin, ymax, i, i1);
        else
            root->left = insertRec(root->left, point, depth + 1, xmin, xmax, ymin, root->point[cd], i, i1);

    }
    else {
        cout << "root point right: " << root->point[cd] << endl;
        if(depth%2 == 0)
            root->right = insertRec(root->right, point, depth + 1, root->point[cd], xmax, ymin, ymax, 0, 0);
        else
            root->right = insertRec(root->right, point, depth + 1, xmin, xmax, root->point[cd], ymax, 0, 0);
    }
    return root;
}


Node *insert(Node *root, int point[], int i, int i1)
{
    return insertRec(root, point, 0, 0, 800, 0, 800, i, i1);
}

int matrix[800][800];


int main()
{
    std::ofstream test;

    test.open ("diagram.ppm");
    test << "P3" << endl << 800 << ' ' << 800 << endl << "1" << endl;

    for (int n = 0; n < 800; ++n) {
        for (int m = 0; m < 800; ++m) {
            matrix[n][m] = 1;
        }
    }

    srand( time(NULL) );

    part4();

    for (int n = 0; n < 800; ++n) {
        for (int m = 0; m < 800; ++m) {

            if(matrix[n][m] == 4)
            {
                test<< 1 << " "<< 0 << " " << 0 << "     ";
            }

            else if(matrix[n][m] == 9)
            {
                test<< 0 << " "<< 0 << " " << 1 << "     ";
            }

            else
                test<< matrix[n][m] << " "<< matrix[n][m] << " " << matrix[n][m] << "     ";          //std::to_string(matrix[n][m]) << " ";

        }
        test << endl;
    }

    test.close();

    return 0;

}

void part4() {
    string option;
    int gen;

    cout << "Would you like points to be generated(type yes or no): ";
    cin >> option;


    if((option == "Yes") || (option == "yes") )
    {
        gen = 1;
        cout << "Generating Points." << endl;
    }

    if((option == "No") || (option == "no") )
    {
        gen = 0;
    }

    //  cout << "Point1!" << endl;
    struct Node *root = nullptr;
//    root->xmax =1;
//    root->ymax =1;
//
//    root->xmin =0;
//    root->ymin =0;
//    root->xory = 0;
    if(gen == 1) {
        generatePoints();
    }
    else
    {
        readPoints();
          cout << "Point2!" << endl;
    }

    int points[n][k];
    for(int i = 0; i< n; i++)
    {
        points[i][0] = (int)(total[i].xval *800);
        points[i][1] = (int)(total[i].yval *800);
          cout << points[i][0] << "  " << points[i][1]<< endl;
    }


    int n = sizeof(points)/sizeof(points[1]);

    for (int i=0; i<n; i++)
    {
        root = insert(root, points[i], points[i][0], points[i][1]);
        drawCircle(2,points[i][0], points[i][1],0 );

    }
    // cout << "Point3!" << endl;

    vector<Point> centroids;
    for (int i = 0; i < 5; ++i) {
        Point p;
        p.xval = points[i][0];
        p.yval = points[i][1];
        centroids.push_back(p);
        root->p.possibleCentroids.push_back(i);
    }


    while (in <= n) {
        //  cout << "Point4!" << endl;
        traverse(root, in, 0);
        in++;
        //cout << "Point5!" << endl;
    }



}



void readPoints() {

    ifstream inputFile("points.txt");

    double x,y;


    while (inputFile>>x>>y) {

        Point generate;
        generate.xval = x;
        generate.yval = y;
        total.push_back(generate);
        n = n+1;

    }


    inputFile.close();

}

void generatePoints() {
    std::ofstream myfile;
    myfile.open ("points.txt");

    for(int i = 0; i < 10; i++)
    {
        Point generate;
        generate.xval =  ((double) rand() / (RAND_MAX));
        generate.yval =  ((double) rand() / (RAND_MAX));

        total.push_back(generate);
        myfile << std::fixed << std::setprecision(23) << generate.xval << "  " << generate.yval<< endl;

        n= n+1;
    }

    myfile.close();

}

void traverse(Node *root, int level, int i)
{
    // << "in traverse!" << endl;
    if (root == NULL)
        return;
    if (level == 1) {
        cout << "Traverse: " << root->point[0] << " " << root->point[1] << endl;

        cout << "X: " << root->xmin << " " << root->xmax << endl;
        cout << "Y: " << root->ymin << " " << root->ymax << endl;

        if(in%2 == 1 && i==0) {
            //drawLine(root->point[0] * 40, 0, root->point[0] * 40, root->point[1] * 40);
            // << "in level1!" << endl;
           // drawExtendedVertical(root->point[0] , root->point[1] );
            drawExtendedVertical(root->ymin, root->ymax, root->point[0]);

        }

        if(in%2 == 1 && i==1) {
            //drawLine(root->point[0] * 40, 800, root->point[0] * 40, root->point[1] * 40);
            drawExtendedVertical(root->ymin, root->ymax, root->point[0]);
        }

        if(in%2 == 0 && i == 0) {
            //drawLine(root->point[0] * 40, root->point[1] * 40, 0, root->point[1] * 40);
            drawExtendedHorizontal(root->xmin, root->xmax, root->point[1]); //blue
        }
        if(in%2 == 0 && i == 1) {
            //drawLine(root->point[0] * 40, root->point[1] * 40, 800, root->point[1] * 40);
            drawExtendedHorizontal(root->xmin, root->xmax, root->point[1]); //blue
        }
    }
    else if (level > 1)
    {

        traverse(root->left, level - 1, 0);
        traverse(root->right, level - 1, 1);
    }
}

void drawExtendedVertical(int y1, int y2, int xval) {


    int x1 = xval;
    int x2 = xval;




    int dx = x2 - x1 ;
    int dy = y2 - y1 ;
    int dx1 = 0, dy1 = 0, dx2 = 0, dy2 = 0 ;
    if (dy<0){
        dy1 = -1 ;
    }
    else if (dy>0){
        dy1 = 1 ;
    }
    if (dx<0) {
        dx1 = -1;
    }
    else if (dx>0){
        dx1 = 1 ;}
    if (dy<0){
        dy2 = -1;
    }
    else if (dy>0){
        dy2 = 1 ;
    }
    int longest = abs(dy) ;
    int shortest = abs(dx) ;
    if (!(longest>shortest)) {
        longest = abs(dx) ;
        shortest = abs(dy) ;
        if (dx<0)
        {
            dx2 = -1 ;
        }

        else if (dx>0) {
            dx2 = 1;
        }
        dy2 = 0 ;
    }
    int numerator = longest >> 1 ;

    for (int i=0;i<=longest;i++) {
        matrix[y1][x1] = 4;
        numerator += shortest ;
        if (!(numerator<longest)) {
            numerator -= longest ;
            x1 += dx1 ;
            y1 += dy1 ;
        } else {
            x1 += dx2 ;
            y1 += dy2 ;
        }
    }

}

void drawExtendedHorizontal(int x2, int x1, int yval) {

    int y2 = yval;
    int y1 = yval;

    int dx = x2 - x1 ;
    int dy = y2 - y1 ;
    int dx1 = 0, dy1 = 0, dx2 = 0, dy2 = 0 ;
    if (dy<0){
        dy1 = -1 ;
    }
    else if (dy>0){
        dy1 = 1 ;
    }
    if (dx<0) {
        dx1 = -1;
    }
    else if (dx>0){
        dx1 = 1 ;}
    if (dy<0){
        dy2 = -1;
    }
    else if (dy>0){
        dy2 = 1 ;
    }
    int longest = abs(dy) ;
    int shortest = abs(dx) ;
    if (!(longest>shortest)) {
        longest = abs(dx) ;
        shortest = abs(dy) ;
        if (dx<0)
        {
            dx2 = -1 ;
        }

        else if (dx>0) {
            dx2 = 1;
        }
        dy2 = 0 ;
    }
    int numerator = longest >> 1 ;

    for (int i=0;i<=longest;i++) {
        matrix[y1][x1] = 9;
        numerator += shortest ;
        if (!(numerator<longest)) {
            numerator -= longest ;
            x1 += dx1 ;
            y1 += dy1 ;
        } else {
            x1 += dx2 ;
            y1 += dy2 ;
        }
    }

}

void drawCircle(double r, double xcenter, double ycenter, int color){
    int x = 0, y, xmax, y2, y2_new, ty;
    xmax = (int) (r * 0.70710678); // maximum x at radius/sqrt(2)

    y = r;

    y2 = y * y;
    ty = (2 * y) - 1;
    y2_new = y2;



    for (x = 0; x <= xmax ; x++) {



        if ((y2 - y2_new) >= ty) {
            y2 -= ty;
            y -= 1;
            ty -= 2;
        }


        if(xcenter + x < 800 && ycenter + y < 800 && xcenter + x > 0 && ycenter + y > 0 )
            matrix[(int) (ycenter + y)][(int) (xcenter + x)] = color;

        if(xcenter + x < 800 && ycenter - y < 800 && xcenter + x > 0 && ycenter - y > 0  )
            matrix[(int) (ycenter - y)][(int) (xcenter + x)] = color;

        if(xcenter - x < 800 && ycenter + y < 800 && xcenter - x > 0 && ycenter + y > 0 )
            matrix[(int) (ycenter + y)][(int) (xcenter - x)] = color;

        if(xcenter - x < 800 && ycenter - y < 800 && xcenter - x > 0 && ycenter - y > 0 )
            matrix[(int) (ycenter - y)][(int) (xcenter - x)] = color;

        if(xcenter + y < 800 && ycenter + x < 800 && xcenter + y > 0 && ycenter + x > 0)
            matrix[(int) (ycenter + x)][(int) (xcenter + y)] = color;

        if(xcenter + y < 800 && ycenter - x < 800 && xcenter + y > 0 && ycenter - x > 0)
            matrix[(int) (ycenter - x)][(int) (xcenter + y)] = color;

        if(xcenter - y < 800 && ycenter + x < 800 && xcenter - y > 0 && ycenter + x > 0  )
            matrix[(int) (ycenter + x)][(int) (xcenter - y)] = color;

        if(xcenter - y < 800 && ycenter - x < 800 && xcenter - y > 0 && ycenter - x > 0 )
            matrix[(int) (ycenter - x)][(int) (xcenter - y)] = color;


        y2_new -= (2 * x) - 3;
    }
}

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
    vector<Point> possibleCentroids;
    vector<Point> eliminated;




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

    void addTo(int x, int y)
    {
        Point add;
        add.xval = x;
        add.yval = y;

        possibleCentroids.push_back(add);

    }

    double distance(Point p) {
        return sqrt((p.xval - xval) * (p.xval - xval) + (p.yval - yval) * (p.yval - yval));
    }

    bool compare (Point first)
    {
        unsigned int i=0;
        if (xval == first.xval && yval == first.yval  ) {
            return true;
        }  else { return false;}
    }



};
vector<Point> centroids;
void drawCircle(double r, double xcenter, double ycenter, int color);


const int k = 2;
int n = 0;
int in = 1;

class Node
{
public:
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
    Node* temp = new Node;

    for (int i=0; i<k; i++)
        temp->point[i] = arr[i];


    temp->left = temp->right = nullptr;
    temp->xmin = xmin;
    temp->xmax = xmax;

    temp->ymin = ymin;
    temp->ymax = ymax;
    temp->p.xval = temp->point[0];
    temp->p.yval = temp->point[1];



    return temp;
}

void part4();

void traverse(Node *root, int level, int i);

void drawExtendedHorizontal(int i, int i1, int i2);

void drawExtendedVertical(int i, int i1, int i2);

void generatePoints();

void readPoints();
std::vector<Point> total;

vector<int> nPoints;
vector<double> sumX, sumY;



void kMeansClustering(vector<Point> *pVector);

int getHeight(Node *pNode);

void traverse2(Node *pNode, int trav);

void traverse3(Node *pNode, int num);

void traverse4(Node *pNode, int next);

void traverse5(Node *pNode, int check);

bool foundIn(vector<Point> vector, Point &point);

Node *insertRec(Node *root, int point[], unsigned int depth, int xmin, int xmax, int ymin, int ymax, int i, int i1)
{

    if (root == nullptr)
        return newNode(point, xmin, xmax, ymin, ymax, i, i1);


    unsigned cd = depth % 2;
    //cout << "depth: " << depth << endl;

    if (point[cd] < (root->point[cd])) {
       // cout << "root point left: " << root->point[cd] << endl;
        if(depth%2 == 0)
            root->left = insertRec(root->left, point, depth + 1, xmin, root->point[cd], ymin, ymax, i, i1);
        else
            root->left = insertRec(root->left, point, depth + 1, xmin, xmax, ymin, root->point[cd], i, i1);

    }
    else {
       // cout << "root point right: " << root->point[cd] << endl;
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
                //og test<< 1 << " "<< 0 << " " << 0 << "     ";
                test<< 0 << " "<< 1 << " " << 0 << "     ";
            }

            else if(matrix[n][m] == 3)
            {
                test<< 1 << " "<< 0 << " " << 0 << "     ";
            }
            else if(matrix[n][m] == 2)
            {

                test<< 0 << " "<< 1 << " " << 1 << "     ";
            }

            else if(matrix[n][m] == 5)
            {

                test<< 1 << " "<< 0 << " " << 1 << "     ";
            }

            else if(matrix[n][m] == 6)
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
     Node *root = nullptr;
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
         // cout << "Point2!" << endl;
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
        //drawCircle(2,points[i][0], points[i][1],0 );

    }
    // cout << "Point3!" << endl;


    for (int i = 0; i < 5; ++i) {
        Point p;
        p.xval = points[i][0];
        p.yval = points[i][1];
        cout << "Centroids: "<<p.xval <<"  " << p.yval << endl;
        centroids.push_back(p);
        //root->p.possibleCentroids.push_back(p);
    }

for(int times = 0; times < 50; times++) {

    for(Point put: centroids)
    {
        root->p.possibleCentroids.push_back(put);
    }

    for (int next =0; next < n; next++) {
        traverse4(root, next);
    }

    while (in <= n) {
        //  cout << "Point4!" << endl;
        traverse(root, in, 0);
        in++;
        //cout << "Point5!" << endl;
    }


    for (int j = 0; j < 5; ++j) {
        nPoints.push_back(0);
        sumX.push_back(0.0);
        sumY.push_back(0.0);
    }
    for (int num =0; num < n; num++) {
        //  cout << "Point4!" << endl;
        traverse3(root, num);

        //cout << "Point5!" << endl;
    }

//        root->p.minDist = __DBL_MAX__;
    int x = 0;
    for (vector<Point>::iterator c = begin(centroids); c != end(centroids); ++c) {
        int clusterId = c - begin(centroids);
        Point gen = centroids[x];

//        cout << "npoints: " << nPoints[clusterId] << endl;
//        cout << "sumx: " << sumX[clusterId] << endl;

        c->xval =(int) (sumX[clusterId] / nPoints[clusterId]);
        c->yval =(int) (sumY[clusterId] / nPoints[clusterId]);
//
////            cout << "CXVAL: " << centroids[x].xval << " CYVAL: " << centroids[x].yval<< endl;
////            cout << "GenXVAL: " << gen.xval << " GENYVAL: " << gen.yval << endl;
////            cout << equalsDouble(gen.xval, centroids[x].xval);

        x++;

    }

}

for(Point p: centroids)
{
    cout << "Centroids: "<<p.xval <<"  " << p.yval << endl;
}

    for (int trav =0; trav < n; trav++) {
        //  cout << "Point4!" << endl;
        traverse2(root, trav);

        //cout << "Point5!" << endl;
    }



}

void traverse4(Node *pNode, int next) {

    if (pNode == NULL)
        return;
    if (next == 1) {

       pNode->p.minDist = __DBL_MAX__;


    } else if (next > 1) {
        traverse4(pNode->left, next - 1);
        traverse4(pNode->right, next - 1);
    }
}

void traverse3(Node *pNode, int num) {
    if (pNode == NULL)
        return;
    if (num == 1) {

        int clusterId = pNode->p.cluster;
        nPoints[clusterId] += 1;
        sumX[clusterId] += pNode->p.xval;
        sumY[clusterId] += pNode->p.yval;


    } else if (num > 1) {
        traverse3(pNode->left, num - 1);
        traverse3(pNode->right, num - 1);
    }

}

void traverse2(Node *pNode, int trav) {

    if (pNode == NULL)
        return;
    if (trav == 1) {
        if(pNode->p.cluster == 0)
        {
            drawCircle(2.00, pNode->p.xval, pNode->p.yval, 6);
        }

        if(pNode->p.cluster == 1)
        {
            drawCircle(2.00, pNode->p.xval, pNode->p.yval , 5);
        }

        if(pNode->p.cluster == 2)
        {
            //matrix[(int)(p.yval *800)][(int)(p.xval *800)] = 4;
            drawCircle(2.00, pNode->p.xval, pNode->p.yval, 2);
        }

        if(pNode->p.cluster == 3)
        {

            drawCircle(2.00, pNode->p.xval, pNode->p.yval, 3);
        }

        if(pNode->p.cluster == 4)
        {

            drawCircle(2.00, pNode->p.xval , pNode->p.yval , 4);
        }

        for(Point x: centroids)
        {
            drawCircle(3.00, x.xval , x.yval , 0);
        }


    } else if (trav > 1) {
        traverse2(pNode->left, trav - 1);
        traverse2(pNode->right, trav - 1);
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


void traverse(Node *root, int level, int i) {
    // << "in traverse!" << endl;
    if (root == NULL)
        return;
    if (level == 1) {


        double dist;
        int val = 0;
        for (Point inn: root->p.possibleCentroids) {
            //cout << inn.xval << "  " << inn.yval << endl;

            dist = inn.distance(root->p);
            cout << "dist: " << dist <<endl;
            if (dist < root->p.minDist) {
                root->p.minDist = dist;
                root->p.cluster = val;
            }
            val++;
        }
//        cout << root->p.xval<< "  " << root->p.yval << endl;
//        cout << root->xmin<< "  " << root->xmax << endl;
//        cout << root->ymin<< "  " << root->ymax << endl;
        cout << "Cluster: " << root->p.cluster << endl;
        cout << "Cluster size: " << root->p.possibleCentroids.size() << endl;
    } else if (level > 1) {

        for(int compare1 = 0; compare1 < root->p.possibleCentroids.size(); compare1++) {
            for (int compare2 = 1; compare2 < root->p.possibleCentroids.size()-1; compare2++) {
                

//                if (root->left != NULL) {
//
//                    Point topLeft;
//                    topLeft.xval = root->left->xmin;
//                    topLeft.yval = root->left->ymin;
//
//                    Point topRight;
//                    topRight.xval = root->left->xmax;
//                    topRight.yval = root->left->ymin;
//
//                    Point btmLeft;
//                    btmLeft.xval = root->left->xmin;
//                    btmLeft.yval = root->left->ymax;
//
//                    Point btmRight;
//                    btmRight.xval = root->left->xmax;
//                    btmRight.yval = root->left->ymax;
//
//                    cout << "comp1: " << root->p.possibleCentroids[compare1].xval << "  " << root->p.possibleCentroids[compare1].yval << endl;
//
//                    if (root->p.possibleCentroids[compare1].distance(topLeft) <
//                        root->p.possibleCentroids[compare2].distance(topLeft)) {
//
//                        if (root->p.possibleCentroids[compare1].distance(topRight) <
//                            root->p.possibleCentroids[compare2].distance(topRight)) {
//
//                            if (root->p.possibleCentroids[compare1].distance(btmLeft) <
//                                root->p.possibleCentroids[compare2].distance(btmLeft)) {
//
//                                if (root->p.possibleCentroids[compare1].distance(btmRight) <
//                                    root->p.possibleCentroids[compare2].distance(btmRight)) {
//
//                                    cout << "HERE2" << endl;
//                                    if(!foundIn(root->left->p.eliminated,root->p.possibleCentroids[compare1]) && !foundIn(root->left->p.eliminated,root->p.possibleCentroids[compare2])) {
//                                        root->left->p.eliminated.push_back(root->p.possibleCentroids[compare2]);
//                                    }
//                                }
//
//                            }
//
//                        }
//
//                    }
 //               }

//                //xmin ymax
//                //xmin ymin done
//                //xmax ymin done
////                //xmax ymax
//
//                if (root->left != NULL) {
//                    if (averageX > root->left->xmin && averageY > root->left->ymin) {
//                        if (averageX > root->left->xmax && averageY > root->left->ymin) {
//                            if (averageX > root->left->xmin && averageY > root->left->ymax) {
//                                if (averageX > root->left->xmax && averageY > root->left->ymax) {
//
//
//
//
//
//                                         if(averageY < root->p.possibleCentroids[compare1].yval)
//                                             root->left->p.eliminated.push_back(root->p.possibleCentroids[compare1]);
//
//
//
//
////                                    else {
////                                        root->left->p.eliminated.push_back(root->p.possibleCentroids[compare2]);
////                                    }
//
//
//
//                                }
//
//                            }
//
//                        }
//                    }
//                }
//
//                    if (averageX < root->left->xmin && averageY < root->left->ymin) {
//                        if (averageX < root->left->xmax && averageY < root->left->ymin) {
//                            if (averageX < root->left->xmin && averageY < root->left->ymax) {
//                                if (averageX < root->left->xmax && averageY < root->left->ymax) {
//                                    Point comp;
//                                    comp.xval = root->left->xmin;
//                                    comp.yval = root->left->ymax;
//                                    //cout << "HERE2" << endl;
//                                    if (root->p.possibleCentroids[compare2].distance(comp) <
//                                        root->p.possibleCentroids[compare1].distance(comp)) {
//                                        root->left->p.eliminated.push_back(root->p.possibleCentroids[compare1]);
//
//                                    } else if (root->p.possibleCentroids[compare2].distance(comp) >
//                                               root->p.possibleCentroids[compare1].distance(comp)) {
//                                        root->left->p.eliminated.push_back(root->p.possibleCentroids[compare2]);
//
//                                    }
//
//                                }
//
//                            }
//
//                        }
//
//                    }
//                }
//
//                if (root->right != NULL) {
//                    if (averageX > root->right->xmin && averageY > root->right->ymin) {
//                        if (averageX > root->right->xmax && averageY > root->right->ymin) {
//                            if (averageX > root->right->xmin && averageY > root->right->ymax) {
//                                if (averageX > root->right->xmax && averageY > root->right->ymax) {
//                                    Point comp;
//                                    comp.xval = root->right->xmin;
//                                    comp.yval = root->right->ymax;
//                                    //cout << "HERE3" << endl;
//                                    if (root->p.possibleCentroids[compare2].distance(comp) <
//                                        root->p.possibleCentroids[compare1].distance(comp)) {
//                                        root->right->p.eliminated.push_back(root->p.possibleCentroids[compare1]);
//
//                                    } else {
//                                        root->right->p.eliminated.push_back(root->p.possibleCentroids[compare2]);
//                                    }
//
//
//                                }
//
//                            }
//
//                        }
//                    }
//
//
//                    if (averageX < root->right->xmin && averageY < root->right->ymin) {
//                        if (averageX < root->right->xmax && averageY < root->right->ymin) {
//                            if (averageX < root->right->xmin && averageY < root->right->ymax) {
//                                if (averageX < root->right->xmax && averageY < root->right->ymax) {
//                                    //cout << "HERE4" << endl;
//                                    Point comp;
//                                    comp.xval = root->right->xmin;
//                                    comp.yval = root->right->ymax;
//
//                                    if (root->p.possibleCentroids[compare2].distance(comp) <
//                                        root->p.possibleCentroids[compare1].distance(comp)) {
//                                        root->right->p.eliminated.push_back(root->p.possibleCentroids[compare1]);
//
//                                    } else if (root->p.possibleCentroids[compare2].distance(comp) >
//                                               root->p.possibleCentroids[compare1].distance(comp)) {
//                                        root->right->p.eliminated.push_back(root->p.possibleCentroids[compare2]);
//
//                                    }
//
//                                }
//
//                            }
//
//                        }
//
//                    }
//
//                }
//
            }
        }

    if(root->left != NULL) {
        root->left->p.possibleCentroids.clear();


        for (Point notin: root->p.possibleCentroids) {
            bool notinn = true;
            for (Point el:root->left->p.eliminated) {
                if (notin.xval == el.xval && notin.yval == el.yval) {

                    notinn = false;

                }

            }

            if (notinn) {
                root->left->p.possibleCentroids.push_back(notin);
            }

        }
        root->left->p.eliminated.clear();

    }
    if(root->right != NULL) {
        root->right->p.possibleCentroids.clear();
        for (Point notin: root->p.possibleCentroids) {
            bool notinn = true;
            for (Point el:root->right->p.eliminated) {
                if (notin.xval == el.xval && notin.yval == el.yval) {

                    notinn = false;

                }

            }

            if (notinn) {
                root->right->p.possibleCentroids.push_back(notin);
            }

        }

        root->right->p.eliminated.clear();

    }

        traverse(root->left, level - 1, 0);
        traverse(root->right, level - 1, 1);
    }
}

bool foundIn(vector<Point> input, Point &check) {
    for(Point p : input){
        if(check.compare(p))
            return true;
    }
    return false;
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


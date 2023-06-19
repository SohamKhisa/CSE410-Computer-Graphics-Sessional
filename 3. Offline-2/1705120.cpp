#include <iostream>
#include <stack>
#include <cmath>
#include <fstream>
#include <iomanip> 
#include <time.h>
#include <vector>
#include <algorithm>

#include "bitmap_image.hpp"

#define pi (2*acos(0.0))

using namespace std;


class Point
{
    double x;
    double y;
    double z;
    double w;
public:
    Point() {
        x = 0;
        y = 0;
        z = 0;
        w = 1;
    }
    Point(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
        w = 1;
    }

    Point(double x, double y, double z, double w) {
        this->x = x;
        this->y = y;
        this->z = z;
        this->w = w;
    }

    void setX(double x) {this->x = x;}
    void setY(double y) {this->y = y;}
    void setZ(double z) {this->z = z;}
    void setW(double w) {this->w = w;}

    double getX() {return this->x;}
    double getY() {return this->y;}
    double getZ() {return this->z;}
    double getW() {return this->w;}

    Point operator+(const Point &point);
    Point operator-(const Point &point);
    double operator*(const Point &point);
    Point operator*(const double &val);
    Point operator/(const Point &point);
    Point operator&(const Point &point);
    Point operator=(const Point &point);
    void scaleW();
    void normalize();
};

Point Point::operator=(const Point &point) {
    this->x = point.x;
    this->y = point.y;
    this->z = point.z;
    this->w = point.w;
    return Point(x, y, z, w);
}

Point Point::operator+(const Point &point) {
    Point res;
    res.x = this->x + point.x;
    res.y = this->y + point.y;
    res.z = this->z + point.z;
    return res;
}

Point Point::operator-(const Point &point) {
    Point res;
    res.x = this->x - point.x;
    res.y = this->y - point.y;
    res.z = this->z - point.z;
    return res;
}

double Point::operator*(const Point &point) {
    double res;
    res = this->x * point.x;
    res += this->y * point.y;
    res += this->z * point.z;
    return res;
}

Point Point::operator*(const double &val) {
    Point res;
    res.x = this->x * val;
    res.y = this->y * val;
    res.z = this->z * val;
    return res;   
}

Point Point::operator/(const Point &point) {
    Point res;
    res.x = this->x / point.x;
    res.y = this->y / point.y;
    res.z = this->z / point.z;
    return res;
}

//CrossProduct
Point Point::operator&(const Point &point) {
    Point res;
    res.x = this->y*point.z - this->z*point.y;
    res.y = this->z*point.x - this->x*point.z;
    res.z = this->x*point.y - this->y*point.x;
    return res;
}

void Point::normalize() {
    double m = sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0));
    x /= m;
    y /= m;
    z /= m;
}

void Point::scaleW() {
    x /= w;
    y /= w;
    z /= w;
    w /= w;
}

//=======================================//
//       TRANSFORMATION STRUCTURE        //
//=======================================//
struct Transformation
{
    double **matrix;
    Transformation() {
        matrix = new double*[4];
        for (int i=0; i<4; i++) {
            matrix[i] = new double[4];
        }
        // Identity matrix
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                if (i == j) matrix[i][j] = 1;
                else    matrix[i][j] = 0;
            }
        }
    }

    // Copy constructor
    Transformation(const Transformation &t) {
        this->matrix = new double*[4];
        for (int i=0; i<4; i++) {
            this->matrix[i] = new double[4];
        }
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                this->matrix[i][j] = t.matrix[i][j];
            }
        }
    }

    Transformation(double **matrix) {
        this->matrix = new double*[4];
        for (int i=0; i<4; i++) {
            this->matrix[i] = new double[4];
        }
        for (int i=0; i<4; i++) {
            for(int j=0; j<4; j++) {
                this->matrix[i][j] = matrix[i][j];
            }
        }
    }

    ~Transformation() {
        for(int i=0; i<4; i++){
            delete[] matrix[i];
        }
        delete[] matrix;
    }

    void translationMatrix(double tx, double ty, double tz);
    void scalingMatrix(double sx, double sy, double sz);
    void rotationMatrix(double angle, double ax, double ay, double az);
    void viewMatrix(Point look, Point eye, Point up);
    void projectionMatrix(double fovY, double aspRatio, double near, double far);
    void operator=(const Transformation &t);
};

Transformation multTrans(Transformation trans1, Transformation trans2);
Point multTransPoint(Transformation trans, Point point);

void Transformation::operator=(const Transformation &t) {
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            this->matrix[i][j] = t.matrix[i][j];
        }
    }
}

void Transformation::translationMatrix(double tx, double ty, double tz) {
    matrix[0][3] = tx;
    matrix[1][3] = ty;
    matrix[2][3] = tz;
}

void Transformation::scalingMatrix(double sx, double sy, double sz) {
    matrix[0][0] = sx;
    matrix[1][1] = sy;
    matrix[2][2] = sz;
}

Point rodriguesFormula(Point x, Point a, double theta) {
    Point res = x*cos(theta*pi/180.0) + a*(a*x)*(1-cos(theta*pi/180.0)) + (a&x)*sin(theta*pi/180.0);
    return res;
}

void Transformation::rotationMatrix(double angle, double ax, double ay, double az) {
    Point a(ax, ay, az);    //to be rotated
    a.normalize();

    Point rx = rodriguesFormula(Point(1, 0, 0), a, angle);
    Point ry = rodriguesFormula(Point(0, 1, 0), a, angle);
    Point rz = rodriguesFormula(Point(0, 0, 1), a, angle);

    matrix[0][0] = rx.getX();
    matrix[1][0] = rx.getY();
    matrix[2][0] = rx.getZ();

    matrix[0][1] = ry.getX();
    matrix[1][1] = ry.getY();
    matrix[2][1] = ry.getZ();

    matrix[0][2] = rz.getX();
    matrix[1][2] = rz.getY();
    matrix[2][2] = rz.getZ();
}

void Transformation::viewMatrix(Point look, Point eye, Point up) {
    Point l, r, u;
    l = look - eye;
    l.normalize();
    r = l & up;
    r.normalize();
    u = r & l;

    Transformation translation, rotation;
    translation.translationMatrix(-eye.getX(), -eye.getY(), -eye.getZ());

    rotation.matrix[0][0] = r.getX();
    rotation.matrix[0][1] = r.getY();
    rotation.matrix[0][2] = r.getZ();

    rotation.matrix[1][0] = u.getX();
    rotation.matrix[1][1] = u.getY();
    rotation.matrix[1][2] = u.getZ();

    rotation.matrix[2][0] = -l.getX();
    rotation.matrix[2][1] = -l.getY();
    rotation.matrix[2][2] = -l.getZ();

    Transformation view;
    view = multTrans(rotation, translation);
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++)
            matrix[i][j] = view.matrix[i][j];
    }
}


void Transformation::projectionMatrix(double fovY, double aspRatio, double near, double far) {
    double fovX = fovY * aspRatio;
    double t = near * tan(fovY/2.0 * pi/180.0);
    double r = near * tan(fovX/2.0 * pi/180.0);

    matrix[0][0] = near/r;
    matrix[1][1] = near/t;
    matrix[2][2] = -(far+near)/(far-near);
    matrix[2][3] = -(2.0*far*near)/(far-near);
    matrix[3][2] = -1;
    matrix[3][3] = 0;
}

Transformation multTrans(Transformation trans1, Transformation trans2) {
    double **temp;
    temp = new double*[4];
    for (int i=0; i<4; i++) {
        temp[i] = new double[4];
    }

    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            temp[i][j] = 0;
            for (int k=0; k<4; k++) {
                temp[i][j] += trans1.matrix[i][k] * trans2.matrix[k][j];
            }
        }
    }

    Transformation result(temp);
    return result;
}

Point multTransPoint(Transformation trans, Point point) {
    double *temp;
    temp = new double[4];

    for (int i=0; i<4; i++) {
        temp[i] = 0;
        for(int j=0; j<4; j++) {
            if (j==0)   temp[i] += trans.matrix[i][j] * point.getX();
            else if(j==1)   temp[i] += trans.matrix[i][j] * point.getY();
            else if(j==2)   temp[i] += trans.matrix[i][j] * point.getZ();
            else    temp[i] += trans.matrix[i][j] * point.getW();
        }
    }
    return Point(temp[0], temp[1], temp[2], temp[3]);
}

struct Color
{
    int red;
    int green;
    int blue;

    void setRandomColor() {
        red = rand()%256;
        green = rand()%256;
        blue = rand()%256;
    }
    void setSpecificColor(int red, int green, int blue) {
        this->red = red;
        this->green = green;
        this->blue = blue;
    }
};

struct Triangle
{
    Point points[3];
    Color color; // Color can be primitive/user-defined
    void randomColor() {
        color.setRandomColor();
    }
    void setPoints(Point point1, Point point2, Point point3) {
        points[0] = point1;
        points[1] = point2;
        points[2] = point3;
    }
};

void printPoint(ofstream &output, Point point) {
    output << fixed << setprecision(7) << point.getX() << " " << point.getY() << " " << point.getZ() << endl;
}

int main()
{
    srand(time(0));

    double eyeX, eyeY, eyeZ;
    double lookX, lookY, lookZ;
    double upX, upY, upZ;
    double fovY, aspRatio, near, far;

    ifstream infile;
    ofstream outfile;

    infile.open("scene.txt");
    if(!infile.is_open()) {
        cout << "Failed to open the input file scene.txt" << endl;
        return 1;
    }

    infile >> eyeX >> eyeY >> eyeZ;
    infile >> lookX >> lookY >> lookZ;
    infile >> upX >> upY >> upZ;
    infile >> fovY >> aspRatio >> near >> far;

    //===============Stage-1:MODELING-TRANSFORMATION================//
    string command;
    stack<Transformation> transformStack;
    int ntriangle = 0;
    
    outfile.open("stage1.txt");
    if(!infile.is_open()) {
        cout << "Failed to open the output file stage1.txt" << endl;
        return 1;
    }

    transformStack.push(Transformation());
    while(true)
    {
        infile >> command;
        if (command == "triangle") {
            ntriangle++;
            double p1x, p1y, p1z;
            double p2x, p2y, p2z;
            double p3x, p3y, p3z;

            infile >> p1x >> p1y >> p1z;
            infile >> p2x >> p2y >> p2z;
            infile >> p3x >> p3y >> p3z;

            Point point1(p1x, p1y, p1z);
            Point point2(p2x, p2y, p2z);
            Point point3(p3x, p3y, p3z);

            point1 = multTransPoint(transformStack.top(), point1);
            point2 = multTransPoint(transformStack.top(), point2);
            point3 = multTransPoint(transformStack.top(), point3);
            
            point1.scaleW();
            point2.scaleW();
            point3.scaleW();

            printPoint(outfile, point1);
            printPoint(outfile, point2);
            printPoint(outfile, point3);
            outfile << endl;

        }
        else if (command == "translate") {
            double tx, ty, tz;  
            infile >> tx >> ty >> tz;

            Transformation translate;
            translate.translationMatrix(tx, ty, tz);

            translate = multTrans(transformStack.top(), translate);
            transformStack.pop();
            transformStack.push(translate);

        }
        else if (command == "scale") {
            double sx, sy, sz;
            infile >> sx >> sy >> sz;

            Transformation scaling;
            scaling.scalingMatrix(sx, sy, sz);

            scaling = multTrans(transformStack.top(), scaling);
            transformStack.pop();
            transformStack.push(scaling);

        }
        else if (command == "rotate") {
            double angle, ax, ay, az;
            infile >> angle >> ax >> ay >> az;

            Transformation rotation;
            rotation.rotationMatrix(angle, ax, ay, az);

            rotation = multTrans(transformStack.top(), rotation);
            transformStack.pop();
            transformStack.push(rotation);

        }
        else if (command == "push") {
            transformStack.push(transformStack.top());
        }
        else if (command == "pop") {
            if(transformStack.size() == 0) {
                cout << "Pop operation on empty stack" << endl;
                infile.close();
                outfile.close();
                exit(1);
            }
            transformStack.pop();
        }
        else if (command == "end") {
            break;
        }
        else {
            cout << "Invalid command: " << command << endl;
            infile.close();
            outfile.close();
            exit(1);
        }
    }
    infile.close();
    outfile.close();

    //===============Stage-2:VIEW-TRANSFORMATION================//
    Point look(lookX, lookY, lookZ);
    Point eye(eyeX, eyeY, eyeZ);
    Point up(upX, upY, upZ);

    infile.open("stage1.txt");
    if(!infile.is_open()) {
        cout << "Failed to open the input file scene.txt" << endl;
        return 1;
    }

    outfile.open("stage2.txt");
    if(!infile.is_open()) {
        cout << "Failed to open the output file stage1.txt" << endl;
        return 1;
    }


    Transformation view;
    view.viewMatrix(look, eye, up);
    for (int i=0; i<ntriangle; i++) {
        double p1x, p1y, p1z;
        double p2x, p2y, p2z;
        double p3x, p3y, p3z;

        infile >> p1x >> p1y >> p1z;
        infile >> p2x >> p2y >> p2z;
        infile >> p3x >> p3y >> p3z;

        Point point1(p1x, p1y, p1z);
        Point point2(p2x, p2y, p2z);
        Point point3(p3x, p3y, p3z);

        point1 = multTransPoint(view, point1);
        point2 = multTransPoint(view, point2);
        point3 = multTransPoint(view, point3);
        
        point1.scaleW();
        point2.scaleW();
        point3.scaleW();

        printPoint(outfile, point1);
        printPoint(outfile, point2);
        printPoint(outfile, point3);
        outfile << endl;
    }
    infile.close();
    outfile.close();

    //================Stage-3:PROJECTION-TRANSFORMATION================//
    infile.open("stage2.txt");
    if(!infile.is_open()) {
        cout << "Failed to open the input file stage2.txt" << endl;
        return 1;
    }

    outfile.open("stage3.txt");
    if(!infile.is_open()) {
        cout << "Failed to open the output file stage3.txt" << endl;
        return 1;
    }

    Transformation projection;
    projection.projectionMatrix(fovY, aspRatio, near, far);
    for (int i=0; i<ntriangle; i++) {
        double p1x, p1y, p1z;
        double p2x, p2y, p2z;
        double p3x, p3y, p3z;

        infile >> p1x >> p1y >> p1z;
        infile >> p2x >> p2y >> p2z;
        infile >> p3x >> p3y >> p3z;

        Point point1(p1x, p1y, p1z);
        Point point2(p2x, p2y, p2z);
        Point point3(p3x, p3y, p3z);

        point1 = multTransPoint(projection, point1);
        point2 = multTransPoint(projection, point2);
        point3 = multTransPoint(projection, point3);
        
        point1.scaleW();
        point2.scaleW();
        point3.scaleW();

        printPoint(outfile, point1);
        printPoint(outfile, point2);
        printPoint(outfile, point3);
        outfile << endl;
    }
    infile.close();
    outfile.close();

    //================Stage-4:CLIPPING & SCAN CONVERSION USING Z-BUFFER ALGORITHM================//
    int screenWidth, screenHeight;
    double xLeftlimit, yBottomlimit;
    double xRightlimit, yToplimit;
    double zFront, zRear;

    infile.open("config.txt");
    if(!infile.is_open()) {
        cout << "Failed to open the input file config.txt" << endl;
        return 1;
    }
    
    infile >> screenWidth >> screenHeight;
    infile >> xLeftlimit >> yBottomlimit;
    infile >> zFront >> zRear;
    xRightlimit = -1.0 * xLeftlimit;
    yToplimit = -1.0 * yBottomlimit;
    infile.close();
    
    infile.open("stage3.txt");
    if(!infile.is_open()) {
        cout << "Failed to open the input file stage3.txt" << endl;
        return 1;
    }

    outfile.open("z-buffer.txt");
    if(!infile.is_open()) {
        cout << "Failed to open the output file z-buffer.txt" << endl;
        return 1;
    }

    Triangle triangle[ntriangle];
    for (int i=0; i<ntriangle; i++) {
        double p1x, p1y, p1z;
        double p2x, p2y, p2z;
        double p3x, p3y, p3z;

        infile >> p1x >> p1y >> p1z;
        infile >> p2x >> p2y >> p2z;
        infile >> p3x >> p3y >> p3z;

        Point point1(p1x, p1y, p1z);
        Point point2(p2x, p2y, p2z);
        Point point3(p3x, p3y, p3z);

        triangle[i].setPoints(point1, point2, point3);
        triangle[i].randomColor();
    }

    double dx, dy;
    double topY, leftX;
    double bottomY, rightX;

    dx = (xRightlimit-xLeftlimit)/(float)screenWidth;
    dy = (yToplimit-yBottomlimit)/(float)screenHeight;
    topY = yToplimit - dy/2.0;
    leftX = xLeftlimit + dx/2.0;
    bottomY = yBottomlimit + dy/2.0;
    rightX = xRightlimit - dx/2.0;

    vector<vector<double>> zbuffer((int)screenHeight , vector<double> ((int)screenWidth, zRear));
    bitmap_image frame(screenWidth, screenHeight);
    for (int i=0; i<screenWidth; i++) {
        for (int j=0; j<screenHeight; j++) {
            frame.set_pixel(i, j, 0, 0, 0);
        }
    }


    for (int i=ntriangle-1; i>=0; i--) {
        //determining and clipping the top-bottom scanline
        double maxY = max(triangle[i].points[0].getY(), max(triangle[i].points[1].getY(), triangle[i].points[2].getY()));
        double minY = min(triangle[i].points[0].getY(), min(triangle[i].points[1].getY(), triangle[i].points[2].getY()));

        double topscanline = min(maxY, topY);
        double bottomscanline = max(minY, bottomY);

        //Finding the co-ordinate of the left-right intersection points for each scanline
        for (double scanline=topscanline; scanline>=bottomscanline; scanline-=dy) {
            double xa, xb;
            double za, zb;
            bool firsttime = true;

            for (int j=0; j<3; j++) {
                int ith = j;
                int jth = (j+1)%3;

                //This following condition will be satisfied exacty two times for each scanline
                //Because the horizontal scanline intersects with exactly two edges
                double x1 = triangle[i].points[ith].getX();
                double y1 = triangle[i].points[ith].getY();
                double z1 = triangle[i].points[ith].getZ();

                double x2 = triangle[i].points[jth].getX();
                double y2 = triangle[i].points[jth].getY();
                double z2 = triangle[i].points[jth].getZ();

                if((scanline>y2 && scanline<y1) || (scanline>y1 && scanline<y2)) {
                    double xtemp = (scanline-y2) * (x1-x2)/(y1-y2) + x2;
                    double ztemp = (scanline-y2) * (z1-z2)/(y1-y2) + z2;
                    if(firsttime) {
                        xa = xtemp;
                        za = ztemp;
                        firsttime = false;
                    }
                    else {
                        xb = xtemp;
                        zb = ztemp;
                    }
                }
                else if(scanline==y1 && scanline==y2) {
                    xa = x1;
                    xb = x2;
                    za = z1;
                    zb = z2;
                }
                else if (scanline == y1) {
                    xa = x1;
                    xb = x1;
                    za = z1;
                    zb = z1;
                }
                else if (scanline == y2) {
                    xa = x2;
                    xb = x2;
                    za = z2;
                    zb = z2;
                }
            }

            int leftcol, rightcol;
            double z;

            // Clipping and mapping the left-right intersection points to cth columns
            if(xa > xb) {
                if(xb < leftX)
                    leftcol = 0;
                else
                    leftcol = round((xb-leftX)/dx);
                
                if(xa > rightX)
                    rightcol = screenWidth-1;
                else
                    rightcol = round((xa-leftX)/dx);
                z = zb;
            }
            else {
                if(xa < leftX)
                    leftcol = 0;
                else
                    leftcol = round((xa-leftX)/dx);
                
                if(xb > rightX)
                    rightcol = screenWidth-1;
                else
                    rightcol = round((xb-leftX)/dx);
                z = za;
            }
            //Mapping the scanline to rth row
            int row = round((topY-scanline)/dy);
            double coefZ = (za-zb)/(xa-xb);

            for (int col=leftcol; col<=rightcol; col++) {
                z += coefZ;
                //considering the zs that are between near and far plane;
                if (z>zFront && z<zRear) {
                    zbuffer[row][col] = z;
                    frame.set_pixel(col, row, triangle[i].color.red, triangle[i].color.green, triangle[i].color.blue);
                }
            }
        }
    }
    //Saving the image
    frame.save_image("out.bmp");
    
    //Saving the z values in z-buffer.txt
    for (int i=0; i<screenHeight; i++) {
        for (int j=0; j<screenWidth; j++) {
            if (zbuffer[i][j] < zRear)
                outfile << fixed << setprecision(6) << zbuffer[i][j] << "\t";
        }
        outfile << endl;
    }
    outfile.close();

    for(int i=0; i<zbuffer.size(); i++) {
        zbuffer[i].clear();
    }
    frame.clear();
    return 0;
}
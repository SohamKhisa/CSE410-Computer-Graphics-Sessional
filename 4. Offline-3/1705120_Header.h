#ifndef _1705120_CLASSES_H
#define _1705120_CLASSES_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <algorithm>
#include <vector>

#include <GL/glut.h>

#define pi (2*acos(0.0))
#define EPSILON 0.000001

//***************************************//
//**************   POINT   **************//
//***************************************//
struct Point
{
    double x, y, z;
    Point() {}
    Point(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    Point operator+(const Point &point);
    Point operator-(const Point &point);
    Point operator*(const double value);
    void normalize();
    void setAll(double x, double y, double z);
};
Point Point::operator+(const Point &point) {
    return Point(this->x+point.x, this->y+point.y, this->z+point.z);
}
Point Point::operator-(const Point &point)  {
    return Point(this->x-point.x, this->y-point.y, this->z-point.z);
}
Point Point::operator*(const double value) {
    return Point(this->x*value, this->y*value, this->z*value);
}
void Point:: setAll(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;
}
void Point::normalize() {
    double m = sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0));
    x /= m;
    y /= m;
    z /= m;
}

Point crossProduct(Point vec1, Point vec2)
{
    Point product;
    product.x = vec1.y*vec2.z - vec2.y*vec1.z;
    product.y = -(vec1.x*vec2.z - vec1.z*vec2.x);
    product.z = vec1.x*vec2.y - vec1.y*vec2.x;
    return product;
}
double dotProduct(Point vec1, Point vec2) {
    return (vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z);
}
Point rotation(Point vec, Point rfvec, std::string direction, double rotationangle)
{
    Point product = crossProduct(vec, rfvec);
    Point newpos;
    double dir = 1;
    if(direction=="anticlock")
        dir = -1;

    double rad = dir * (pi/(double)180)*rotationangle;
    newpos.x = vec.x*cos(rad) + product.x*sin(rad);
    newpos.y = vec.y*cos(rad) + product.y*sin(rad);
    newpos.z = vec.z*cos(rad) + product.z*sin(rad);

    return newpos;
}
double getDistance(Point vec1, Point vec2) {
    double x = pow(vec1.x-vec2.x, 2.0);
    double y = pow(vec1.y-vec2.y, 2.0);
    double z = pow(vec1.z-vec2.z, 2.0);
    return sqrt(x+y+z);
}



//***************************************//
//**************   COLOR   **************//
//***************************************//
struct Color
{
    double red, green, blue;
    Color() {}
    Color(double red, double green, double blue) {
        this->red = red;
        this->green = green;
        this->blue = blue;
    }
    Color operator*(const double value);
    Color operator+(const Color &color);
};
Color Color::operator*(const double value) {
    return Color(std::min(this->red*value, 1.0), std::min(this->green*value, 1.0), std::min(this->blue*value, 1.0));
}
Color Color::operator+(const Color &color) {
    return Color(this->red+color.red, this->green+color.green, this->blue+color.blue);
}



//********************************************//
//*************   COEFFICIENTS   *************//
//********************************************//
struct CoEfficients
{
    double ambient, diffuse, specular, recur_reflec_coef;
    CoEfficients() {}
    CoEfficients(double ambient, double diffuse, double specular, double recur_reflec_coef) {
        this->ambient = ambient;
        this->diffuse = diffuse;
        this->specular = specular;
        this->recur_reflec_coef = recur_reflec_coef;
    }
};



//*****************************************//
//****************   RAY   ****************//
//*****************************************//
class Ray
{
    Point start;
    Point dir; // normalize for easier calculations
    //write appropriate constructor
public:
    Ray(Point start, Point dir) {
        this->start = start;
        this->dir = dir;
        this->dir.normalize();
    }
    Point getStart();
    Point getDir();
};
Point Ray::getStart() {
    return start;
}
Point Ray::getDir() {
    return dir;
}


//**************************************//
//************ POINT LIGHT *************//
//**************************************//
class PointLight
{
    Point light_pos;
    Color color;
public:
    PointLight() {}
    PointLight(Point position) {
        light_pos = position;
    }
    void setColor(Color color);
    void draw(bool show);
    Point getLightPos();
    Color getColor();
};
void PointLight::setColor(Color color) {
    this->color = color;
}
void PointLight::draw(bool show) {
    if (show) {
        glColor3f(color.red, color.green, color.blue);
        glPointSize(5);
        glBegin(GL_POINTS);
        {
            glVertex3f(light_pos.x, light_pos.y, light_pos.z);
        }
        glEnd();
    }
}
Point PointLight::getLightPos() {
    return light_pos;
}
Color PointLight::getColor() {
    return color;
}


//**************************************//
//************* SPOT LIGHT *************//
//**************************************//
class SpotLight
{
    PointLight point_light;
    double cutoff;
    Point dir;
public:
    SpotLight(Point point_light) {
        this->point_light = point_light;
    }
    void setColor(Color color);
    void setDir(Point dir);
};
void SpotLight::setColor(Color color) {
    point_light.setColor(color);
}
void SpotLight::setDir(Point dir) {
    this->dir = dir;
}


//**************************************//
//*************   OBJECT   *************//
//**************************************//
class Object
{
protected:
    struct Point ref_point;
    double height, width, length;
    struct Color color;
    CoEfficients coeff;     // ambient, diffuse, specular, reflection-coefficients
    int shine;
public:
    Object() {}
    virtual void draw() {}
    void setRefPoint(Point p);
    void setColor(Color color);
    void setShine(int shine);
    void setCoEfficients(CoEfficients coeff);
    void setHeight(double height);
    void setWidth(double width);
    void setLength(double length);
    Color getColor();
    int getShine();
    CoEfficients getCoEfficient();
    virtual std::string getName() {return "object";}
    virtual double intersect(Ray ray, Color &color, int level) {return -1;};
};
void Object::setRefPoint(Point p) {
    ref_point = p;
}
void Object::setColor(Color color) {
    this->color = color;
}
void Object::setShine(int shine) {
    this->shine = shine;
}
void Object::setCoEfficients(CoEfficients coeff) {
    this->coeff = coeff;
}
void Object::setHeight(double height) {
    this->height = height;
}
void Object::setWidth(double width) {
    this->width = width;
}
void Object::setLength(double length) {
    this->length = length;
}
Color Object::getColor() {
    return color;
}
CoEfficients Object::getCoEfficient() {
    return coeff;
}
int Object::getShine() {
    return shine;
}


//================================= ILLUMINATE ======================================//
void illuminate(Color &color, Object *obj, Point intersection, Point normal, Ray ray);
//===================================================================================//
//================================== REFLECT ========================================//
void reflect(Color &color, Object *obj, Point intersection, Point normal, Ray ray, int level);
//===================================================================================//

//============= LIST, OBJECTS AND VARS FROM 1705120_main.cpp ===============//
std::vector <Object*> objects;
std::vector <PointLight> pointLights;
std::vector <SpotLight> spotLights;

struct Point pos, u, r, l;
int recursion_level, pixels;
//==========================================================================//


//**************************************//
//*************   SPHERE   *************//
//**************************************//
class Sphere: public Object
{
public:
    Sphere(Point center, double radius) {
        ref_point = center;
        length = radius;
    }
    void draw();
    std::string getName() {return "sphere";}
    double intersect(Ray ray, Color &color, int level);
};
void Sphere::draw() {
	struct Point points[100][100];
	int i,j;
	double h,r;
    int stacks = 30;
    int slices = 30;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h = length*sin(((double)i/(double)stacks)*(pi/2));
		r = length*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        glColor3f(color.red, color.green, color.blue);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f((ref_point+points[i][j]).x, (ref_point+points[i][j]).y, (ref_point+points[i][j]).z);
				glVertex3f((ref_point+points[i][j+1]).x, (ref_point+points[i][j+1]).y, (ref_point+points[i][j+1]).z);
				glVertex3f((ref_point+points[i+1][j+1]).x, (ref_point+points[i+1][j+1]).y, (ref_point+points[i+1][j+1]).z);
				glVertex3f((ref_point+points[i+1][j]).x, (ref_point+points[i+1][j]).y, (ref_point+points[i+1][j]).z);
                //lower hemisphere
                glVertex3f((ref_point+points[i][j]).x, (ref_point+points[i][j]).y, (ref_point-points[i][j]).z);
				glVertex3f((ref_point+points[i][j+1]).x, (ref_point+points[i][j+1]).y, (ref_point-points[i][j+1]).z);
				glVertex3f((ref_point+points[i+1][j+1]).x, (ref_point+points[i+1][j+1]).y, (ref_point-points[i+1][j+1]).z);
				glVertex3f((ref_point+points[i+1][j]).x, (ref_point+points[i+1][j]).y, (ref_point-points[i+1][j]).z);
			}glEnd();
		}
	}
}
double Sphere::intersect(Ray ray, Color &color, int level) {
    double a, b, c, d, t1, t2, t;
    a = 1;
    b = 2.0*dotProduct((ray.getStart()-ref_point), ray.getDir());
    c = dotProduct((ray.getStart()-ref_point), (ray.getStart()-ref_point)) - pow(length, 2.0);
    d = pow(b, 2.0) - 4*a*c;
    
    if (d < 0) return -1;
    
    d = sqrt(d);

    t1 = (-b + d)/(2*a);
    t2 = (-b - d)/(2*a);
    if (t1>0 && t2>0)   t = std::min(t1, t2);
    else if (t1<0 && t2<0)  t = -1;
    else if (t1 < 0)    t = t2;
    else if (t2 < 0)    t = t1;

    if (level == 0)   return t;

    //Illumination
    if (fabs(t+1)<0.0001)  t = INT32_MAX;
    Point intersection = ray.getStart()+ray.getDir()*t;
    Point normal = intersection-ref_point;
    normal.normalize();
    if (getDistance(ref_point, pos) < length || fabs(getDistance(ref_point, pos)-length) < 0.00001)
        normal = normal * -1.0;
    illuminate(color, this, intersection, normal, ray);
    if (level < recursion_level)
        reflect(color, this, intersection, normal, ray, level);
    return t;
}



//**************************************//
//************   TRIANGLE   ************//
//**************************************//
class Triangle: public Object
{
    Point point1, point2, point3;
public:
    Triangle(Point point1, Point point2, Point point3) {
        this->point1 = point1;
        this->point2 = point2;
        this->point3 = point3;
    }
    void draw();
    std::string getName() {return "triangle";};
    double intersect(Ray ray, Color &color, int level);
};
void Triangle::draw() {
    glColor3f(color.red, color.green, color.blue);
    glBegin(GL_TRIANGLES);
    {
        glVertex3f(point1.x, point1.y, point1.z);
        glVertex3f(point2.x, point2.y, point2.z);
        glVertex3f(point3.x, point3.y, point3.z);
    }
    glEnd();
}
double Triangle::intersect(Ray ray, Color &color, int level) {
    double deterA, deterBeta, deterGamma, deterT;
    
    deterA = (point1.x-point2.x) * ((point1.y-point3.y)*ray.getDir().z - (point1.z-point3.z)*ray.getDir().y);
    deterA += -(point1.x-point3.x) * ((point1.y-point2.y)*ray.getDir().z - (point1.z-point2.z)*ray.getDir().y);
    deterA += ray.getDir().x * ((point1.y-point2.y)*(point1.z-point3.z) - (point1.z-point2.z)*(point1.y-point3.y));

    deterBeta = (point1.x-ray.getStart().x) * ((point1.y-point3.y)*ray.getDir().z - (point1.z-point3.z)*ray.getDir().y);
    deterBeta += -(point1.x-point3.x) * ((point1.y-ray.getStart().y)*ray.getDir().z - (point1.z-ray.getStart().z)*ray.getDir().y);
    deterBeta += ray.getDir().x * ((point1.y-ray.getStart().y)*(point1.z-point3.z) - (point1.z-ray.getStart().z)*(point1.y-point3.y));

    deterGamma = (point1.x-point2.x) * ((point1.y-ray.getStart().y)*ray.getDir().z - (point1.z-ray.getStart().z)*ray.getDir().y);
    deterGamma += -(point1.x-ray.getStart().x) * ((point1.y-point2.y)*ray.getDir().z - (point1.z-point2.z)*ray.getDir().y);
    deterGamma += ray.getDir().x * ((point1.y-point2.y)*(point1.z-ray.getStart().z) - (point1.z-point2.z)*(point1.y-ray.getStart().y));

    deterT = (point1.x-point2.x) * ((point1.y-point3.y)*(point1.z-ray.getStart().z) - (point1.z-point3.z)*(point1.y-ray.getStart().y));
    deterT += -(point1.x-point3.x) * ((point1.y-point2.y)*(point1.z-ray.getStart().z) - (point1.z-point2.z)*(point1.y-ray.getStart().y));
    deterT += (point1.x-ray.getStart().x) * ((point1.y-point2.y)*(point1.z-point3.z) - (point1.z-point2.z)*(point1.y-point3.y));

    double beta, gamma, t, res;

    if (deterA == 0) res = -1;
    beta = deterBeta/deterA;
    gamma = deterGamma/deterA;
    t = deterT/deterA;

    if (beta+gamma<1 && beta>0 && gamma>0 && t>0) {
        res = t;
    }
    else res = -1;

    if (level == 0) return res;

    if (fabs(res+1)<0.00001) t=INT32_MAX;

    //Illumination
    Point intersection = ray.getStart() + ray.getDir()*t;
    Point normal = crossProduct((point2-point1), (point3-point1));
    normal.normalize();
    if (dotProduct((ray.getDir()*(-1.0)), normal) < 0 || fabs(dotProduct((ray.getDir()*(-1.0)), normal)-0) < 0.00001)
        normal = normal * -1.0;
    illuminate(color, this, intersection, normal, ray);
    if (level < recursion_level)
        reflect(color, this, intersection, normal, ray, level);
    return res;
}


//==========================================//
//================ GENERAL =================//
//==========================================//
class QuadricSurface: public Object
{
    double a, b, c, d, e, f, g, h, i, j;
public:
    QuadricSurface() {}
    QuadricSurface(double a, double b, double c, double d, double e, double f, double g, double h, double i, double j) {
        this->a = a;
        this->b = b;
        this->c = c;
        this->d = d;
        this->e = e;
        this->f = f;
        this->g = g;
        this->h = h;
        this->i = i;
        this->j = j;
    }
    double intersect(Ray ray, Color &color, int level);
};
double QuadricSurface::intersect(Ray ray, Color &color, int level) {
    double a, b, c, t1, t2, t=INT32_MAX; // ax^2 + bx + c = 0, solve x
    a = this->a*pow(ray.getDir().x, 2.0) + this->b*pow(ray.getDir().y, 2.0) + this->c*pow(ray.getDir().z, 2.0);
    a += d*ray.getDir().x*ray.getDir().y + e*ray.getDir().y*ray.getDir().z + f*ray.getDir().x*ray.getDir().z;

    b = 2*this->a*ray.getStart().x*ray.getDir().x + 2*this->b*ray.getStart().y*ray.getDir().y + 2*this->c*ray.getStart().z*ray.getDir().z;
    b += d*ray.getDir().x*ray.getStart().y + d*ray.getStart().x*ray.getDir().y + e*ray.getDir().y*ray.getStart().z;
    b += e*ray.getStart().y*ray.getDir().z + f*ray.getDir().x*ray.getStart().z + f*ray.getStart().x*ray.getDir().z;
    b += g*ray.getDir().x + h*ray.getDir().y + i*ray.getDir().z;

    c = this->a*pow(ray.getStart().x, 2.0) + this->b*pow(ray.getStart().y, 2.0) + this->c*pow(ray.getStart().z, 2.0);
    c += d*ray.getStart().x*ray.getStart().y + e*ray.getStart().y*ray.getStart().z + f*ray.getStart().x*ray.getStart().z;
    c += g*ray.getStart().x + h*ray.getStart().y + i*ray.getStart().z + j;

    double discrim = b*b-4*a*c;
    if(discrim < 0) return -1;
    discrim = sqrt(discrim);

    t1 = (-b + discrim)/(2.0*a);
    t2 = (-b - discrim)/(2.0*a);
    
    if (a==0 && b!=0) t=-c/b;
    else if (t2 > 0)    t = t2;
    else if (t1 > 0)    t = t1;

    // Clipping
    if (t2 > 0) {
        Point intersection = ray.getStart() + ray.getDir()*t2 - ref_point;
        if ((length!=0.0 && (ref_point.x>intersection.x || ref_point.x+length<intersection.x)) || (width!=0.0 && (ref_point.y>intersection.y || ref_point.y+width<intersection.y)) || (height!=0.0 && (ref_point.z>intersection.z || ref_point.z+height<intersection.z)))
            t2 = INT32_MAX;
    }
    if (t1 > 0) {
        Point intersection = ray.getStart() + ray.getDir()*t1 - ref_point;
        if ((length!=0.0 && (ref_point.x>intersection.x || ref_point.x+length<intersection.x)) || (width!=0.0 && (ref_point.y>intersection.y || ref_point.y+width<intersection.y)) || (height!=0.0 && (ref_point.z>intersection.z || ref_point.z+height<intersection.z)))
            t1 = INT32_MAX;
    }
    t = std::min(t1, t2);
    if (fabs(t+1)<0.00001 || t<-1)    t = INT32_MAX;
    if (level==0) return t;

    // Illumination
    Point normal(a, b, c);
    normal.normalize();
    if (dotProduct(ray.getDir()*(-1.0), normal)<0 || fabs(dotProduct(ray.getDir()*(-1.0), normal)-0)<0.00001)
        normal = normal*(-1.0);
    Point intersection = ray.getStart() + ray.getDir()*t - ref_point;
    illuminate(color, this, intersection, normal, ray);
    if (level < recursion_level)
       reflect(color, this, intersection, normal, ray, level);
    return t;
}


//*************************************//
//************    FLOOR    ************//
//*************************************//
class Floor: public Object
{
public:
    Floor(double floorWidth, double tileWidth){
        ref_point = Point(-floorWidth/2.0, -floorWidth/2.0, 0);
        length = tileWidth;
        width = floorWidth;
    }
    void draw();
    double getFloorWidth();
    double getTileWidth();
    std::string getName() {return "Floor";}
    double intersect(Ray ray, Color &color, int level);
};
double Floor::getFloorWidth() {
    return width;
}
double Floor::getTileWidth() {
    return length;
}
void Floor::draw() {
    int ntiles = width/length;
    for (int i=0; i<ntiles; i++) {
        for (int j=0; j<ntiles; j++) {
            glColor3f((i+j)%2, (i+j)%2, (i+j)%2);
            glBegin(GL_QUADS);
            {
                glVertex3f(ref_point.x+i*length, ref_point.y+j*length, ref_point.z);
                glVertex3f(ref_point.x+(i+1)*length, ref_point.y+j*length, ref_point.z);
                glVertex3f(ref_point.x+(i+1)*length, ref_point.y+(j+1)*length, ref_point.z);
                glVertex3f(ref_point.x+i*length, ref_point.y+(j+1)*length, ref_point.z);
            }
            glEnd();
        }
    }
}
double Floor::intersect(Ray ray, Color &color, int level) {
    
    Point normal(0, 0, 1);
    double t = INT32_MAX;

    if(dotProduct(normal, pos) == 0.0)
        normal = normal * (-1.0);
    if (dotProduct(normal, ray.getDir()) == 0.0)
        return -1;
    
    t = (-1.0)*(dotProduct(normal, ray.getStart()))/(dotProduct(normal, ray.getDir()));
    Point intersection = ray.getStart() + ray.getDir()*t;
    if (t>0 && t<INT32_MAX) {
        int tx = (intersection-ref_point).x/length;
        int ty = (intersection-ref_point).y/length;

        double c = (tx+ty)%2;
        this->color = Color(c, c, c);
        // Check if the intersecting floor is within range
        // The ref_point is -floorwith/2, -floorwidth/2, 0
        if (!((intersection.x>ref_point.x && intersection.x<-ref_point.x) && (intersection.y>ref_point.y && intersection.y<-ref_point.y)))
            t = -1;
    }

    if (level == 0) return t;
    
    //Illumination
    if (fabs(t+1) < 0.00001)    t = INT32_MAX;
    intersection = ray.getStart() + ray.getDir()*t;
    illuminate(color, this, intersection, normal, ray);
    if (level < recursion_level)
        reflect(color, this, intersection, normal, ray, level);
    return t;
}


void illuminate(Color &color, Object *obj, Point intersection, Point normal, Ray ray)
{
    color = obj->getColor() * obj->getCoEfficient().ambient;
    
    for (auto light: pointLights) {
        Ray incident(light.getLightPos(), intersection-light.getLightPos());
        double tMin = INT32_MAX;

        for (auto object: objects) {
            double t = object->intersect(incident, color, 0);
            if (t>0 & t<tMin) tMin = t;
        }

        Point shadowpoint = incident.getStart() + incident.getDir()*tMin;
        double dist1 = getDistance(incident.getStart(), intersection) - EPSILON;
        double dist2 = getDistance(incident.getStart(), shadowpoint) - EPSILON;

        // intersection point not in the shadow
        if (dist1 < dist2 || fabs(dist1-dist2)<0.00001) {
            Point reflectdir = incident.getDir() - normal*(dotProduct(normal, incident.getDir())*2.0);
            Ray reflectray(intersection, reflectdir);
            double phong = dotProduct(ray.getDir()*(-1.0), reflectray.getDir());
            double lambert = dotProduct(normal, incident.getDir()*(-1.0));

            color.red += light.getColor().red * obj->getColor().red * obj->getCoEfficient().specular * pow(std::max(phong, 0.0), obj->getShine());
            color.green += light.getColor().green * obj->getColor().green * obj->getCoEfficient().specular * pow(std::max(phong, 0.0), obj->getShine());
            color.blue += light.getColor().blue * obj->getColor().blue * obj->getCoEfficient().specular * pow(std::max(phong, 0.0), obj->getShine());

            color.red += light.getColor().red * obj->getColor().red * obj->getCoEfficient().diffuse * std::max(lambert, 0.0);
            color.green += light.getColor().green * obj->getColor().green * obj->getCoEfficient().diffuse * std::max(lambert, 0.0);
            color.blue += light.getColor().blue * obj->getColor().blue * obj->getCoEfficient().diffuse * std::max(lambert, 0.0);
        }
    }
    //printf("color after: r:%f, g:%f, b:%f\n", obj->getColor().red, obj->getColor().green, obj->getColor().blue);
}

void reflect(Color &color, Object *obj, Point intersection, Point normal, Ray ray, int level)
{
    Point reflection = ray.getDir() - normal*dotProduct(normal, ray.getDir())*2.0;
    reflection.normalize();
    Ray reflectray(intersection+reflection, reflection);
    Color reflectcolor;

    int nearest = -1;
	double t, tMin = INT32_MAX;
    for (int i=0; i<objects.size(); i++) {
        t = objects.at(i)->intersect(reflectray, color, 0);
        if (t>0 && t<tMin) {
            nearest = i;
            tMin = t;
        }
    }
    if(nearest != -1)    t = objects.at(nearest)->intersect(reflectray, reflectcolor, level+1);

    color = color + reflectcolor * obj->getCoEfficient().recur_reflec_coef;

    if (color.red > 1)  color.red = 1;
    else if (color.red < 0) color.red = 0;
    if (color.green > 1)    color.green = 1;
    else if (color.green < 0)  color.green = 0;
    if (color.blue > 1) color.blue = 1;
    else if (color.blue < 0)    color.blue = 0;
}

#endif

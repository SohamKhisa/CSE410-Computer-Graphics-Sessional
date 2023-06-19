#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>

#include <GL/glut.h>
#include "1705120_Header.h"
#include "bitmap_image.hpp"

using namespace std;

#define pi (2*acos(0.0))

double angle;
double viewAngle;
double rotationangle;

// Object and variables for offline-3 ray tracing
extern vector <Object*> objects;
extern vector <PointLight> pointLights;
extern vector <SpotLight> spotLights;
ifstream input;
int bitmapImageCount = 10;

//pos => position of the camera
//u, r, l => unit vectors and all perpendicular to each other
extern struct Point pos, u, r, l;
extern int recursion_level, pixels;
int windowHeight=500, windowWidth=500;

//***************************************************//
//************ RAY TRACING IMPLEMENTATION ***********//
//***************************************************//

//===================================================//
//=================== LOAD-DATA =====================//
//===================================================//
void loadData() {
    string command;
    Object *temp;
    double red, green, blue;
    double ambient, diffuse, specular, recur_reflec_coef;
    double shininess;
    int nobject, npointlight, nspotlight;

    input.open("scene_test.txt");
    if(!input.is_open()) {
        printf("Failed to open the input file scene.txt\n");
        return;
    }

    temp = new Floor(1000, 20);
    temp->setColor(Color(0, 0, 0));
    temp->setCoEfficients(CoEfficients(0.35, 0.25, 0.3, 0.27));
    temp->setShine(3);
    objects.push_back(temp);
    input >> recursion_level;
    input >> pixels;
    input >> nobject;

    while (nobject)
    {
        input >> command;
        if (command == "sphere") {
            double xcenter, ycenter, zcenter;
            double radius;
            
            input >> xcenter >> ycenter >> zcenter;
            input >> radius;
            input >> red >> green >> blue;
            input >> ambient >> diffuse >> specular >> recur_reflec_coef;
            input >> shininess;

            temp = new Sphere(Point(xcenter, ycenter, zcenter), radius);
            temp->setColor(Color(red, green, blue));
            temp->setCoEfficients(CoEfficients(ambient, diffuse, specular, recur_reflec_coef));
            temp->setShine(shininess);

            objects.push_back(temp);
        }
        else if (command == "triangle") {
			Point point1, point2, point3;

			input >> point1.x >> point1.y >> point1.z;
			input >> point2.x >> point2.y >> point2.z;
			input >> point3.x >> point3.y >> point3.z;
			input >> red >> green >> blue;
			input >> ambient >> diffuse >> specular >> recur_reflec_coef;
			input >> shininess;
			
			temp = new Triangle(point1, point2, point3);
			temp->setColor(Color(red, green, blue));
			temp->setCoEfficients(CoEfficients(ambient, diffuse, specular, recur_reflec_coef));
            temp->setShine(shininess);

			objects.push_back(temp);
        }
        else if (command == "general") {
			double A, B, C, D, E, F, G, H, I, J;
			double length, width, height;
			Point cube_ref_point;

			input >> A >> B >> C >> D >> E >> F >> G >> H >> I >> J;
			input >> cube_ref_point.x >> cube_ref_point.y >> cube_ref_point.z;
			input >> length >> width >> height;
			input >> red >> green >> blue;
			input >> ambient >> diffuse >> specular >> recur_reflec_coef;
			input >> shininess;

			temp = new QuadricSurface(A, B, C, D, E, F, G, H, I, J);
			temp->setRefPoint(cube_ref_point);
			temp->setColor(Color(red, green, blue));
			temp->setHeight(height);
			temp->setWidth(width);
			temp->setLength(length);
			temp->setCoEfficients(CoEfficients(ambient, diffuse, specular, recur_reflec_coef));
			temp->setShine(shininess);

			objects.push_back(temp);
        }
        nobject--;
    }

	input >> npointlight;
    while (npointlight) {
        double xpos, ypos, zpos;
		input >> xpos >> ypos >> zpos;
		input >> red >> green >> blue;

        PointLight pointlight(Point(xpos, ypos, zpos));
		pointlight.setColor(Color(red, green, blue));
        pointLights.push_back(pointlight);
        npointlight--;
    }

	input >> nspotlight;
	while(nspotlight) {
		Point source, dir;
		double cutoff;

		input >> source.x >> source.y >> source.z;
		input >> red >> green >> blue;
		input >> dir.x >> dir.y >> dir.z;
		input >> cutoff;

		SpotLight spotLight(source);
		spotLight.setColor(Color(red, green, blue));
		spotLight.setDir(dir);
		spotLights.push_back(spotLight);
		nspotlight--;
	}

    // Here take input of spotlight objects
	input.close();
}

//===================================================//
//==================== CAPTURE ======================//
//===================================================//
void capture() {
	bitmap_image frame(pixels, pixels);
    for (int i=0; i<pixels; i++) {
        for (int j=0; j<pixels; j++) {
            frame.set_pixel(i, j, 0, 0, 0);
        }
    }
	double planeDistance = (windowHeight)/(2.0*tan((viewAngle/2.0)*(pi/180.0)));
	Point topleft = pos + l*planeDistance - r*(windowWidth/2.0) + u*(windowHeight/2.0);		//Can create problem

	double du = (double)windowWidth/pixels;
	double dv = (double)windowHeight/pixels;

	topleft = topleft + r*(du*0.5) - u*(dv*0.5);
	int nearest = -1;
	double t, tMin = INT32_MAX;
	
	for (int i=0; i<pixels; i++) {
		for (int j=0; j<pixels; j++) {
			Point curPixel = topleft + r*((double)i*du) - u*((double)j*dv);
			Ray ray(pos, curPixel-pos);
			Color color;
			
			for (int k=0; k<objects.size(); k++) {
				
				t = objects.at(k)->intersect(ray, color, 0);
				if(t>0 && t<tMin) {
					nearest = k;
					tMin = t;
				}
			}

			if (nearest != -1) {
				t = objects.at(nearest)->intersect(ray, color, 1);
				frame.set_pixel(i, j, color.red*255, color.green*255, color.blue*255);
			}
			nearest = -1;
			tMin = INT32_MAX;
		}
	}
	std::stringstream BitmapCount;
	bitmapImageCount++;
    BitmapCount << bitmapImageCount;
	frame.save_image("output_" + BitmapCount.str() + ".bmp");
	frame.clear();
}



void keyboardListener(unsigned char key, int x,int y){
	switch(key){
		case '0':
			capture();
		case '1':
            r = rotation(r, u, "anticlock", rotationangle);
            l = rotation(l, u, "anticlock", rotationangle);
			break;
        case '2':
			r = rotation(r, u, "clock", rotationangle);
			l = rotation(l, u, "clock", rotationangle);
			break;
        case '3':
			u = rotation(u, r, "anticlock", rotationangle);
            l = rotation(l, r, "anticlock", rotationangle);
			break;
        case '4':
			u = rotation(u, r, "clock", rotationangle);
            l = rotation(l, r, "clock", rotationangle);
			break;
        case '5':
			r = rotation(r, l, "anticlock", rotationangle);
			u = rotation(u, l, "anticlock", rotationangle);
			break;
        case '6':
            r = rotation(r, l, "clock", rotationangle);
			u = rotation(u, l, "clock", rotationangle);
		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow keys
            pos = pos - l*3;
			break;
		case GLUT_KEY_UP:		// up arrow key
			pos = pos + l*3;
			break;

		case GLUT_KEY_RIGHT:
			pos = pos + r*3;
			break;
		case GLUT_KEY_LEFT:
			pos = pos - r*3;
			break;

		case GLUT_KEY_PAGE_UP:
		    pos = pos + u;
			break;
		case GLUT_KEY_PAGE_DOWN:
		    pos = pos - u;
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}

void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				
			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}


void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	gluLookAt(pos.x,pos.y,pos.z,	pos.x+l.x,pos.y+l.y,pos.z+l.z,	u.x,u.y,u.z);    //-> 0,0,200,  0,0,0,  u.x,u.y,u.z


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects
	for (auto o: objects) {
		if(o->getName() != "general")
			o->draw();
	}
	for (auto p: pointLights) {
		p.draw(true);
	}

    //glColor3f(1,0,0);
    //drawSquare(10);

    //drawSS();
    //drawOneEightSphere(70);
    //drawOneForthCylinder(30, 150);

    //drawCone(20,50,24);

	//drawSphere(30,24,20);




	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	angle=0;
	viewAngle = 80;

    pos.setAll(150, 170, 50);
    u.setAll(0, 0, 1);
    r.setAll(-1/(double)sqrt(2), 1/(double)sqrt(2), 0);
    l.setAll(-1/(double)sqrt(2), -1/(double)sqrt(2), 0);
    rotationangle = 3;

	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(viewAngle,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(windowHeight, windowWidth);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("Ray-Tracing");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);


    loadData();
	glutMainLoop();		//The main loop of OpenGL

	objects.clear();
	pointLights.clear();
	spotLights.clear();
	return 0;
}


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>

//#include <windows.h>
#include <GL/glut.h>

#define pi (2*acos(0.0))

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;
double rotationangle;

struct point
{
	double x,y,z;
};

void assignVal(struct point *p, double x, double y, double z)
{
    p->x = x;
    p->y = y;
    p->z = z;
}

void addition(struct point *p1, struct point *p2)
{
    p1->x += p2->x;
    p1->y += p2->y;
    p1->z += p2->z;
}

void subtraction(struct point *p1, struct point *p2)
{
    p1->x -= p2->x;
    p1->y -= p2->y;
    p1->z -= p2->z;
}

struct point crossProduct(struct point vec1, struct point vec2)
{
    struct point product;
    product.x = vec1.y*vec2.z - vec2.y*vec1.z;
    product.y = -(vec1.x*vec2.z - vec1.z*vec2.x);
    product.z = vec1.x*vec2.y - vec1.y*vec2.x;
    return product;
}

struct point rotation(struct point vec, struct point rfvec, std::string direction)
{
    struct point product = crossProduct(vec, rfvec);
    struct point newpos;
    double dir = 1;
    if(direction=="anticlock")
        dir = -1;

    double rad = dir * (pi/(double)180)*rotationangle;
    newpos.x = vec.x*cos(rad) + product.x*sin(rad);
    newpos.y = vec.y*cos(rad) + product.y*sin(rad);
    newpos.z = vec.z*cos(rad) + product.z*sin(rad);

    return newpos;
}
//pos => position of the camera
//u, r, l => unit vectors and all perpendicular to each other
struct point pos, u, r, l;
double cradius, cheight, a;



void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);

			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
		}glEnd();
	}
}


void drawGrid()
{
	int i;
	if(drawgrid==1)
	{
		glColor3f(0.6, 0.6, 0.6);	//grey
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -90, 0);
				glVertex3f(i*10,  90, 0);

				//lines parallel to X-axis
				glVertex3f(-90, i*10, 0);
				glVertex3f( 90, i*10, 0);
			}
		}glEnd();
	}
}

void drawSquare(double a)
{
    //glColor3f(1.0,0.0,0.0);
	glBegin(GL_QUADS);{
		glVertex3f( a, a,2);
		glVertex3f( a,-a,2);
		glVertex3f(-a,-a,2);
		glVertex3f(-a, a,2);
	}glEnd();
}


void drawCircle(double radius,int segments)
{
    int i;
    struct point points[100];
    glColor3f(0.7,0.7,0.7);
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

void drawCone(double radius,double height,int segments)
{
    int i;
    double shade;
    struct point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw triangles using generated points
    for(i=0;i<segments;i++)
    {
        //create shading effect
        if(i<segments/2)shade=2*(double)i/(double)segments;
        else shade=2*(1.0-(double)i/(double)segments);
        glColor3f(shade,shade,shade);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(0,0,height);
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}


void drawSphere(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
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
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();
		}
	}
}

void drawOneEightSphere(double radius)
{
    glColor3f (1, 0, 0);
    struct point points[100][100];
	int i,j;
	double h,r;
	int stacks = 90, slices = 90;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h = radius*sin(((double)i/(double)stacks)*(pi/2));
		r = radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*(pi/2));
			points[i][j].y=r*sin(((double)j/(double)slices)*(pi/2));
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
			}glEnd();
		}
	}
}

void drawOneForthCylinder(double radius, double height)
{
    glColor3f(0, 1, 0);
    int slice = 500;
    point points[501];

    //generate points
    for(int i=0; i<=slice; i++)
    {
        points[i].x = radius*cos(((double)i/(double)slice)*(pi/2));
        points[i].y = radius*sin(((double)i/(double)slice)*(pi/2));
    }

    for(int i=0; i<slice; i++)
    {
        glBegin(GL_QUADS);
        {
            glVertex3f(points[i].x, points[i].y, 0);
            glVertex3f(points[i+1].x, points[i+1].y, 0);
            glVertex3f(points[i+1].x, points[i+1].y, height);
            glVertex3f(points[i].x, points[i].y, height);
        }
        glEnd();
    }
}

void drawCubicSquare(double edge)
{
    glColor3f(1, 1, 1);
    glBegin(GL_QUADS);
    {
        glVertex3f(edge/(double)2, 0, edge/(double)2);
        glVertex3f(-edge/(double)2, 0, edge/(double)2);
        glVertex3f(-edge/(double)2, 0, -edge/(double)2);
        glVertex3f(edge/(double)2, 0, -edge/(double)2);
    }
    glEnd();
}

void drawCubicSphere()
{
    /// 12 -> 1/4 * Cylinders
    glPushMatrix();
    {
        //1
        glTranslatef(cheight/(double)2, cheight/(double)2, -cheight/(double)2);
        drawOneForthCylinder(cradius, cheight);
    }
    glPopMatrix();
    //2
    glPushMatrix();
    {
        glTranslatef(cheight/(double)2, cheight/(double)2, -cheight/(double)2);
        glTranslatef(0, 0, cheight);
        glRotatef(-90, 0, 1, 0);
        drawOneForthCylinder(cradius, cheight);
    }
    glPopMatrix();

    //3
    glPushMatrix();
    {
        glTranslatef(cheight/(double)2, cheight/(double)2, -cheight/(double)2);
        glTranslatef(0, 0, cheight);
        glRotatef(90, 1, 0, 0);
        drawOneForthCylinder(cradius, cheight);
    }
    glPopMatrix();

    //4
    glPushMatrix();
    {
        glTranslatef(cheight/(double)2, cheight/(double)2, -cheight/(double)2);
        glRotatef(-90, 1, 0, 0);
        glRotatef(-90, 0, 1, 0);
        drawOneForthCylinder(cradius, cheight);
    }
    glPopMatrix();

    //5
    glPushMatrix();
    {
        glTranslatef(cheight/(double)2, cheight/(double)2, -cheight/(double)2);
        glRotatef(90, 0, 1, 0);
        glRotatef(90, 1, 0, 0);
        drawOneForthCylinder(cradius, cheight);
    }
    glPopMatrix();

    //6
    glPushMatrix();
    {
        glTranslatef(cheight/(double)2, cheight/(double)2, -cheight/(double)2);
        glTranslatef(-cheight, 0, 0);
        glRotatef(90, 0, 0, 1);
        drawOneForthCylinder(cradius, cheight);
    }
    glPopMatrix();

    //7
    glPushMatrix();
    {
        glTranslatef(cheight/(double)2, cheight/(double)2, -cheight/(double)2);
        glTranslatef(-cheight, 0, 0);
        glTranslatef(0, 0, cheight);
        glRotatef(-90, 0, 1, 0);
        glRotatef(90, 1, 0, 0);
        drawOneForthCylinder(cradius, cheight);
    }
    glPopMatrix();

    //8
    glPushMatrix();
    {
        glTranslatef(cheight/(double)2, cheight/(double)2, -cheight/(double)2);
        glTranslatef(-cheight, 0, 0);
        glRotatef(180, 0, 1, 0);
        glRotatef(90, 1, 0, 0);
        drawOneForthCylinder(cradius, cheight);
    }
    glPopMatrix();

    //9
    glPushMatrix();
    {
        glTranslatef(cheight/(double)2, cheight/(double)2, -cheight/(double)2);
        glTranslatef(-cheight, -cheight, 0);
        glRotatef(180, 0, 0, 1);
        drawOneForthCylinder(cradius, cheight);
    }
    glPopMatrix();

    //10
    glPushMatrix();
    {
        glTranslatef(cheight/(double)2, cheight/(double)2, -cheight/(double)2);
        glTranslatef(0, -cheight, 0);
        glTranslatef(0, 0, cheight);
        glRotatef(90, 1, 0, 0);
        glRotatef(-90, 0, 1, 0);
        drawOneForthCylinder(cradius, cheight);
    }
    glPopMatrix();

    //11
    glPushMatrix();
    {
        glTranslatef(cheight/(double)2, cheight/(double)2, -cheight/(double)2);
        glTranslatef(0, -cheight, 0);
        glRotatef(-180, 1, 0, 0);
        glRotatef(-90, 0, 1, 0);
        drawOneForthCylinder(cradius, cheight);
    }
    glPopMatrix();

    //12
    glPushMatrix();
    {
        glTranslatef(cheight/(double)2, cheight/(double)2, -cheight/(double)2);
        glTranslatef(0, -cheight, 0);
        glRotatef(-90, 0, 0, 1);
        drawOneForthCylinder(cradius, cheight);
    }
    glPopMatrix();


    /// 8 -> 1/8 * Sphere
    glPushMatrix();
    {
        glTranslatef(cheight/(double)2, cheight/(double)2, -cheight/(double)2);
        //a
        glPushMatrix();
        {
            glTranslatef(0, 0, cheight);
            drawOneEightSphere(cradius);
        }
        glPopMatrix();

        //b
        glPushMatrix();
        {
            glRotatef(-180, 1, 1, 0);
            drawOneEightSphere(cradius);
        }
        glPopMatrix();

        //c
        glPushMatrix();
        {
            glTranslatef(-cheight, 0, 0);
            glRotatef(90, 0, 0, 1);
            glTranslatef(0, 0, cheight);
            drawOneEightSphere(cradius);
        }
        glPopMatrix();

        //d
        glPushMatrix();
        {
            glTranslatef(-cheight, 0, 0);
            glRotatef(90, 0, 0, 1);
            glRotatef(-180, 1, 1, 0);
            drawOneEightSphere(cradius);
        }
        glPopMatrix();

        //e
        glPushMatrix();
        {
            glTranslatef(-cheight, -cheight, cheight);
            glRotatef(180, 0, 0, 1);
            drawOneEightSphere(cradius);
        }
        glPopMatrix();

        //f
        glPushMatrix();
        {
            glTranslatef(-cheight, -cheight, 0);
            glRotatef(-90, 0, 0, 1);
            glRotatef(180, 1, 0, 0);
            drawOneEightSphere(cradius);
        }
        glPopMatrix();

        //g
        glPushMatrix();
        {
            glTranslatef(0, -cheight, cheight);
            glRotatef(-90, 0, 0, 1);
            drawOneEightSphere(cradius);
        }
        glPopMatrix();

        //h
        glPushMatrix();
        {
            glTranslatef(0, -cheight, 0);
            glRotatef(180, 1, 0, 0);
            drawOneEightSphere(cradius);
        }
        glPopMatrix();
    }
    glPopMatrix();

    /// Squared screen
    // i
    glPushMatrix();
    {
        glTranslatef(0, (cheight/(double)2)+cradius, 0);
        drawCubicSquare(cheight);
    }
    glPopMatrix();

    // ii
    glPushMatrix();
    {
        glTranslatef(-(cheight/(double)2)-cradius, 0, 0);
        glRotatef(-90, 0, 0, 1);
        drawCubicSquare(cheight);
    }
    glPopMatrix();

    // iii
    glPushMatrix();
    {
        glTranslatef((cheight/(double)2)+cradius, 0, 0);
        glRotatef(90, 0, 0, 1);
        drawCubicSquare(cheight);
    }
    glPopMatrix();

    // iv
    glPushMatrix();
    {
        glTranslatef(0, -(cheight/(double)2)-cradius, 0);
        drawCubicSquare(cheight);
    }
    glPopMatrix();

    // v
    glPushMatrix();
    {
        glTranslatef(0, 0, (cheight/(double)2)+cradius);
        glRotatef(90, 1, 0, 0);
        drawCubicSquare(cheight);
    }
    glPopMatrix();

    //vi
    glPushMatrix();
    {
        glTranslatef(0, 0, -(cheight/(double)2)-cradius);
        glRotatef(90, 1, 0, 0);
        drawCubicSquare(cheight);
    }
    glPopMatrix();
}


struct point path[500];
void elipseRotate(double a, double b)
{
    double theta = 30;
    for (int i=0; i<500; i++) {
        path[i].x = a*cos((i/(double)500) * 2 * pi);
        path[i].y = a*sin((i/(double)500) * 2 * pi);
        path[i].z = 0;
    }
}


void drawSS()
{
    glColor3f(1,0,0);
    drawSphere(30, 90, 90);

    glPushMatrix();
    {
        glRotatef(5*angle,0,0,1);
        glTranslatef(70,0,0);
        glRotatef(10*angle,0,0,1);
        glColor3f(0,1,0);
        drawSphere(5, 90, 90);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glRotatef(3*angle,0,0,1);
        glTranslatef(130,0,0);
        glRotatef(3*angle,0,0,1);
        glColor3f(0,0,1);
        drawSphere(15, 90, 90);
    }
    glPopMatrix();

    glRotatef(2*angle,0,0,1);
    glTranslatef(250,0,0);
    glRotatef(2*angle,0,0,1);
    glColor3f(0,1,0);
    drawSphere(10, 90, 90);

    glPushMatrix();
    {
        glRotatef(5*angle,0,0,1);
        glTranslatef(30,0,0);
        glRotatef(8*angle,0,0,1);
        glColor3f(0,0,1);
        drawSphere(3, 90, 90);
    }
    glPopMatrix();
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case '1':
            r = rotation(r, u, "anticlock");
            l = rotation(l, u, "anticlock");
			break;
        case '2':
			r = rotation(r, u, "clock");
			l = rotation(l, u, "clock");
			break;
        case '3':
			u = rotation(u, r, "anticlock");
            l = rotation(l, r, "anticlock");
			break;
        case '4':
			u = rotation(u, r, "clock");
            l = rotation(l, r, "clock");
			break;
        case '5':
			r = rotation(r, l, "anticlock");
			u = rotation(u, l, "anticlock");
			break;
        case '6':
            r = rotation(r, l, "clock");
			u = rotation(u, l, "clock");
		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
            subtraction(&pos, &l);
			break;
		case GLUT_KEY_UP:		// up arrow key
			addition(&pos, &l);
			break;

		case GLUT_KEY_RIGHT:
			addition(&pos, &r);
			break;
		case GLUT_KEY_LEFT:
			subtraction(&pos, &r);
			break;

		case GLUT_KEY_PAGE_UP:
		    addition(&pos, &u);
			break;
		case GLUT_KEY_PAGE_DOWN:
		    subtraction(&pos, &u);
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
            cradius += 1;
            if(cradius>=a/(double)2)  cradius=a/(double)2;
            cheight = a-2*cradius;
            if(cheight<=0) cheight = 0;
			break;
		case GLUT_KEY_END:
		    cradius -= 1;
            if(cradius<=0)  cradius=0;
		    cheight = a-2*cradius;
		    if(cheight>=60)  cheight = 60;
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
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

	drawAxes();
	drawGrid();

    //glColor3f(1,0,0);
    //drawSquare(10);

    drawSS();

    // drawCubicSphere();
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
	drawgrid=0;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;

	assignVal(&pos, 100, 100, 0);
	assignVal(&u, 0, 0, 1);
    assignVal(&r, -1/(double)sqrt(2), 1/(double)sqrt(2), 0);
    assignVal(&l, -1/(double)sqrt(2), -1/(double)sqrt(2), 0);
    rotationangle = 3;

    a = 60;
    cradius = 10;
    cheight = a-2*cradius;

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
	gluPerspective(80,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("Task-1 & Task-2");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
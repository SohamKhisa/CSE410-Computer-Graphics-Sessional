#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//#include <windows.h>
#include <GL/glut.h>

#define pi (2*acos(0.0))

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;
double dis;
double rotangle;
double centerotangle;
double wheelradius;

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

struct point pos;

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
			for(i=-20;i<=20;i++){

				//if(i==0)
					//continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -200, 0);
				glVertex3f(i*10,  200, 0);

				//lines parallel to X-axis
				glVertex3f(-200, i*10, 0);
				glVertex3f( 200, i*10, 0);
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
        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
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


void drawRectangle(double p)
{
    glColor3f(160/(double)255, 160/(double)255, 160/(double)255);
    glBegin(GL_QUADS);
    {
        glVertex3f(0, p, wheelradius);
        glVertex3f(0, -p, wheelradius);
        glVertex3f(0, -p, -wheelradius);
        glVertex3f(0, p, -wheelradius);
    }
    glEnd();
}


void drawRim(double width)
{
    struct point points[500];
    int segment = 499;
    double h = width;
    for(int i=0; i<=segment; i++) {
        points[i].x = wheelradius*cos((double)i/(double)segment * (pi/2));
        points[i].y = h;
        points[i].z = wheelradius*sin((double)i/(double)segment * (pi/2));
    }
    for(int j=0; j<segment; j++)
    {
        glColor3f(j/(double)segment, j/(double)segment, j/(double)segment);
        glBegin(GL_QUADS);{
            //1st quadrant
            glVertex3f(points[j].x, -points[j].y, points[j].z);
            glVertex3f(points[j+1].x, -points[j+1].y, points[j+1].z);
            glVertex3f(points[j+1].x, points[j+1].y, points[j+1].z);
            glVertex3f(points[j].x, points[j].y, points[j].z);

            //2nd
            glVertex3f(-points[j].x, -points[j].y, points[j].z);
            glVertex3f(-points[j+1].x, -points[j+1].y, points[j+1].z);
            glVertex3f(-points[j+1].x, points[j+1].y, points[j+1].z);
            glVertex3f(-points[j].x, points[j].y, points[j].z);

            //3rd
            glVertex3f(-points[j].x, -points[j].y, -points[j].z);
            glVertex3f(-points[j+1].x, -points[j+1].y, -points[j+1].z);
            glVertex3f(-points[j+1].x, points[j+1].y, -points[j+1].z);
            glVertex3f(-points[j].x, points[j].y, -points[j].z);

            //4th
            glVertex3f(points[j].x, -points[j].y, -points[j].z);
            glVertex3f(points[j+1].x, -points[j+1].y, -points[j+1].z);
            glVertex3f(points[j+1].x, points[j+1].y, -points[j+1].z);
            glVertex3f(points[j].x, points[j].y, -points[j].z);
        }glEnd();
    }
}

void drawWheel(double halfwidth)
{
    glPushMatrix();
    {
        drawRim(halfwidth);
    }
    glPopMatrix();

    glPushMatrix();
    {
        drawRectangle(halfwidth);
    }
    glPopMatrix();


    glPushMatrix();
    {
        glRotatef(-90, 0, 1, 0);
        drawRectangle(halfwidth);
    }
    glPopMatrix();
}


void drawPanzer(double halfwidth, double len, double widt)
{
	glPushMatrix();
	{
		glTranslatef(0, 0, wheelradius);
    glPushMatrix();
    {
        glTranslatef(len/(double)2, widt/(double)2, 0);
		glTranslatef(pos.x, pos.y, pos.z);  // forward/backward movement
        glRotatef(rotangle, 0, 0, 1);       // wheel rotation a/d
        glRotatef(centerotangle, 0, 1, 0); // wheel rotation while moving forward or backward
        drawWheel(halfwidth);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslatef(len/(double)2, -widt/(double)2, 0);
		glTranslatef(pos.x, pos.y, pos.z);  // forward/backward movement
        glRotatef(rotangle, 0, 0, 1);       // wheel rotation a/d
        glRotatef(centerotangle, 0, 1, 0); // wheel rotation while moving forward or backward
        drawWheel(halfwidth);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslatef(-len/(double)2, widt/(double)2, 0);
		glTranslatef(pos.x, pos.y, pos.z);  // forward/backward movement
        glRotatef(rotangle, 0, 0, 1);       // wheel rotation a/d
        glRotatef(centerotangle, 0, 1, 0); // wheel rotation while moving forward or backward
        drawWheel(halfwidth);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslatef(-len/(double)2, -widt/(double)2, 0);
		glTranslatef(pos.x, pos.y, pos.z);  // forward/backward movement
        glRotatef(rotangle, 0, 0, 1);       // wheel rotation a/d
        glRotatef(centerotangle, 0, 1, 0); // wheel rotation while moving forward or backward
        drawWheel(halfwidth);
    }
    glPopMatrix();

    glPushMatrix();
    {
		glTranslatef(pos.x, pos.y, pos.z);
        glBegin(GL_QUADS);
        {
            glVertex3f(len/(double)2, widt/(double)2, wheelradius);
            glVertex3f(len/(double)2, -widt/(double)2, wheelradius);
            glVertex3f(-len/(double)2, -widt/(double)2, wheelradius);
            glVertex3f(-len/(double)2, widt/(double)2, wheelradius);
        }
        glEnd();
    }
    glPopMatrix();

	glPushMatrix();
	{
		glTranslatef(pos.x, pos.y, pos.z);
		glRotatef(rotangle, 0, 0, 1);
		glColor3f(1, 0, 0);
		glBegin(GL_LINES);
		{
			glVertex3f(0, 0, wheelradius);
			glVertex3f(0, 0, wheelradius+10);
			glVertex3f(0, 0, wheelradius+10);
			glVertex3f(len+10, 0, wheelradius+10);
		}
		glEnd();
	}
	glPopMatrix();
	}
	glPopMatrix();
}


void moveforward()
{
    pos.x = pos.x - dis*cos((rotangle)*pi/(double)180);
    pos.y = pos.y - dis*sin((rotangle)*pi/(double)180);
}

void movebackward()
{
    pos.x = pos.x + dis*cos((rotangle)*pi/(double)180);
    pos.y = pos.y + dis*sin((rotangle)*pi/(double)180);
}

void wheelmovement(double halfwidth)
{
    glPushMatrix();
    {
        glTranslatef(0, 0, wheelradius);    // positioning the wheel on the xy axis
        glTranslatef(pos.x, pos.y, pos.z);  // forward/backward movement
        glRotatef(rotangle, 0, 0, 1);       // wheel rotation a/d
        glRotatef(centerotangle, 0, 1, 0); // wheel rotation while moving forward or backward
        drawWheel(halfwidth);
    }
    glPopMatrix();
}

void drawSS()
{
    glColor3f(1,0,0);
    drawSquare(20);

    glRotatef(angle,0,0,1);
    glTranslatef(110,0,0);
    glRotatef(2*angle,0,0,1);
    glColor3f(0,1,0);
    drawSquare(15);

    glPushMatrix();
    {
        glRotatef(angle,0,0,1);
        glTranslatef(60,0,0);
        glRotatef(2*angle,0,0,1);
        glColor3f(0,0,1);
        drawSquare(10);
    }
    glPopMatrix();

    glRotatef(3*angle,0,0,1);
    glTranslatef(40,0,0);
    glRotatef(4*angle,0,0,1);
    glColor3f(1,1,0);
    drawSquare(5);
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case 'w':
		    moveforward();
		    centerotangle -= dis/(2*pi*wheelradius) * 360;
			break;
        case 's':
            movebackward();
            centerotangle += dis/(2*pi*wheelradius) * 360;
            break;
        case 'a':
            rotangle+=2;
            break;
        case 'd':
            rotangle-=2;
            break;
		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			cameraHeight -= 3.0;
			break;
		case GLUT_KEY_UP:		// up arrow key
			cameraHeight += 3.0;
			break;

		case GLUT_KEY_RIGHT:
			cameraAngle += 0.03;
			break;
		case GLUT_KEY_LEFT:
			cameraAngle -= 0.03;
			break;

		case GLUT_KEY_PAGE_UP:
			break;
		case GLUT_KEY_PAGE_DOWN:
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
	gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	//gluLookAt(0,0,200,	0,0,0,	0,1,0);


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

    // wheelmovement(5);
    drawPanzer(5, 60, 60);
    
    //drawSS();

    //drawCircle(30,24);

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
	drawgrid=1;
	drawaxes=0;
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;

	rotangle=0;
	dis=2;
	centerotangle=0;
	wheelradius=20;

	assignVal(&pos, 0, 0, 0);

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
	glutInitWindowSize(600, 600);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

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
#include "CEarth.h"

//旋转地球及坐标系用变量
GLfloat rotx = -90.0f;
GLfloat roty = 0.0f;
////////////////////////////////////////////////////////////////////////////////
//class CEarth
CEarth::CEarth()
{
	radius = 63.71f;//地球半径单位为100 km
	rotateSpeed = 0.5f;//设置旋转速度

}

AUX_RGBImageRec* earthMap[2];	//用于存贮地球纹理

AUX_RGBImageRec* LoadBMP(char* Filename)				// Loads A Bitmap Image
{
	FILE* File = NULL;									// File Handle

	if (!Filename)										// Make Sure A Filename Was Given
	{
		return NULL;									// If Not Return NULL
	}

	fopen_s(&File, Filename, "r");							// Check To See If The File Exists

	if (File)											// Does The File Exist?
	{
		fclose(File);									// Close The Handle
		return auxDIBImageLoadA(Filename);				// Load The Bitmap And Return A Pointer
	}

	return NULL;										// If Load Failed Return NULL
}

void CEarth::CreateTexture()
{
	if (earthMap[1] = LoadBMP("./data/earthmap.bmp"))
	{
		int sizeX = earthMap[1]->sizeX;
		int sizeY = earthMap[1]->sizeY;

		glGenTextures(1, &texture[0]);					// Create Three Textures

														// Create Linear Filtered Texture
		glBindTexture(GL_TEXTURE_2D, texture[0]);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexImage2D(GL_TEXTURE_2D, 0, 3, sizeX, sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, earthMap[1]->data);

		free(earthMap[1]->data);
	}
}

void CEarth::render()
{

	glPushMatrix();
	glRotatef(rotx, 1.0f, 0.0f, 0.0f);

	gluQuadricDrawStyle(earthModel, GL_FILL);
	gluQuadricTexture(earthModel, GL_TRUE);
	gluSphere(earthModel, radius, 100, 100);


	glPopMatrix();

	draw3DAxis();
	drawEquator();

}

void CEarth::draw3DAxis()
{

	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glEnable(GL_BLEND);
	glEnable(GL_LINE_SMOOTH);

	float lwidth = 2.0;	//set line width
	float axisLength = 300.0f;

	glBegin(GL_LINES);
	glLineWidth(lwidth);

	//请注意，下面按地固坐标系绘制的三轴，分别对应OpenGL的坐标轴有所不同
	glColor3f(1.0f, 0.0f, 0.0f);
	//此处绘制的X轴，在OpenGL的Z轴上
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, axisLength);

	glColor3f(0.0f, 1.0f, 0.0f);
	//此处绘制的Y轴，在OpenGL的X轴上
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(axisLength, 0.0f, 0.0f);

	glColor3f(0.0f, 0.0f, 1.0f);
	//此处绘制的Z轴，在OpenGL的Y轴上
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, axisLength, 0.0f);

	glEnd();
	glPopAttrib();

}

void CEarth::drawEquator()
{
	glPushAttrib(GL_ALL_ATTRIB_BITS);

	glEnable(GL_BLEND);
	glEnable(GL_LINE_SMOOTH);

	float lwidth = 2.0;	//set line width
	float slid = 100.0f;

	glLineWidth(lwidth);
	glColor3f(1.0f, 1.0f, 0.0f);

	glBegin(GL_LINE_LOOP);
	float rd = radius + 0.1f;//增加一点赤道半径，从地表突出显示
	for (int i = 0; i < slid; ++i)
		glVertex3f(rd * cos(2.0f * PI / slid * i), 0.0f, rd * sin(2.0f * PI / slid * i));

	glEnd();

	glPopAttrib();

}
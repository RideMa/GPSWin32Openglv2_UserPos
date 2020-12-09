#include "CSatellite.h"

////////////////////////////////////////////////////////////////////////////////
CSatellite::CSatellite()
{
	curPosition.x = -1.0;
	curPosition.y = -1.0;
	curPosition.z = -1.0;

	sateModel = gluNewQuadric();
	datahasRead = false;
}

void CSatellite::render(float* clr)
{
	if (datahasRead)
	{
		orbit.ephe = ephemeris.e;
		orbit.ephrA = ephemeris.rA;

		orbit.render(clr);

		glPushMatrix();

		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glColor3f(1.0f, 0.0f, 1.0f);

		//注意：此处OpenGL坐标系与地球坐标系的X、Y、Z轴需要交换，参照函数draw3DAxis中的坐标交换顺序
		glTranslatef(curPosition.y / SCALE, curPosition.z / SCALE, curPosition.x / SCALE);

		gluQuadricDrawStyle(sateModel, GL_FILL);
		gluQuadricNormals(sateModel, GLU_SMOOTH);
		gluSphere(sateModel, 6, 10, 10);
		glPopAttrib();

		glPopMatrix();
	}
	else
	{
		//下面仅为未读取数据时，绘制一个默认的卫星及轨道

		orbit.render(clr);

		glPushMatrix();

		glPushAttrib(GL_ALL_ATTRIB_BITS);

		//如果卫星位置没有读取计算，则因使用按事先给定的位置轨道数据，故需进行旋转
		glRotatef(60.0, 0.0f, 1.0f, 0.0f);
		glRotatef(55.0, 0.0f, 0.0f, 1.0f);
		glColor3f(0.0f, 1.0f, 1.0f);

		//注意：此处OpenGL坐标系与地球坐标系的X、Y、Z轴需要交换，参照函数draw3DAxis中的坐标交换顺序
		glTranslatef(curPosition.y / SCALE, curPosition.z / SCALE, curPosition.x / SCALE);

		gluQuadricDrawStyle(sateModel, GL_FILL);
		gluQuadricNormals(sateModel, GLU_SMOOTH);
		gluSphere(sateModel, 6, 10, 10);
		glPopAttrib();

		glPopMatrix();
	}
}

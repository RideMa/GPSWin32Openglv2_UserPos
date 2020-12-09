#include "COrbit.h"

////////////////////////////////////////////////////////////////////////////////
COrbit::COrbit()
{
	I = 0;
	L = 0;

	ephe = 0.006164086284f;//轨道偏心率值
	ephrA = 5153.69580f;//轨道半径平方根

}


static GLfloat color[][3] = { { 0.0, 1.0, 0.0 },{ 1.0, 0.0, 0.0 },{ 0.2, 0.5, 1.0 },
{ 1.0, 1.0, 0.0 },{ 0.0, 1.0, 1.0 },{ 1.0, 0.7, 0.5 }, };

void COrbit::render(float* clr)
{
	//如果轨道倾角与赤经值均为0，说明程序尚未进行计算，否则我们认为已经完成计算
	if (I == 0 && L == 0)//如果轨道未完成真实轨道计算，则依实际比例与大致参数绘制一个示意性的椭圆轨道
	{
		double ra = ephrA * ephrA / SCALE; //计算轨道半径，并将其单位由米转化为公里，再缩小100倍，即尺度单位为 "百公里"
		double rb = ra * (1 - ephe);
		double cf = sqrt(ra * ra - rb * rb);//椭圆焦点
		renderwithParm(cf, ra, rb, 55.0f, 60.0f, clr);
	}
	else
	{
		renderwithParm2(I, L, clr);
	}

}


void COrbit::renderwithParm(double cf, double ra, double rb, float i, float L, float* clr)
{
	float slid = 100.0f;
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glEnable(GL_BLEND);
	glEnable(GL_LINE_SMOOTH);

	float lwidth = 1.0;	//set line width
	glLineWidth(lwidth);
	glColor3f(clr[0], clr[1], clr[2]);

	glPushMatrix();
	glRotatef(L, 0.0f, 1.0f, 0.0f);
	glRotatef(i, 0.0f, 0.0f, 1.0f);

	glBegin(GL_LINE_LOOP);
	for (int i = 0; i < slid; ++i)
	{
		float sx = ra * cos(2.0f * PI / slid * i) + 0.2f - cf;
		float sy = rb * sin(2.0f * PI / slid * i) + 0.2f;

		glVertex3f(sx, 0.0f, sy);
	}
	glEnd();

	glPopMatrix();
	glPopAttrib();

}

void COrbit::renderwithParm2(float i, float L, float* clr)
{
	double A = ephrA * ephrA;
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glEnable(GL_BLEND);
	glEnable(GL_LINE_SMOOTH);

	float lwidth = 1.0;	//set line width
	glLineWidth(lwidth);
	glColor3f(clr[0], clr[1], clr[2]);

	glBegin(GL_LINE_LOOP);
	for (float En = 0; En < 2 * PI; En += 0.05f)
	{
		float r = A * (1 - ephe * cos(En));
		float x = r * cos(En);
		float y = r * sin(En);

		float px = x * cos(L) - y * cos(i) * sin(L);
		float py = x * sin(L) + y * cos(i) * cos(L);
		float pz = y * sin(i);

		glVertex3f(py / SCALE, pz / SCALE, px / SCALE);
	}
	glEnd();

	glPopAttrib();
}

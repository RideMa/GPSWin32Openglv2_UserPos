#include "COrbit.h"

////////////////////////////////////////////////////////////////////////////////
COrbit::COrbit()
{
	I = 0;
	L = 0;

	ephe = 0.006164086284f;//���ƫ����ֵ
	ephrA = 5153.69580f;//����뾶ƽ����

}


static GLfloat color[][3] = { { 0.0, 1.0, 0.0 },{ 1.0, 0.0, 0.0 },{ 0.2, 0.5, 1.0 },
{ 1.0, 1.0, 0.0 },{ 0.0, 1.0, 1.0 },{ 1.0, 0.7, 0.5 }, };

void COrbit::render(float* clr)
{
	//�����������ྭֵ��Ϊ0��˵��������δ���м��㣬����������Ϊ�Ѿ���ɼ���
	if (I == 0 && L == 0)//������δ�����ʵ������㣬����ʵ�ʱ�������²�������һ��ʾ���Ե���Բ���
	{
		double ra = ephrA * ephrA / SCALE; //�������뾶�������䵥λ����ת��Ϊ�������С100�������߶ȵ�λΪ "�ٹ���"
		double rb = ra * (1 - ephe);
		double cf = sqrt(ra * ra - rb * rb);//��Բ����
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

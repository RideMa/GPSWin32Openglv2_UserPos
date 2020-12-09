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

		//ע�⣺�˴�OpenGL����ϵ���������ϵ��X��Y��Z����Ҫ���������պ���draw3DAxis�е����꽻��˳��
		glTranslatef(curPosition.y / SCALE, curPosition.z / SCALE, curPosition.x / SCALE);

		gluQuadricDrawStyle(sateModel, GL_FILL);
		gluQuadricNormals(sateModel, GLU_SMOOTH);
		gluSphere(sateModel, 6, 10, 10);
		glPopAttrib();

		glPopMatrix();
	}
	else
	{
		//�����Ϊδ��ȡ����ʱ������һ��Ĭ�ϵ����Ǽ����

		orbit.render(clr);

		glPushMatrix();

		glPushAttrib(GL_ALL_ATTRIB_BITS);

		//�������λ��û�ж�ȡ���㣬����ʹ�ð����ȸ�����λ�ù�����ݣ����������ת
		glRotatef(60.0, 0.0f, 1.0f, 0.0f);
		glRotatef(55.0, 0.0f, 0.0f, 1.0f);
		glColor3f(0.0f, 1.0f, 1.0f);

		//ע�⣺�˴�OpenGL����ϵ���������ϵ��X��Y��Z����Ҫ���������պ���draw3DAxis�е����꽻��˳��
		glTranslatef(curPosition.y / SCALE, curPosition.z / SCALE, curPosition.x / SCALE);

		gluQuadricDrawStyle(sateModel, GL_FILL);
		gluQuadricNormals(sateModel, GLU_SMOOTH);
		gluSphere(sateModel, 6, 10, 10);
		glPopAttrib();

		glPopMatrix();
	}
}

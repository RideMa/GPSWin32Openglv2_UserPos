#pragma once
#include "myWork.h"
#include "COrbit.h"

class CSatellite {
public:
	short id;
	CSatellite();
	double xc, yc;//��ά���ƽ������
	POINT3D curPosition;//��ǰʱ������λ��

	COrbit orbit;
	Ephemeris ephemeris;//������ǰ��������
	bool datahasRead;

	GLUquadric* sateModel;//�������������ģ��

public:
	void SetPos(double a, double b, double c)
	{
		curPosition.x = a; curPosition.y = b; curPosition.z = c;
	}

	//Ϊ��ʹ���û��۲�ʱ�䣬������Ҫ���ú����Ƶ�CUser����
	//void CalculateSatPosition();//calculate satellite position by ephemeris data
	void render(float* clr);

};

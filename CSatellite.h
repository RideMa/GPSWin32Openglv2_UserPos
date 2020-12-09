#pragma once
#include "myWork.h"
#include "COrbit.h"

class CSatellite {
public:
	short id;
	CSatellite();
	double xc, yc;//二维轨道平面坐标
	POINT3D curPosition;//当前时刻卫星位置

	COrbit orbit;
	Ephemeris ephemeris;//存贮当前卫星星历
	bool datahasRead;

	GLUquadric* sateModel;//按球体绘制卫星模型

public:
	void SetPos(double a, double b, double c)
	{
		curPosition.x = a; curPosition.y = b; curPosition.z = c;
	}

	//为了使用用户观测时间，我们需要将该函数移到CUser类中
	//void CalculateSatPosition();//calculate satellite position by ephemeris data
	void render(float* clr);

};

#pragma once
#include "myWork.h"
//new add class
class COrbit
{

public:
	COrbit();

	void renderwithParm(double cf, double ra, double rb, float i, float L, float* clr);
	void renderwithParm2(float i, float L, float* clr);
	void render(float* clr);

public:
	double I;//轨道倾角
	double L;//轨道赤经

	double ephrA;//轨道半径
	double ephe;//轨道偏心率

};



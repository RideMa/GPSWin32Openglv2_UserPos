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
	double I;//������
	double L;//����ྭ

	double ephrA;//����뾶
	double ephe;//���ƫ����

};



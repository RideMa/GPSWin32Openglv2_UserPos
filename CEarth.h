#pragma once
#include "myWork.h"

//new add class
class CEarth {

	GLUquadric* earthModel;

public:
	CEarth();

	void CreateEarth() //initialize earth model
	{
		CreateTexture();
		earthModel = gluNewQuadric();
	}

	float radius;
	bool bRotate;
	float rotateSpeed;

	bool bSwitchCoordinates;//�����жϵع�����ϵ����������ϵ֮����л�

	void render();//���Ƶ���
	void draw3DAxis();//������ά����ϵ
	void drawEquator();//���Ƴ��
	void CreateTexture();//������������

};


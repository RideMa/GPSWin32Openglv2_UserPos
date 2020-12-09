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

	bool bSwitchCoordinates;//用于判断地固坐标系与天球坐标系之间的切换

	void render();//绘制地球
	void draw3DAxis();//绘制三维坐标系
	void drawEquator();//绘制赤道
	void CreateTexture();//创建地球纹理

};


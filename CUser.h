#pragma once
#include "myWork.h"
#include "CSatellite.h"

class CUser {

public:
	CUser();
	int satNum;//save satellite number, however, its taken as a constant in current case
	ObserveData observeData;//观测数值

public:
	CSatellite* sate;//定义为指针便于操作，后续可用于多卫星数据处理
	void readEphemeris();//此函数读取星历文件
	bool datahasRead;
	bool bGotStationPos;//判断用户或接收机位置计算是否完成
	bool bDrawsigLines;//判断是否绘制表达卫星信号的直线，以直观可视化接收机观测卫星信号的状态

   //我们将定位过程分两步进行，计算卫星轨道位置和计算接收机位置，因而将原来的
   //void processData();  函数拆分为：
	void processSateData();
	void processUserData();

	void render();

	//计算卫星在轨位置，为计算方便从CSatellite类移到此处，并作为一定修改
	void CalculateSatPosition(CSatellite*);


	//计算接收机位置
	double CalculateUserPosition(double* seudoDis);
	void XYZ2LambdaPhiH(double x, double y, double z, double& lambda, double& phi, double& h);
	void LambdaPhiH2ENU(double dx, double dy, double dz, double lambda, double phi, double h, double& E, double& N, double& U);
	//误差改正函数
	double Clockerr(short sateID);
	double Clockrelative_err(short sateID);
	double Tropospheric_Hopfield_err(short sateID);
	double Ionospheric_Klobuchar_err(short sateID);

	//误差处理总调用函数
	double processErr(short sateID);//for all data process
									//关于坐标转换的函数，作为内联函数放在cpp文件


};

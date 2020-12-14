#pragma once

#include <windows.h>		// Header File For Windows
#include <gl\gl.h>			// Header File For The OpenGL32 Library
#include <gl\glu.h>			// Header File For The GLu32 Library
#include <gl\GLAux.h>		// Header File For The Glaux Library

//该文件头用于解决文本字符集因编译环境不同引起的兼容性问题
//如提示“不能从const char *转换为LPCWSTR”时，可使用 _T()转换字符串
#include <tchar.h>			

#include <math.h>

#include <iostream>
#include <fstream>

//加入Eigen 函数库，用于处理向量与矩阵计算
#include <Eigen/Dense>
using namespace Eigen;

using namespace std;


//define some constant for calculation
#define Omega_dote	7.2921151467e-5	//WGS-84 value of the earths roration rate
#define GM	3.986005e14 //地球质量引力常数
#define PI 3.14159265758f //圆周率
#define F_RELATIVITY -4.442807633e-10 // 相对论误差中，椭圆轨道导致的误差公式里面的常数
#define SCALE 100000 //可视化时，将实际尺度缩小10万倍，如轨道半径由原来的米变为"百公里"
#define LIGHTSPEED 299792458.0

//此处应添加其它常数，如光速，相对论改正常数等

struct Ephemeris {
	double a0, a1, a2;
	double M0;
	double rA;//orbit raduis
	double e;
	double omega;
	double I0;
	double omega_dot;
	double I_dot;
	double Tgd;
	double dn;
	double Crc, Crs;
	double Cuc, Cus;
	double Cic, Cis;
	double w;
	double toe;
	double toc;
	double observt;//observe time
};

struct POINT3D
{
	double x;
	double y;
	double z;
};

//需定义观测数据结构，如：
struct ObserveData {
	double t; //用户观测时间
	POINT3D upos; //接收机位置，或观测位置
	float T, P, Pva; //气象参数：温度、大气压与水汽压
	double Latitude, longitude; //接收机或观测位置处的纬度与经度
	double altitude; //接收机或观测位置处的高程
	double iono_A[4]; //电离层参数，来自导航电文
	double iono_B[4]; //电离层参数，来自导航电文
	double pseudo_range[32]; //伪距观测值,来自测站观测文件
}; 


//new add variables
static unsigned int texture[3];	// 纹理数组



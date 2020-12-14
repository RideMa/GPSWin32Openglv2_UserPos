#pragma once

#include <windows.h>		// Header File For Windows
#include <gl\gl.h>			// Header File For The OpenGL32 Library
#include <gl\glu.h>			// Header File For The GLu32 Library
#include <gl\GLAux.h>		// Header File For The Glaux Library

//���ļ�ͷ���ڽ���ı��ַ�������뻷����ͬ����ļ���������
//����ʾ�����ܴ�const char *ת��ΪLPCWSTR��ʱ����ʹ�� _T()ת���ַ���
#include <tchar.h>			

#include <math.h>

#include <iostream>
#include <fstream>

//����Eigen �����⣬���ڴ���������������
#include <Eigen/Dense>
using namespace Eigen;

using namespace std;


//define some constant for calculation
#define Omega_dote	7.2921151467e-5	//WGS-84 value of the earths roration rate
#define GM	3.986005e14 //����������������
#define PI 3.14159265758f //Բ����
#define F_RELATIVITY -4.442807633e-10 // ���������У���Բ������µ���ʽ����ĳ���
#define SCALE 100000 //���ӻ�ʱ����ʵ�ʳ߶���С10�򱶣������뾶��ԭ�����ױ�Ϊ"�ٹ���"
#define LIGHTSPEED 299792458.0

//�˴�Ӧ�����������������٣�����۸���������

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

//�趨��۲����ݽṹ���磺
struct ObserveData {
	double t; //�û��۲�ʱ��
	POINT3D upos; //���ջ�λ�ã���۲�λ��
	float T, P, Pva; //����������¶ȡ�����ѹ��ˮ��ѹ
	double Latitude, longitude; //���ջ���۲�λ�ô���γ���뾭��
	double altitude; //���ջ���۲�λ�ô��ĸ߳�
	double iono_A[4]; //�������������Ե�������
	double iono_B[4]; //�������������Ե�������
	double pseudo_range[32]; //α��۲�ֵ,���Բ�վ�۲��ļ�
}; 


//new add variables
static unsigned int texture[3];	// ��������



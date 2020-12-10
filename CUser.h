#pragma once
#include "myWork.h"
#include "CSatellite.h"

class CUser {

public:
	CUser();
	int satNum;//save satellite number, however, its taken as a constant in current case
	ObserveData observeData;//�۲���ֵ

public:
	CSatellite* sate;//����Ϊָ����ڲ��������������ڶ��������ݴ���
	void readEphemeris();//�˺�����ȡ�����ļ�
	bool datahasRead;
	bool bGotStationPos;//�ж��û�����ջ�λ�ü����Ƿ����
	bool bDrawsigLines;//�ж��Ƿ���Ʊ�������źŵ�ֱ�ߣ���ֱ�ۿ��ӻ����ջ��۲������źŵ�״̬

   //���ǽ���λ���̷��������У��������ǹ��λ�úͼ�����ջ�λ�ã������ԭ����
   //void processData();  �������Ϊ��
	void processSateData();
	void processUserData();

	void render();

	//���������ڹ�λ�ã�Ϊ���㷽���CSatellite���Ƶ��˴�������Ϊһ���޸�
	void CalculateSatPosition(CSatellite*);


	//������ջ�λ��
	double CalculateUserPosition(double* seudoDis);
	void XYZ2LambdaPhiH(double x, double y, double z, double& lambda, double& phi, double& h);
	void LambdaPhiH2ENU(double dx, double dy, double dz, double lambda, double phi, double h, double& E, double& N, double& U);
	//����������
	double Clockerr(short sateID);
	double Clockrelative_err(short sateID);
	double Tropospheric_Hopfield_err(short sateID);
	double Ionospheric_Klobuchar_err(short sateID);

	//�����ܵ��ú���
	double processErr(short sateID);//for all data process
									//��������ת���ĺ�������Ϊ������������cpp�ļ�


};

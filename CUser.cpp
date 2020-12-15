#include "CUser.h"


CUser::CUser()
{

	satNum = SAT_NUM; //������������ֻ����6�����ǵ�����
	sate = new CSatellite[satNum];

	//��ʼ��α�����
	fstream file;
	file.open(DIST_PATH);
	if (file.is_open())
	{

		double temp = 0;
		const int monthdays[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
		int y, mn, d, h, min, sec;
		file >> y >> mn >> d >> h >> min >> sec;
		// toc
		int sumDay = y * 365;
		for (int i = 0; i < mn - 1; i++)
			sumDay += monthdays[i];

		if (mn > 2) sumDay += y / 4;
		else sumDay += (y - 1) / 4;
		sumDay += (d - 1);
		int ResDays = sumDay % 7;
		observeData.t = ResDays * 86400 + h * 3600 + min * 60 + sec;

		//observeData.pseudo_range = new double[satNum];
		for (int i = 0; i < satNum; ++i)
			file >> observeData.pseudo_range[i] >> temp >> temp >> temp;

		observeData.T = 5 + 273.15;
		observeData.P = 101.3;
		observeData.Pva = 0.6;
		file.close();
	}

	//��ʼ���û���ʼλ��
	observeData.upos = *new POINT3D();
	observeData.upos.x = -2165390.324400646;
	observeData.upos.y = 4380659.056451826;
	observeData.upos.z = 4098945.436700276;

	for (int i = 0; i < satNum; i++)
		sate[i].id = i;

	bDrawsigLines = true;
	bGotStationPos = false;
	datahasRead = false;

}

void CUser::readEphemeris()
{
	if (datahasRead)
		return;


	// ��д�Լ��Ĵ���
	//��ȡn�����ǵ���������
	fstream file(EPHE_PATH);
	if (file.is_open())
	{
		for (int i = 0; i < 4; ++i)
			file >> observeData.iono_A[i];

		for (int i = 0; i < 4; ++i)
			file >> observeData.iono_B[i];

		file >> satNum;

		for (int i = 0; i < satNum; ++i)
		{
			file >> sate[i].id;
			if (sate[i].id >= 0)
			{
				double temp = 0;

				const int monthdays[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
				int y, mn, d, h, min;
				double sec;
				file >> y >> mn >> d >> h >> min >> sec;
				// toc
				int sumDay = y * 365;
				for (int i = 0; i < mn - 1; i++)
					sumDay += monthdays[i];

				if (mn > 2) sumDay += y / 4;
				else sumDay += (y - 1) / 4;
				sumDay += (d - 1);
				int ResDays = sumDay % 7;
				

				Ephemeris ep{};
				if (file.is_open())
				{
					file >> ep.a0 >> ep.a1 >> ep.a2
						>> temp >> ep.Crs >> ep.dn >> ep.M0
						>> ep.Cus >> ep.e >> ep.Cuc >> ep.rA
						>> ep.toe >> ep.Cic >> ep.omega >> ep.Cis
						>> ep.I0 >> ep.Crc >> ep.w >> ep.omega_dot
						>> ep.I_dot >> temp >> temp >> temp
						>> temp >> temp >> ep.Tgd >> temp
						>> ep.toc >> temp;
				}
				sate[i].ephemeris = ep;
				sate[i].ephemeris.toc = ResDays * 86400 + h * 3600 + min * 60 + sec;
				sate[i].ephemeris.observt = observeData.t;
			}
			else exit(0);
			sate[i].datahasRead = true;
		}
		file.close();
	}

	datahasRead = true;//data successfully has read

}

//kepler function get root with Newton itertion
//we start with E0 = M
inline double calculateE(double e, double m, double En1)
{
	return En1 - (En1 - e * sin(En1) - m) / (1 - e * cos(En1));

}

void CUser::CalculateSatPosition(CSatellite* sate)// get coordinate of current satellite
{
	if (!datahasRead)
		return;
	Ephemeris e = sate->ephemeris;
	//����������

	// ������Ľ��ٶ�n
	double n0 = sqrt(GM) / (pow(e.rA, 3));
	double n = n0 + e.dn;

	// ƽ�����M

	double delta_t = e.observt - e.toe;
	double M = e.M0 + n * delta_t;

	// ����ƫ�����E
	double E = M;
	double deltaE = 0.000001, dE = 1;
	while (abs(dE) > abs(deltaE))
	{
		double En = calculateE(e.e, M, E);
		dE = En - E;
		E = En;
	}

	// ����������f
	double f = atan2(sqrt(1 - e.e * e.e) * sin(E), (cos(E) - e.e));

	// ����Ǿ�ȷu���뾶rt�����it
	double raw_u = f + e.w;
	double rt = e.rA * e.rA * (1 - e.e * cos(E));
	double it = e.I0 + e.I_dot * delta_t;

	// ������
	double delta_u = e.Cuc * cos(2 * raw_u) + e.Cus * sin(2 * raw_u);
	double delta_r = e.Crc * cos(2 * raw_u) + e.Crs * sin(2 * raw_u);
	double delta_i = e.Cic * cos(2 * raw_u) + e.Cis * sin(2 * raw_u);

	// ��ȷr��u��i
	double r = rt + delta_r, u = raw_u + delta_u, i = it + delta_i;

	// ƽ��x��y
	double x = r * cos(u), y = r * sin(u);

	// ����L
	double L = e.omega + e.omega_dot * delta_t - Omega_dote * e.observt;

	double X = x * cos(L) - y * cos(i) * sin(L);
	double Y = x * sin(L) + y * cos(i) * cos(L);
	double Z = y * sin(i);

	// �洢
	sate->orbit.I = i;
	sate->orbit.L = L;
	sate->curPosition.x = X;
	sate->curPosition.y = Y;
	sate->curPosition.z = Z;
	return;
}


//�˺���һ��������й��������뽫���ݶ�ȡ��λ�ü���ֿ����ɲ����øú���
void CUser::processSateData()
{
	readEphemeris();//һ���Զ�ȡ������������

	for (int i = 0; i < satNum; i++)
	{
		CalculateSatPosition(&sate[i]);
	}

	datahasRead = true;

}


//����������
double CUser::processErr(short sateID)
{

	//���ʱ��ȡ������ע�ͣ����峣��LIGHTSPEED������ʹ���������


	double dtrop = Tropospheric_Hopfield_err(sateID);//���������

	double dion = Ionospheric_Klobuchar_err(sateID) * LIGHTSPEED;//��������

	double dclk = Clockerr(sateID) * LIGHTSPEED;//ʱ�����

	double drclk = Clockrelative_err(sateID) * LIGHTSPEED;//��������


	double sumerr = 0;
	sumerr = - dtrop - dion + dclk;// +drclk;
	return sumerr;

}

double CUser::Clockerr(short sateID)
{
	double clk_offset = 0.0;
	Ephemeris e = sate[sateID].ephemeris;
	clk_offset = e.a0 + e.a1 * (observeData.t - e.toc) + e.a2 * pow(observeData.t - e.toc, 2) - e.Tgd;
	return clk_offset;
}


double CUser::Clockrelative_err(short sateID)
{

	double clk_relati = 0.0;
	Ephemeris e = sate[sateID].ephemeris;
	double E = e.M0;
	double deltaE = 0.000001, dE = 1;
	while (abs(dE) > abs(deltaE))
	{
		double En = calculateE(e.e, e.M0, E);
		dE = En - E;
		E = En;
	}
	clk_relati = F_RELATIVITY * sate[sateID].ephemeris.rA * sate[sateID].ephemeris.e * sin(E);

	return clk_relati;

}
/* 
* ���´���ʵ��Klobchar����Ƶ��������ģ�ͺ�Hopefield�����������ģ��
* ���ڸ����Ĳ���������ȷ������û��ֱ�ӽ���Ҫ�Ĳ��������CUser���У�ֻ���Ժ�����������ʽд������������
* Klobuchar��Ҫ��
* 		��վ��ά���꣬��������һ����ʽ��ĿǰĬ�ϸ���������άֱ������XYZ���ں�������Ux��Uy��Uz��
* 		������ά���꣬���ｫ���Ƕ������˺�����Ϊ������
* 		����ʱ��������ù۲�ʱ�����3600�ٶ�24��ģ���㣻
* 		�����������������㹫ʽ�в���alpha��beta�����������ṹ���ж�ȡ��
* Hopefield��Ҫ��
* 		��վ��ά���꣬��������һ����ʽ��ĿǰĬ�ϸ���������άֱ������XYZ���ں�������Ux��Uy��Uz��
* 		������ά���꣬���ｫ���Ƕ������˺�����Ϊ������
* 		��վ��ѹ��ˮ��ѹ���¶ȣ����ﶼ��Ϊ�������룻
* ��ν�����������������뻹��Ҫ�������Ĳ������룻
*/

// ���Ƕ�ת��Ϊ����
double Angle2Arc(double angle)
{
	return PI * angle / 180;
}
// ������ת��Ϊ�Ƕ�
double Arc2Angle(double arc)
{
	return 180 * arc / PI;
}

// ����������׼������վ���꣬���ò�վ������֮������Ĳ�ֵ�Ͳ�վ��γ�ȸ߳̾Ϳ���ͨ������˷��õ����
// �������ղ�վ�����ǵ���ά������վ��γ�ȸ̡߳�ENU�������������
void CUser::LambdaPhiH2ENU(double dx, double dy, double dz, double lambda, double phi, double h, double& E, double& N, double& U)
{
	lambda = Angle2Arc(lambda);
	phi = Angle2Arc(phi);
	// ������Ϊ����˷��������������ʽ������������ֱ�����ó˷������
	E = -sin(lambda) * dx + cos(lambda) * dy;
	N = -sin(phi) * cos(lambda) * dx - -sin(phi) * sin(lambda) * dy + cos(phi) * dz;
	U = cos(phi) * cos(lambda) * dx + cos(phi) * sin(lambda) * dy + sin(phi) * dz;
}

// ��XYZ����ת��Ϊ��γ�ȸ߳����꣬����XYZ����;�γ�ȸ̵߳�����
// ʹ����ʦ��PPT�и����ļ��㷽����ʡȥ�˵������鷳
void CUser::XYZ2LambdaPhiH(double x, double y, double z, double& lambda, double& phi, double& h)
{
	// ��Ҫ��������
	double a = 6378137; // ����������볤��
	double e = 8.1819190842622e-2; // ����������ƫ����
	double b = sqrt(a * a * (1 - e * e)); // �̰���
	double ep = sqrt((a * a - b * b) / (b * b)); // �ڶ�ƫ����
	double p = sqrt(x * x + y * y); // ����

	lambda = atan2(y, x); // ���㾭��

	double th = atan2((a * z) , (b * p));// ����
	phi = atan2(z + ep * ep * b * pow(sin(th), 3), p - e * e * a * pow(cos(th), 3)); // ����γ��

	double N = a / sqrt(1 - e * e * pow(sin(phi), 2)); // ����
	h = p / cos(phi) - N; // ����߳�
	lambda = Arc2Angle(lambda);
	phi = Arc2Angle(phi);
}

double CUser::Ionospheric_Klobuchar_err(short sateID)
{
	double Iono_delay;
	CSatellite sate = this->sate[sateID];

	// ���Ǵ��������ȡ�����ģ���������P��A
	double* alpha = (observeData.iono_A), * beta = (observeData.iono_B);

	// ��վ�����������
	double dx = sate.curPosition.x - observeData.upos.x, dy = sate.curPosition.y - observeData.upos.y, dz = sate.curPosition.z - observeData.upos.z;

	// ��γ�ȸ߳�����
	double lambda = observeData.longitude, phi = observeData.Latitude, h = observeData.altitude;

	double UT = lambda * 240 + (int)observeData.t % 86400; // ���㵱ǰ��׼����ʱ

	// ����γ�ȸ߳�����׼��ΪENU���ꡣҲ���ǽ���������ת������վ����ϵ��
	double E = 0, N = 0, U = 0;
	LambdaPhiH2ENU(dx, dy, dz, lambda, phi, h, E, N, U);

	// ��Ҫ�������������еĽǶȶ�ʹ�û����ƣ���Ϊmath.hʹ�õĶ��ǻ����ƣ�����Ҫ��ʱ���ٽ���ת��
	double norm = sqrt(E * E + N * N + U * U); // ģ
	double el = asin(U / norm);    // ���Ǹ߶Ƚ�
	double Azim = atan2(E / norm, N / norm);   // ���̵������λ�÷�λ��
	double EA = 445 / (Arc2Angle(el) + 20) - 4; // ����λ�úͲ�վλ����Ե��ļн�
	double phiN = phi + EA * cos(Azim);   // ���̵����γ��
	double lambdaN = lambda + EA * sin(Azim) / cos(Angle2Arc(phi)); // ���̵���ľ���
	double phiM = Angle2Arc(phiN + 10.07 * cos(Angle2Arc(lambdaN - 288.04))); // ���̵�N�ĵش�γ�� 

	double t = lambda * 240 + (int)observeData.t % 86400;; // ���̵�N�ĵط�ʱ
	if (t > 86400) t -= 86400; // �������������24Сʱ�����ȥ24h������
	else if (t < 0) t += 86400; // ���������С��0����ӦתΪ��ֵ


	// ��������չ�ʽ�е����������������������и�����alpha��beta
	double A = 0, P = 0;
	for (int i = 0; i < 4; ++i)
	{
		A += alpha[i] * pow(phiM, i);
		P += beta[i] * pow(phiM, i);
	}

	// ��������̵��춥������ź�ʱ��
	double Tg = 5e-9 + A * cos((2 * (double)PI / P) * (t - 50400));
	// �ź�·��������춥�������ϵ��
	double secZ = 1 + (96 - Arc2Angle(el)) / 45;
	// �����źŷ���ʱ��
	Iono_delay = Tg * secZ;

	return Iono_delay;
}

double CUser::Tropospheric_Hopfield_err(short sateID)
{
	double Trop_delay;
	CSatellite sate = this->sate[sateID];
	// ��վ�����������
	double dx = sate.curPosition.x - observeData.upos.x, dy = sate.curPosition.y - observeData.upos.y, dz = sate.curPosition.z - observeData.upos.z;

	// ��γ�ȸ߳�����
	double lambda = observeData.longitude, phi = observeData.Latitude, h = observeData.altitude;

	// ����γ�ȸ߳�����׼��ΪENU���ꡣҲ���ǽ���������ת������վ����ϵ��
	double E = 0, N = 0, U = 0;
	LambdaPhiH2ENU(dx, dy, dz, lambda, phi, h, E, N, U);

	double norm = sqrt(E * E + N * N + U * U); // ģ
	double el = Arc2Angle(asin(U / norm)); // ���Ǹ߶Ƚ�

	double hw = 11000; // ʪ�����ֵĸ߶�
	double hd = 40136 + 148.72 * (observeData.T - 273.16); // �������ֵĸ߶�
	double Kd = 155.2e-7 * observeData.P * (hd - h) / observeData.T;  // �������ֵ��ӳ�
	double Kw = 155.2e-7 * 4810 * observeData.Pva * (hw - h) / observeData.T / observeData.T; // ʪ�����ֵ��ӳ�

	// �����ڸ�ʪ������ͶӰ��ɵļӺͣ�Ҳ���Ƕ�������źŵ��ӳ�
	Trop_delay = (Kd / sin(Angle2Arc(sqrt(el * el + 6.25)))) + (Kw / sin(Angle2Arc(sqrt(el * el + 2.25))));

	return Trop_delay;
}

void CUser::processUserData()
{
	//��ÿһ�����ǵ������������
	double* seudoRange = new double[satNum];

	double residualErr = 100; //����ѭ���������
	double delta = 1.0e-3;
	do {
		//ÿ�����ǽ��������㣬������α��������
		for (int i = 0; i < satNum; i++) {
			double err = processErr(i);
			seudoRange[i] = observeData.pseudo_range[i];// +err;
		}
		//������̣��õ����Ž⣬���زв�
		residualErr = CalculateUserPosition(seudoRange);
		XYZ2LambdaPhiH(observeData.upos.x, observeData.upos.y, observeData.upos.z, observeData.longitude, observeData.Latitude, observeData.altitude);//���µ�XYZд�뾭γ����
	} while (residualErr > delta);
	bGotStationPos = true;//��ɼ���󣬸�����Ϣ��������ȷ���û�ͼ����
	delete[] seudoRange;
	//�����γ������
	FILE* file=fopen("./position.txt", "w");
	fprintf(file,"x:%.10lf,y:%.10lf,z:%.10lf\n", observeData.longitude, observeData.Latitude, observeData.altitude);
	fclose(file);
}

double CUser::CalculateUserPosition(double* seudoDis)
{
	//����C++��������Eigen����д�û�λ�ü������
	double* approxRou = new double[satNum];
	VectorXd L(satNum);
	MatrixXd A(satNum, 4);

	for (int i = 0; i < satNum; ++i)
	{
		approxRou[i] = sqrt(pow(sate[i].curPosition.x - observeData.upos.x, 2) + pow(sate[i].curPosition.y - observeData.upos.y, 2) + pow(sate[i].curPosition.z - observeData.upos.z, 2));
		L(i) = seudoDis[i] - approxRou[i];
		A(i, 0) = -(sate[i].curPosition.x - observeData.upos.x) / approxRou[i];
		A(i, 1) = -(sate[i].curPosition.y - observeData.upos.y) / approxRou[i];
		A(i, 2) = -(sate[i].curPosition.z - observeData.upos.z) / approxRou[i];
		A(i, 3) = -1;
	}
	MatrixXd N = A.transpose() * A;
	MatrixXd K = A.transpose() * L;
	VectorXd v = N.inverse() * K;
	observeData.upos.x += v(0);
	observeData.upos.y += v(1);
	observeData.upos.z += v(2);
	observeData.t += v(3) / LIGHTSPEED;
	delete[] approxRou;
	double err = v(0) * v(0) + v(1) * v(1) + v(2) * v(2) + (v(3) / LIGHTSPEED) * (v(3) / LIGHTSPEED);
	return err;
}

void CUser::render()
{

	for (int i = 0; i < satNum; ++i)
	{
		float clr[3] = { rand(), rand(), rand() };
		sate[i].render(clr);
	}

	if (bGotStationPos)//�����ɼ��㣬������û�λ��
	{
		// render user position

		if (bDrawsigLines)
		{
			glPushAttrib(GL_ALL_ATTRIB_BITS);
			glEnable(GL_LINE_SMOOTH);

			float lwidth = 0.5f;
			glLineWidth(lwidth);
			glColor3f(1.0f, 1.0f, 0.0f);

			glPushMatrix();
			//glTranslatef(observeData.upos.y / SCALE, observeData.upos.z / SCALE, observeData.upos.x / SCALE);
			auxWireOctahedron(7);
			glPopMatrix();

			if (bDrawsigLines)
			{
				glPushMatrix();
				glBegin(GL_LINES);
				for (int i = 0; i < satNum; ++i)
				{
					glVertex3f(observeData.upos.y / SCALE, observeData.upos.z / SCALE, observeData.upos.x / SCALE);
					glVertex3f(sate[i].curPosition.y / SCALE, sate[i].curPosition.z / SCALE, sate[i].curPosition.x / SCALE);
				}
				glEnd();
				glPopMatrix();
			}
			glPopAttrib();
		}
	}


}
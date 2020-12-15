#include "CUser.h"


CUser::CUser()
{

	satNum = SAT_NUM; //本程序中我们只处理6颗卫星的数据
	sate = new CSatellite[satNum];

	//初始化伪距变量
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

	//初始化用户初始位置
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


	// 编写自己的代码
	//读取n颗卫星的星历数据
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
	//处理轨道计算

	// 修正后的角速度n
	double n0 = sqrt(GM) / (pow(e.rA, 3));
	double n = n0 + e.dn;

	// 平近点角M

	double delta_t = e.observt - e.toe;
	double M = e.M0 + n * delta_t;

	// 计算偏近点角E
	double E = M;
	double deltaE = 0.000001, dE = 1;
	while (abs(dE) > abs(deltaE))
	{
		double En = calculateE(e.e, M, E);
		dE = En - E;
		E = En;
	}

	// 计算真近点角f
	double f = atan2(sqrt(1 - e.e * e.e) * sin(E), (cos(E) - e.e));

	// 计算非精确u，半径rt，倾角it
	double raw_u = f + e.w;
	double rt = e.rA * e.rA * (1 - e.e * cos(E));
	double it = e.I0 + e.I_dot * delta_t;

	// 修正量
	double delta_u = e.Cuc * cos(2 * raw_u) + e.Cus * sin(2 * raw_u);
	double delta_r = e.Crc * cos(2 * raw_u) + e.Crs * sin(2 * raw_u);
	double delta_i = e.Cic * cos(2 * raw_u) + e.Cis * sin(2 * raw_u);

	// 精确r，u，i
	double r = rt + delta_r, u = raw_u + delta_u, i = it + delta_i;

	// 平面x，y
	double x = r * cos(u), y = r * sin(u);

	// 计算L
	double L = e.omega + e.omega_dot * delta_t - Omega_dote * e.observt;

	double X = x * cos(L) - y * cos(i) * sin(L);
	double Y = x * sin(L) + y * cos(i) * cos(L);
	double Z = y * sin(i);

	// 存储
	sate->orbit.I = i;
	sate->orbit.L = L;
	sate->curPosition.x = X;
	sate->curPosition.y = Y;
	sate->curPosition.z = Z;
	return;
}


//此函数一次完成所有工作，如想将数据读取与位置计算分开，可不调用该函数
void CUser::processSateData()
{
	readEphemeris();//一次性读取所有卫星数据

	for (int i = 0; i < satNum; i++)
	{
		CalculateSatPosition(&sate[i]);
	}

	datahasRead = true;

}


//误差处理主函数
double CUser::processErr(short sateID)
{

	//编程时可取消下面注释，定义常数LIGHTSPEED后，正常使用下面代码


	double dtrop = Tropospheric_Hopfield_err(sateID);//对流层误差

	double dion = Ionospheric_Klobuchar_err(sateID) * LIGHTSPEED;//电离层误差

	double dclk = Clockerr(sateID) * LIGHTSPEED;//时钟误差

	double drclk = Clockrelative_err(sateID) * LIGHTSPEED;//相对论误差


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
* 以下代码实现Klobchar法单频电离层改正模型和Hopefield法对流层改正模型
* 由于给出的参数并不明确，这里没有直接将想要的参数添加在CUser类中，只是以函数参数的形式写再来函数体中
* Klobuchar需要：
* 		测站三维坐标，不论是哪一种形式，目前默认给定的是三维直角坐标XYZ，在函数中是Ux，Uy，Uz；
* 		卫星三维坐标，这里将卫星对象传入了函数作为参数；
* 		世界时，这里采用观测时间除以3600再对24求模计算；
* 		星历给出的用来计算公式中参数alpha，beta，这里从行李结构体中读取；
* Hopefield需要：
* 		测站三维坐标，不论是哪一种形式，目前默认给定的是三维直角坐标XYZ，在函数中是Ux，Uy，Uz；
* 		卫星三维坐标，这里将卫星对象传入了函数作为参数；
* 		测站气压、水汽压、温度，这里都作为参数传入；
* 如何将计算结果加入整体代码还需要更清晰的参数传入；
*/

// 将角度转换为弧度
double Angle2Arc(double angle)
{
	return PI * angle / 180;
}
// 将弧度转换为角度
double Arc2Angle(double arc)
{
	return 180 * arc / PI;
}

// 将卫星坐标准换到测站坐标，利用测站和卫星之间坐标的差值和测站经纬度高程就可以通过矩阵乘法得到结果
// 函数接收测站和卫星的三维坐标差、测站经纬度高程、ENU坐标变量的引用
void CUser::LambdaPhiH2ENU(double dx, double dy, double dz, double lambda, double phi, double h, double& E, double& N, double& U)
{
	lambda = Angle2Arc(lambda);
	phi = Angle2Arc(phi);
	// 这里因为矩阵乘法可以用三个表达式计算出结果，就直接利用乘法解决了
	E = -sin(lambda) * dx + cos(lambda) * dy;
	N = -sin(phi) * cos(lambda) * dx - -sin(phi) * sin(lambda) * dy + cos(phi) * dz;
	U = cos(phi) * cos(lambda) * dx + cos(phi) * sin(lambda) * dy + sin(phi) * dz;
}

// 将XYZ坐标转换为经纬度高程坐标，接收XYZ坐标和经纬度高程的引用
// 使用老师在PPT中给出的计算方法，省去了迭代的麻烦
void CUser::XYZ2LambdaPhiH(double x, double y, double z, double& lambda, double& phi, double& h)
{
	// 必要参数计算
	double a = 6378137; // 地球椭球体半长轴
	double e = 8.1819190842622e-2; // 地球椭球体偏心率
	double b = sqrt(a * a * (1 - e * e)); // 短半轴
	double ep = sqrt((a * a - b * b) / (b * b)); // 第二偏心率
	double p = sqrt(x * x + y * y); // 参数

	lambda = atan2(y, x); // 计算经度

	double th = atan2((a * z) , (b * p));// 参数
	phi = atan2(z + ep * ep * b * pow(sin(th), 3), p - e * e * a * pow(cos(th), 3)); // 计算纬度

	double N = a / sqrt(1 - e * e * pow(sin(phi), 2)); // 参数
	h = p / cos(phi) - N; // 计算高程
	lambda = Arc2Angle(lambda);
	phi = Arc2Angle(phi);
}

double CUser::Ionospheric_Klobuchar_err(short sateID)
{
	double Iono_delay;
	CSatellite sate = this->sate[sateID];

	// 这是从星历里读取出来的，用作计算P和A
	double* alpha = (observeData.iono_A), * beta = (observeData.iono_B);

	// 测站和卫星坐标差
	double dx = sate.curPosition.x - observeData.upos.x, dy = sate.curPosition.y - observeData.upos.y, dz = sate.curPosition.z - observeData.upos.z;

	// 经纬度高程坐标
	double lambda = observeData.longitude, phi = observeData.Latitude, h = observeData.altitude;

	double UT = lambda * 240 + (int)observeData.t % 86400; // 计算当前标准世界时

	// 将经纬度高程坐标准换为ENU坐标。也就是将卫星坐标转换到测站坐标系下
	double E = 0, N = 0, U = 0;
	LambdaPhiH2ENU(dx, dy, dz, lambda, phi, h, E, N, U);

	// 必要参数，这里所有的角度都使用弧度制，因为math.h使用的都是弧度制，在需要的时候再进行转换
	double norm = sqrt(E * E + N * N + U * U); // 模
	double el = asin(U / norm);    // 卫星高度角
	double Azim = atan2(E / norm, N / norm);   // 穿刺点的星下位置方位角
	double EA = 445 / (Arc2Angle(el) + 20) - 4; // 卫星位置和测站位置相对地心夹角
	double phiN = phi + EA * cos(Azim);   // 穿刺点地心纬度
	double lambdaN = lambda + EA * sin(Azim) / cos(Angle2Arc(phi)); // 穿刺点地心经度
	double phiM = Angle2Arc(phiN + 10.07 * cos(Angle2Arc(lambdaN - 288.04))); // 穿刺点N的地磁纬度 

	double t = lambda * 240 + (int)observeData.t % 86400;; // 穿刺点N的地方时
	if (t > 86400) t -= 86400; // 计算结果如果大于24小时，则减去24h的秒数
	else if (t < 0) t += 86400; // 计算结果如果小于0，由应转为正值


	// 计算出最终公式中的两个参数，利用了星历中给出的alpha和beta
	double A = 0, P = 0;
	for (int i = 0; i < 4; ++i)
	{
		A += alpha[i] * pow(phiM, i);
		P += beta[i] * pow(phiM, i);
	}

	// 计算出穿刺点天顶方向的信号时延
	double Tg = 5e-9 + A * cos((2 * (double)PI / P) * (t - 50400));
	// 信号路径方向和天顶方向改正系数
	double secZ = 1 + (96 - Arc2Angle(el)) / 45;
	// 返回信号方向时延
	Iono_delay = Tg * secZ;

	return Iono_delay;
}

double CUser::Tropospheric_Hopfield_err(short sateID)
{
	double Trop_delay;
	CSatellite sate = this->sate[sateID];
	// 测站和卫星坐标差
	double dx = sate.curPosition.x - observeData.upos.x, dy = sate.curPosition.y - observeData.upos.y, dz = sate.curPosition.z - observeData.upos.z;

	// 经纬度高程坐标
	double lambda = observeData.longitude, phi = observeData.Latitude, h = observeData.altitude;

	// 将经纬度高程坐标准换为ENU坐标。也就是将卫星坐标转换到测站坐标系下
	double E = 0, N = 0, U = 0;
	LambdaPhiH2ENU(dx, dy, dz, lambda, phi, h, E, N, U);

	double norm = sqrt(E * E + N * N + U * U); // 模
	double el = Arc2Angle(asin(U / norm)); // 卫星高度角

	double hw = 11000; // 湿气部分的高度
	double hd = 40136 + 148.72 * (observeData.T - 273.16); // 干气部分的高度
	double Kd = 155.2e-7 * observeData.P * (hd - h) / observeData.T;  // 干气部分的延长
	double Kw = 155.2e-7 * 4810 * observeData.Pva * (hw - h) / observeData.T / observeData.T; // 湿气部分的延长

	// 返回在干湿方向上投影完成的加和，也就是对流层对信号的延迟
	Trop_delay = (Kd / sin(Angle2Arc(sqrt(el * el + 6.25)))) + (Kw / sin(Angle2Arc(sqrt(el * el + 2.25))));

	return Trop_delay;
}

void CUser::processUserData()
{
	//对每一颗卫星单独计算其误差
	double* seudoRange = new double[satNum];

	double residualErr = 100; //控制循环计算过程
	double delta = 1.0e-3;
	do {
		//每颗卫星进行误差计算，改正其伪距测量误差
		for (int i = 0; i < satNum; i++) {
			double err = processErr(i);
			seudoRange[i] = observeData.pseudo_range[i];// +err;
		}
		//求解误差方程，得到最优解，返回残差
		residualErr = CalculateUserPosition(seudoRange);
		XYZ2LambdaPhiH(observeData.upos.x, observeData.upos.y, observeData.upos.z, observeData.longitude, observeData.Latitude, observeData.altitude);//将新的XYZ写入经纬度中
	} while (residualErr > delta);
	bGotStationPos = true;//完成计算后，给出信息，方便正确调用绘图函数
	delete[] seudoRange;
	//输出经纬度坐标
	FILE* file=fopen("./position.txt", "w");
	fprintf(file,"x:%.10lf,y:%.10lf,z:%.10lf\n", observeData.longitude, observeData.Latitude, observeData.altitude);
	fclose(file);
}

double CUser::CalculateUserPosition(double* seudoDis)
{
	//基于C++矩阵计算库Eigen，编写用户位置计算代码
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

	if (bGotStationPos)//如果完成计算，则绘制用户位置
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
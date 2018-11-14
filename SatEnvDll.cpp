#include "SatEnvDll.h"

#define PI 3.141592653589793/*圆周率*/
#define RAD2DEG 57.295779513082321/*弧度转换度数的量*/
#define DEG2RAD 0.01745329251994329/*度数转换弧度的量*/

//slfmath	数学函数

/**
* [FunctionName] [矩阵复制]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

pMtx1 = pMtx2

* [output] [pMtx1] [输出矩阵]
* [input]  [pMtx2] [输入矩阵]
* [input]  [row] [行数]
* [input]  [col] [列数]
**/
int mtxCpy(double *pMtx1, double *pMtx2, int row, int col)
{
	int i, j;

	for (i = 0; i < row; i++)
		for (j = 0; j < col; j++)
			*(pMtx1 + i * col + j) = *(pMtx2 + i * col + j);

	return 1;
}

/**
* [FunctionName] [矩阵乘法]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

pMtx0 = pMtx1 × pMtx2

* [output] [pMtx0] [输出矩阵]
* [input]  [pMtx1] [输入矩阵1]
* [input]  [row1] [行数1]
* [input]  [col1] [列数1]
* [input]  [pMtx2] [输入矩阵2]
* [input]  [row2] [行数2]
* [input]  [col2] [列数2]
**/
int mtxMtp(double *pMtx0, double *pMtx1, int row1, int col1, double *pMtx2, int row2, int col2)
{
	int i, j, k;

	if ((row1 < 2) && (col1 < 2))
	{
		for (i = 0; i < row2; i++)
			for (j = 0; j < col2; j++)
				*(pMtx0 + i * col2 + j) = (*(pMtx1)) * (*(pMtx2 + i * col2 + j));
		return 1;
	}
	if ((row2 < 2) && (col2 < 2))
	{
		for (i = 0; i < row1; i++)
			for (j = 0; j < col1; j++)
			{
				*(pMtx0 + i * col1 + j) = (*(pMtx1 + i * col1 + j)) *(*(pMtx2));
			}
		return 1;
	}
	if (fabs(col1 - row2) < 1)
	{
		for (i = 0; i < row1; i++)
			for (j = 0; j < col2; j++)
			{
				*(pMtx0 + i * col2 + j) = 0;
				for (k = 0; k < row2; k++)
					*(pMtx0 + i * col2 + j) += (*(pMtx1 + i * col1 + k)) * (*(pMtx2 + k * col2 + j));
			}
		return 1;
	}

	return 0;
}

/**
* [FunctionName] [向量求模]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

返回值 = norm(pMtx0)

* [output] [返回值] [向量的模]
* [input]  [pMtx0] [输入向量]
* [input]  [cnt] [向量维数]
**/
double norm(double *pMtx0, int cnt)
{
	int i;
	double tmp;

	tmp = 0;
	for (i = 0; i < cnt; i++)
		tmp += (*(pMtx0 + i))*(*(pMtx0 + i));
	tmp = sqrt(tmp);

	return tmp;
}

/**
* [FunctionName] [矩阵加法]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

pMtx0 = pMtx1 + pMtx2

* [output] [pMtx0] [输出矩阵]
* [input]  [pMtx1] [输入矩阵1]
* [input]  [pMtx2] [输入矩阵2]
* [input]  [row] [行数]
* [input]  [col] [列数]
**/
void mtxAdd(double *pMtx0, double *pMtx1, double *pMtx2, int row, int col)
{
	int i, j;

	for (i = 0; i < row; i++)
		for (j = 0; j < col; j++)
		{
			*(pMtx0 + i * col + j) = (*(pMtx1 + i * col + j)) + (*(pMtx2 + i * col + j));

		}


	return;
}

/**
* [FunctionName] [三维向量叉乘]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

pMtx0 = pMtx1 × pMtx2

* [output] [pMtx0] [输出向量]
* [input]  [pMtx1] [输入向量1]
* [input]  [pMtx2] [输入向量2]
**/
void VecCross(double *pMtx0, double *pMtx1, double *pMtx2)
{
	double Mtx[3][3];

	Mtx[0][0] = 0; Mtx[0][1] = -(*(pMtx1 + 2)); Mtx[0][2] = (*(pMtx1 + 1));
	Mtx[1][0] = (*(pMtx1 + 2)); Mtx[1][1] = 0; Mtx[1][2] = -(*pMtx1);
	Mtx[2][0] = -(*(pMtx1 + 1)); Mtx[2][1] = (*pMtx1); Mtx[2][2] = 0;

	mtxMtp(pMtx0, Mtx[0], 3, 3, pMtx2, 3, 1);

	return;
}

/**
* [FunctionName] [矩阵转置]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

pMtx0 = pMtx1’

* [output] [pMtx0] [输出矩阵]
* [input]  [pMtx1] [输入矩阵]
* [input]  [row] [行数]
* [input]  [col] [列数]
**/
void mtxT(double *pMtx0, double *pMtx1, int row, int col)
{
	int i, j;

	for (i = 0; i < row; i++)
		for (j = 0; j < col; j++)
			*(pMtx0 + i * col + j) = *(pMtx1 + j * row + i);

	return;
}

/**
* [FunctionName] [高斯随机数生成]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

生成均值为0，方差为1的高斯随机数

* [output] [返回值] [高斯随机数]
* [input]  [param] [definition]
**/
double randn()
{
	int i;
	double result;

	result = 0;
	for (i = 0; i < 12; i++)
		result += fmod(rand(), 360) / 359.0;

	result -= 6;

	return result;
}

/**
* [FunctionName] [均匀随机数生成]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

生成0~1均匀随机数

* [output] [返回值] [均匀随机数]
* [input]  [param] [definition]
**/
double rand01()
{
	double result;

	result = fmod(rand(), 360) / 359.0;

	return result;
}

/**
* [FunctionName] [矩阵减法]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

pMtx0 = pMtx1 - pMtx2

* [output] [pMtx0] [输出矩阵]
* [input]  [pMtx1] [输入矩阵1]
* [input]  [pMtx2] [输入矩阵2]
* [input]  [row] [行数]
* [input]  [col] [列数]
**/
void mtxSub(double *pMtx0, double *pMtx1, double *pMtx2, int row, int col)
{
	int i, j;

	for (i = 0; i < row; i++)
		for (j = 0; j < col; j++)
			*(pMtx0 + i * col + j) = (*(pMtx1 + i * col + j)) - (*(pMtx2 + i * col + j));

	return;
}

/**
* [FunctionName] [三阶矩阵求逆]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

a = a^-1

* [output] [a] [输出矩阵]
* [input]  [a] [输出矩阵]
* [input]  [n] [阶数]
**/
int mtxInv(double *a, int n)
{
	int /**is,*js,*/i, j, k, l, u, v;
	double d, p;
	double js[3], is[3];
	/*is = (int *)malloc(n*sizeof(int));
	js = (int *)malloc(n*sizeof(int));*/
	for (k = 0; k <= n - 1; k++)
	{
		d = 0.0;
		for (i = k; i <= n - 1; i++)
			for (j = k; j <= n - 1; j++)
			{
				l = i * n + j; p = fabs(*(a + l));
				if (p > d) { d = p; is[k] = i; js[k] = j; }
			}
		if (fabs(d) < (1e-6))
		{
			/* free(is); free(js);;*/
			return(0);
		}
		if (is[k] != k)
			for (j = 0; j <= n - 1; j++)
			{
				u = k * n + j; v = is[k] * n + j;
				p = *(a + u); *(a + u) = *(a + v); *(a + v) = p;
			}
		if (js[k] != k)
			for (i = 0; i <= n - 1; i++)
			{
				u = i * n + k; v = i * n + js[k];
				p = *(a + u); *(a + u) = *(a + v); *(a + v) = p;
			}
		l = k * n + k;
		*(a + l) = 1.0 / (*(a + l));
		for (j = 0; j <= n - 1; j++)
			if (j != k)
			{
				u = k * n + j; *(a + u) = (*(a + u)) * (*(a + l));
			}
		for (i = 0; i <= n - 1; i++)
			if (i != k)
				for (j = 0; j <= n - 1; j++)
					if (j != k)
					{
						u = i * n + j;
						*(a + u) = (*(a + u)) - (*(a + i * n + k))*(*(a + k * n + j));
					}
		for (i = 0; i <= n - 1; i++)
			if (i != k)
			{
				u = i * n + k; *(a + u) = -(*(a + u))*(*(a + l));
			}
	}
	for (k = n - 1; k >= 0; k--)
	{
		if (js[k] != k)
			for (j = 0; j <= n - 1; j++)
			{
				u = k * n + j; v = js[k] * n + j;
				p = *(a + u); *(a + u) = *(a + v); *(a + v) = p;
			}
		if (is[k] != k)
			for (i = 0; i <= n - 1; i++)
			{
				u = i * n + k; v = i * n + is[k];
				p = *(a + u); *(a + u) = *(a + v); *(a + v) = p;
			}
	}
	/* free(is); free(js);*/
	return(1);
}

//attitudeLib	姿态函数库

/**
* [FunctionName] [姿态初始化]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

由欧拉角和欧拉角速度初始化姿态四元数

* [output] [attitude] [四元数 + 角速度rad/s]
* [input]  [orbit] [惯性系下 位置m 速度m/s]
* [input]  [angle] [三轴欧拉角rad]
* [input]  [angleRate] [三轴角速度rad/s]
**/
void attInit(double attitude[7], double orbit[6], double angle[3], double angleRate[3])
{
	double MtxJtoO[3][3], MtxOtoB[3][3], MtxJtoB[3][3];

	attitude[4] = angleRate[0];
	attitude[5] = angleRate[1];
	attitude[6] = angleRate[2];
	/* 设置卫星的角速率 */
	MtxJtoOGetG(MtxJtoO, orbit); /*利用轨道信息获得轨道系相对于惯性系的姿态矩阵*/
	MtxOtoBGet(MtxOtoB, angle); /*计算本体系相对于轨道系的姿态矩阵*/
	mtxMtp(MtxJtoB[0], MtxOtoB[0], 3, 3, MtxJtoO[0], 3, 3);
	/* 矩阵相乘，获得本体系相对惯性系的姿态矩阵 */
	QuatGetG(attitude, MtxJtoB); /*从姿态矩阵中获得姿态四元数*/

	return;
}

/**
* [FunctionName] [获取本体相对惯性系姿态四元数]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

从姿态转换矩阵获取姿态四元数

* [output] [attitude] [四元数 + 角速度rad/s]
* [input]  [MtxJtoB] [惯性系到本体系姿态转换矩阵]
**/
void QuatGetG(double attitude[7], double MtxJtoB[3][3])
{
	double tmp;
	attitude[3] = sqrt(1 + MtxJtoB[0][0] + MtxJtoB[1][1] + MtxJtoB[2][2]) / 2; /*q4*/
	attitude[2] = (MtxJtoB[0][1] - MtxJtoB[1][0]) / 4 / attitude[3]; /*q3*/
	attitude[1] = (MtxJtoB[2][0] - MtxJtoB[0][2]) / 4 / attitude[3]; /*q2*/
	attitude[0] = (MtxJtoB[1][2] - MtxJtoB[2][1]) / 4 / attitude[3]; /*q1*/
	tmp = 1.0 / norm(attitude, 4);
	mtxMtp(attitude, attitude, 4, 1, &tmp, 1, 1);
	return;
}

/**
* [FunctionName] [姿态递推]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [attitude] [四元数 + 角速度rad/s]
* [input]  [JM] [转动惯量矩阵]
* [input]  [JMI] [转动惯量矩阵的逆]
* [input]  [deltaMom] [三轴角动量]
* [input]  [Tq] [三轴控制力矩]
* [input]  [ha] [迭代步长]
**/
void attProp(double attitude[7], double JM[3][3], double JMI[3][3], double deltaMom[3], double Tq[3], double ha)
{
	int i;
	double A11[4][4], TransA[4][4], attTmp[7], tmp;
	double FVec[3], FV[3], attk[4][3];

	/*由四元数和姿态角的关系进行姿态运动方程递推*/
	mtxCpy(attTmp, attitude, 7, 1);

	A11[0][0] = 0;
	A11[0][1] = attTmp[6];
	A11[0][2] = -attTmp[5];
	A11[0][3] = attTmp[4];

	A11[1][0] = -attTmp[6];
	A11[1][1] = 0;
	A11[1][2] = attTmp[4];
	A11[1][3] = attTmp[5];

	A11[2][0] = attTmp[5];
	A11[2][1] = -attTmp[4];
	A11[2][2] = 0;
	A11[2][3] = attTmp[6];

	A11[3][0] = -attTmp[4];
	A11[3][1] = -attTmp[5];
	A11[3][2] = -attTmp[6];
	A11[3][3] = 0;
	/*dq=0.5All.*q(四元数和姿态角的关系）*/
	tmp = 0.5*ha;
	mtxMtp(TransA[0], A11[0], 4, 4, &tmp, 1, 1);/* TransA=0.05All */
	for (i = 0; i < 4; i++)
	{
		TransA[i][i] += 1;
	}
	mtxMtp(attTmp, TransA[0], 4, 4, attitude, 4, 1);
	tmp = 1.0 / norm(attTmp, 4);
	mtxMtp(attitude, attTmp, 4, 1, &tmp, 1, 1);

	/*-----------------------四元数传播（一阶RK算法）-----------------------*/
	/*-----------------------姿态动力学递推-----------------------*/
	//角速度微分方程：w' = I^-1*(T - w × (Iw + h) - h')
	//动量轮转速不变，h' == 0

	// RK1	
	AngRateDiff(FV, &attTmp[4], JM, JMI, deltaMom, Tq); //计算 w'
	mtxMtp(attk[0], FV, 3, 1, &ha, 1, 1);		//计算K1

												// RK2
	tmp = 0.5;
	mtxMtp(FVec, attk[0], 1, 3, &tmp, 1, 1);			//计算 0.5 * K1
	mtxAdd(&attTmp[4], &attitude[4], FVec, 1, 3);		//计算 w + 0.5 * K1
	AngRateDiff(FV, &attTmp[4], JM, JMI, deltaMom, Tq);	//计算 w'
	mtxMtp(attk[1], FV, 3, 1, &ha, 1, 1);		//计算K2

												// RK3
	tmp = 0.5;
	mtxMtp(FVec, attk[1], 1, 3, &tmp, 1, 1);			//计算 0.5 * K2
	mtxAdd(&attTmp[4], &attitude[4], FVec, 1, 3);		//计算 w + 0.5 * K2
	AngRateDiff(FV, &attTmp[4], JM, JMI, deltaMom, Tq);	//计算 w'
	mtxMtp(attk[2], FV, 3, 1, &ha, 1, 1);		//计算K3

	mtxAdd(&attTmp[4], &attitude[4], attk[2], 1, 3);	//计算 w + K3
	AngRateDiff(FV, &attTmp[4], JM, JMI, deltaMom, Tq);	//计算 w'
	mtxMtp(attk[3], FV, 3, 1, &ha, 1, 1);		//计算K4

	mtxAdd(FV, attk[0], attk[1], 1, 3);
	mtxAdd(FV, FV, attk[2], 1, 3);
	mtxAdd(FV, FV, attk[3], 1, 3);
	mtxAdd(FV, FV, attk[1], 1, 3);
	mtxAdd(FV, FV, attk[2], 1, 3);
	/*以上三行计算h(K1+2K2+2K3+K4)*/
	tmp = 1.0 / 6.0;
	mtxMtp(FVec, FV, 1, 3, &tmp, 1, 1);
	mtxAdd((attitude + 4), (attitude + 4), FVec, 3, 1);/*y(i+1)=y(i)+1/6h(k1+2k2+2k3+k4)*/

													   /*四元数归一化*/
	tmp = 1.0 / norm(attitude, 4);
	mtxMtp(attitude, attitude, 4, 1, &tmp, 1, 1);

	return;
}

/**
* [FunctionName] [四元数转欧拉角]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [Eulera] [欧拉角 [roll pitch yaw] ]
* [input]  [Quat] [四元数 [q1 q2 q3 q0] ]
**/
void Quat2Eulera(double Quat[4], double Eulera[3])
{
	Eulera[0] = atan2(2 * Quat[3] * Quat[0] + 2 * Quat[1] * Quat[2],
		Quat[3] * Quat[3] - Quat[0] * Quat[0] - Quat[1] * Quat[1] + Quat[2] * Quat[2]);
	Eulera[1] = asin(2 * Quat[3] * Quat[1] - 2 * Quat[0] * Quat[2]);
	Eulera[2] = atan2(2 * Quat[3] * Quat[2] + 2 * Quat[0] * Quat[1],
		Quat[3] * Quat[3] + Quat[0] * Quat[0] - Quat[1] * Quat[1] - Quat[2] * Quat[2]);
}

/**
* [FunctionName] [获取旋转矩阵]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [RotMtx] [旋转矩阵3×3]
* [input]  [Angle] [旋转角度rad]
* [input]  [Axis] [旋转轴 X/Y/Z ]
**/
void RotMtx(double RotMtx[3][3], double Angle, char Axis)
{
	if (Axis == 'X')
	{
		RotMtx[0][0] = 1; RotMtx[0][1] = 0; RotMtx[0][2] = 0;
		RotMtx[1][0] = 0; RotMtx[1][1] = cos(Angle); RotMtx[1][2] = sin(Angle);
		RotMtx[2][0] = 0; RotMtx[2][1] = -1.0*sin(Angle); RotMtx[2][2] = cos(Angle);
	}
	if (Axis == 'Y')
	{
		RotMtx[0][0] = cos(Angle); RotMtx[0][1] = 0; RotMtx[0][2] = -1.0*sin(Angle);
		RotMtx[1][0] = 0; RotMtx[1][1] = 1; RotMtx[1][2] = 0;
		RotMtx[2][0] = sin(Angle); RotMtx[2][1] = 0; RotMtx[2][2] = cos(Angle);
	}
	if (Axis == 'Z')
	{
		RotMtx[0][0] = cos(Angle); RotMtx[0][1] = sin(Angle); RotMtx[0][2] = 0;
		RotMtx[1][0] = -1.0*sin(Angle); RotMtx[1][1] = cos(Angle); RotMtx[1][2] = 0;
		RotMtx[2][0] = 0; RotMtx[2][1] = 0; RotMtx[2][2] = 1;
	}
}

/**
* [FunctionName] [惯性系下四元数转换为轨道系下四元数]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [QuatO] [轨道系下四元数]
* [input]  [QuatJ] [惯性系下四元数]
* [input]  [Orb] [轨道位置、速度信息]
**/
void QuatJ2O(double QuatO[4], double QuatJ[4], double Orb[6])
{
	int i;
	double MtxJtoO[3][3];
	double tmpM[4][4];
	double q[4];
	double tmp;

	/*获得轨道系相对于惯性系的姿态矩阵*/
	MtxJtoOGetG(MtxJtoO, Orb);
	q[3] = sqrt(1 + MtxJtoO[0][0] + MtxJtoO[1][1] + MtxJtoO[2][2]) / 2; /*q4 */
	q[2] = (MtxJtoO[0][1] - MtxJtoO[1][0]) / 4 / q[3];   /*q3 */
	q[1] = (MtxJtoO[2][0] - MtxJtoO[0][2]) / 4 / q[3];   /*q2 */
	q[0] = (MtxJtoO[1][2] - MtxJtoO[2][1]) / 4 / q[3];	 /*q1 */

	for (i = 0; i < 3; i++)
	{
		q[i] = -q[i];
	}

	tmpM[0][0] = q[3];
	tmpM[0][1] = -q[2];
	tmpM[0][2] = q[1];
	tmpM[0][3] = q[0];

	tmpM[1][0] = q[2];
	tmpM[1][1] = q[3];
	tmpM[1][2] = -q[0];
	tmpM[1][3] = q[1];

	tmpM[2][0] = -q[1];
	tmpM[2][1] = q[0];
	tmpM[2][2] = q[3];
	tmpM[2][3] = q[2];

	tmpM[3][0] = -q[0];
	tmpM[3][1] = -q[1];
	tmpM[3][2] = -q[2];
	tmpM[3][3] = q[3];
	/*获得相对于轨道系的姿态四元数*/
	mtxMtp(q, tmpM[0], 4, 4, QuatJ, 4, 1);
	tmp = 1.0 / norm(q, 4);
	mtxMtp(QuatO, q, 4, 1, &tmp, 1, 1);
}

/**
* [FunctionName] [计算当前状态角速度的微分]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [dw] [角速度微分]
* [input]  [w] [当前角速度]
* [input]  [J] [转动惯量矩阵]
* [input]  [JI] [转动惯量矩阵逆]
* [input]  [h] [动量轮角动量]
* [input]  [T] [控制力矩]
**/
void AngRateDiff(double dw[3], double w[3], double J[3][3], double JI[3][3], double h[3], double T[3])
{
	double tmp = -1.0;
	double Vec[3], Vec1[3];

	mtxMtp(Vec, (double*)J, 3, 3, w, 3, 1);//计算 Iw
	mtxAdd(Vec, Vec, h, 3, 1);				//计算 (Iw + h)
	VecCross(Vec1, w, Vec);			//计算 w × (Iw + h)
	tmp = -1;
	mtxMtp(Vec1, Vec1, 3, 1, &tmp, 1, 1);		//计算 - w × (Iw + h)
	mtxAdd(Vec, Vec1, T, 3, 1);				//计算 T - w × (Iw + h)
	mtxMtp(dw, (double*)JI, 3, 3, Vec, 3, 1);		//计算 w'
}

//calTransMtx	转换矩阵计算


/**
* [FunctionName] [根据轨道信息获得惯性系到轨道系的姿态矩阵]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [MtxJtoO] [惯性系到轨道系的姿态矩阵]
* [input]  [orbit] [轨道信息]
**/
void MtxJtoOGetG(double MtxJtoO[3][3], double orbit[6])
{
	double XVec[3], YVec[3], ZVec[3], VVec[3];
	double tmp, tmpM[3][3];

	tmp = -1.0 / norm(orbit, 3);           /*norm为向量求模函数*/
	mtxMtp(ZVec, orbit, 3, 1, &tmp, 1, 1); /*矩阵乘法函数*/
										   /* 获得轨道系Z在惯性系的单位矢量 */
	tmp = 1.0 / norm(&orbit[3], 3);
	mtxMtp(VVec, &orbit[3], 3, 1, &tmp, 1, 1);

	tmpM[0][0] = 0;
	tmpM[0][1] = -ZVec[2];
	tmpM[0][2] = ZVec[1];

	tmpM[1][0] = ZVec[2];
	tmpM[1][1] = 0;
	tmpM[1][2] = -ZVec[0];

	tmpM[2][0] = -ZVec[1];
	tmpM[2][1] = ZVec[0];
	tmpM[2][2] = 0;

	mtxMtp(YVec, tmpM[0], 3, 3, VVec, 3, 1);
	tmp = 1.0 / norm(YVec, 3);
	mtxMtp(YVec, YVec, 3, 1, &tmp, 1, 1);
	/* 获得轨道系Y在惯性系的单位矢量 */
	tmpM[0][0] = 0;
	tmpM[0][1] = -YVec[2];
	tmpM[0][2] = YVec[1];

	tmpM[1][0] = YVec[2];
	tmpM[1][1] = 0;
	tmpM[1][2] = -YVec[0];

	tmpM[2][0] = -YVec[1];
	tmpM[2][1] = YVec[0];
	tmpM[2][2] = 0;

	mtxMtp(XVec, tmpM[0], 3, 3, ZVec, 3, 1);
	/* 获得轨道系X在惯性系的单位矢量 */
	mtxCpy(MtxJtoO[0], XVec, 1, 3);
	mtxCpy(MtxJtoO[1], YVec, 1, 3);
	mtxCpy(MtxJtoO[2], ZVec, 1, 3);
	/* 利用XYZ矢量构成姿态矩阵 */
	return;
}

/**
* [FunctionName] [根据设定的321欧拉角，计算从轨道系到本体系的姿态矩阵]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [MtxOtoB] [解算后的轨道系到本体系的转换矩阵]
* [input]  [angle] [欧拉角rad]
**/
void MtxOtoBGet(double MtxOtoB[3][3], double angle[3])
{
	double ZMtx[3][3], YMtx[3][3], XMtx[3][3], XYMtx[3][3];
	double heading, pitch, roll, tmpsin, tmpcos;

	heading = angle[2];
	pitch = angle[1];
	roll = angle[0];

	/*heading = angle[0];
	pitch = angle[1];
	roll = angle[2];*/

	tmpsin = sin(heading);
	tmpcos = cos(heading);

	ZMtx[0][0] = tmpcos;
	ZMtx[0][1] = tmpsin;
	ZMtx[0][2] = 0;

	ZMtx[1][0] = -tmpsin;
	ZMtx[1][1] = tmpcos;
	ZMtx[1][2] = 0;

	ZMtx[2][0] = 0;
	ZMtx[2][1] = 0;
	ZMtx[2][2] = 1;
	/* 绕Z的旋转矩阵 */
	tmpsin = sin(pitch);
	tmpcos = cos(pitch);

	YMtx[0][0] = tmpcos;
	YMtx[0][1] = 0;
	YMtx[0][2] = -tmpsin;

	YMtx[1][0] = 0;
	YMtx[1][1] = 1;
	YMtx[1][2] = 0;

	YMtx[2][0] = tmpsin;
	YMtx[2][1] = 0;
	YMtx[2][2] = tmpcos;
	/* 绕Y的旋转矩阵 */
	tmpsin = sin(roll);
	tmpcos = cos(roll);

	XMtx[0][0] = 1;
	XMtx[0][1] = 0;
	XMtx[0][2] = 0;

	XMtx[1][0] = 0;
	XMtx[1][1] = tmpcos;
	XMtx[1][2] = tmpsin;

	XMtx[2][0] = 0;
	XMtx[2][1] = -tmpsin;
	XMtx[2][2] = tmpcos;
	/* 绕x的旋转矩阵 */
	mtxMtp((double *)XYMtx, (double *)XMtx, 3, 3, (double *)YMtx, 3, 3);
	mtxMtp((double *)MtxOtoB, (double *)XYMtx, 3, 3, (double *)ZMtx, 3, 3);
	/* 利用3个旋转构成姿态矩阵 */
	return;
}

/**
* [FunctionName] [实时计算从轨道系到本体系的姿态矩阵]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [MtxOtoB] [由轨道系到本体系的转换矩阵]
* [input]  [angle] [姿态角]
**/
void MtxOtoBGet_realtime(double MtxOtoB[3][3], double angle[3])
{
	double ZMtx[3][3], YMtx[3][3], XMtx[3][3], XYMtx[3][3];
	double heading, pitch, roll, tmpsin, tmpcos;
	/*初始姿态角角度到弧度的转换*/
	heading = angle[2]; /* 30 0 -25*/
	pitch = angle[1];
	roll = angle[0];

	tmpsin = sin(heading);
	tmpcos = cos(heading);

	ZMtx[0][0] = tmpcos;
	ZMtx[0][1] = tmpsin;
	ZMtx[0][2] = 0;

	ZMtx[1][0] = -tmpsin;
	ZMtx[1][1] = tmpcos;
	ZMtx[1][2] = 0;

	ZMtx[2][0] = 0;
	ZMtx[2][1] = 0;
	ZMtx[2][2] = 1;
	/* 绕Z的旋转矩阵 */
	tmpsin = sin(pitch);
	tmpcos = cos(pitch);

	YMtx[0][0] = tmpcos;
	YMtx[0][1] = 0;
	YMtx[0][2] = -tmpsin;

	YMtx[1][0] = 0;
	YMtx[1][1] = 1;
	YMtx[1][2] = 0;

	YMtx[2][0] = tmpsin;
	YMtx[2][1] = 0;
	YMtx[2][2] = tmpcos;
	/* 绕Y的旋转矩阵 */
	tmpsin = sin(roll);
	tmpcos = cos(roll);

	XMtx[0][0] = 1;
	XMtx[0][1] = 0;
	XMtx[0][2] = 0;

	XMtx[1][0] = 0;
	XMtx[1][1] = tmpcos;
	XMtx[1][2] = tmpsin;

	XMtx[2][0] = 0;
	XMtx[2][1] = -tmpsin;
	XMtx[2][2] = tmpcos;
	/* 绕x的旋转矩阵 */
	mtxMtp(XYMtx[0], XMtx[0], 3, 3, YMtx[0], 3, 3);
	mtxMtp(MtxOtoB[0], XYMtx[0], 3, 3, ZMtx[0], 3, 3);
	/* 利用3个旋转构成姿态矩阵 */
	return;
}

/**
* [FunctionName] [获得J2000到WGS84坐标的转换矩阵]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [JtoWGS] [J2000到WGS84坐标的转换矩阵]
* [input]  [JD] [儒略日]
* [input]  [date] [当前日期]
**/
void cordMtxJToWGSGetG(double JtoWGS[3][3], double JD, iTime_struct date)
{
	double Tu, worldTime; /*worldTime为世界时UT换成h单位时的值*/
	double tmpx, tmpy, tmpz, tmpsin, tmp, tmpcos;
	double i;
	/*计算儒略世纪=========(t - 2451545.0) / 36525.0*/
	Tu = JD - 2451545.0;
	Tu = Tu / 36525.0;

	tmpx = Tu;
	tmpy = Tu * Tu;
	tmpz = tmpy * Tu;

	/*求格林尼治世界时零时的恒星时
	ThetaG0 = 100.4606184 + 36000.77004 * Tu + 0.000387933 * Tu^2 - 2.583 * 10^-8 * Tu^3*/
	tmp = 100.4606184 + 36000.77004 * tmpx +
		0.000387933 * tmpy - 2.583 * pow(10.0, -8.0) * tmpz;
	if (tmp < 0)
	{
		for (i = 1.0; i <= 300.0; i++)
		{
			tmp = tmp + i * 360.0;
			if ((tmp >= 0) && (tmp <= 360))
			{
				break;
			}
		}
	}
	else if (tmp > 360)
	{
		for (i = 1.0; i <= 300.0; i++)
		{
			tmp = tmp - i * 360.0;
			if ((tmp >= 0) && (tmp <= 360))
			{
				break;
			}
		}
	}
	worldTime = date.hour + date.minute / 60.0 +
		date.second / 3600.0 + date.uSecond / 3600000000.0;
	/*计算其他任意世界时(UT)时刻的格林威治恒星时*/
	tmp += 360.98564724 * worldTime / 24.0; /*此为任意世界时UT的格林威治恒星时，单位为度*/
											/*再次限定范围为0-360度*/
	if (tmp < 0)
	{
		for (i = 1.0; i <= 300.0; i++)
		{
			tmp = tmp + i * 360.0;
			if ((tmp >= 0) && (tmp <= 360))
			{
				break;
			}
		}
	}
	else if (tmp > 360)
	{
		for (i = 1.0; i <= 300.0; i++)
		{
			tmp = tmp - i * 360.0;
			if ((tmp >= 0) && (tmp <= 360))
			{
				break;
			}
		}
	}
	tmpcos = cos(tmp);
	tmpsin = sin(tmp);

	JtoWGS[0][0] = tmpcos;
	JtoWGS[0][1] = tmpsin;
	JtoWGS[0][2] = 0;

	JtoWGS[1][0] = -tmpsin;
	JtoWGS[1][1] = tmpcos;
	JtoWGS[1][2] = 0;

	JtoWGS[2][0] = 0;
	JtoWGS[2][1] = 0;
	JtoWGS[2][2] = 1;
	return;
}

/**
* [FunctionName] [获得J2000到84坐标的转换矩阵]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [JtoWGS] [J2000到84坐标的转换矩阵]
* [input]  [JToFJ] [J2000到瞬时平赤道坐标系的转换矩阵]
* [input]  [tTime] [儒略日]
**/
void cordMtxJToWGSGet_Jin(double JtoWGS[3][3], double JToFJ[3][3], double tTime)
{
	//double PI = 3.1415926535897932384626433832795;
	double Tu, Td, DA, ZA, TA, KesiMean, OmegaMoon, LongiSun, DeltaKesi, DeltaFai, KesiReal;
	double tmpx, tmpy, tmpz, tmpsin, tmp, tmpcos, DATrans[3][3], ZATrans[3][3], TATrans[3][3];
	double TransKm[3][3], TransKr[3][3], TransDf[3][3];
	double tmpTrans[3][3], TransX[3][3], TransY[3][3], TransZ[3][3];

	Tu = tTime - 2451545;
	Td = Tu;
	Tu = Tu / 36525.0;

	tmpx = 180.0 * 3600 / PI;
	tmpy = Tu * Tu;
	tmpz = tmpy * Tu;
	DA = (2306.2181 * Tu + 0.30188 * tmpy + 0.017998 * tmpz) / tmpx;
	ZA = (2306.2181 * Tu + 1.09468 * tmpy + 0.018203 * tmpz) / tmpx;
	TA = (2004.3109 * Tu - 0.42665 * tmpy - 0.041833 * tmpz) / tmpx;
	tmpcos = cos(DA);
	tmpsin = sin(DA);
	DATrans[0][0] = tmpcos;
	DATrans[0][1] = -tmpsin;
	DATrans[0][2] = 0;
	DATrans[1][0] = tmpsin;
	DATrans[1][1] = tmpcos;
	DATrans[1][2] = 0;
	DATrans[2][0] = 0;
	DATrans[2][1] = 0;
	DATrans[2][2] = 1;
	tmpcos = cos(ZA);
	tmpsin = sin(ZA);
	ZATrans[0][0] = tmpcos;
	ZATrans[0][1] = -tmpsin;
	ZATrans[0][2] = 0;
	ZATrans[1][0] = tmpsin;
	ZATrans[1][1] = tmpcos;
	ZATrans[1][2] = 0;
	ZATrans[2][0] = 0;
	ZATrans[2][1] = 0;
	ZATrans[2][2] = 1;
	tmpcos = cos(TA);
	tmpsin = sin(TA);
	TATrans[0][0] = tmpcos;
	TATrans[0][1] = 0;
	TATrans[0][2] = -tmpsin;
	TATrans[1][0] = 0;
	TATrans[1][1] = 1;
	TATrans[1][2] = 0;
	TATrans[2][0] = tmpsin;
	TATrans[2][1] = 0;
	TATrans[2][2] = tmpcos;
	mtxMtp((double *)tmpTrans, (double *)ZATrans, 3, 3, (double *)TATrans, 3, 3);
	mtxMtp((double *)TransX, (double *)tmpTrans, 3, 3, (double *)DATrans, 3, 3);
	mtxCpy((double *)JToFJ, (double *)TransX, 3, 3);
	/* 从J2000到瞬时平赤道坐标系的转换关系 */
	tmp = 180.0 / PI;
	KesiMean = (84381.448 - 46.815 * Tu - 0.00059 * tmpy + 0.001813 * tmpz) / tmpx;
	OmegaMoon = (125.04 - 0.052954 * Td) / tmp;
	LongiSun = (280.47 + 0.98565 * Td) / tmp;
	DeltaKesi = (9.2025 * cos(OmegaMoon) + 0.5736 * cos(2 * LongiSun)) / tmpx;
	DeltaFai = (-17.1996 * sin(OmegaMoon) - 1.3187 * sin(2 * LongiSun)) / tmpx;
	KesiReal = KesiMean + DeltaKesi;
	tmpcos = cos(KesiMean);
	tmpsin = sin(KesiMean);
	TransKm[0][0] = 1;
	TransKm[0][1] = 0;
	TransKm[0][2] = 0;
	TransKm[1][0] = 0;
	TransKm[1][1] = tmpcos;
	TransKm[1][2] = tmpsin;
	TransKm[2][0] = 0;
	TransKm[2][1] = -tmpsin;
	TransKm[2][2] = tmpcos;
	tmpcos = cos(KesiReal);
	tmpsin = sin(KesiReal);
	TransKr[0][0] = 1;
	TransKr[0][1] = 0;
	TransKr[0][2] = 0;
	TransKr[1][0] = 0;
	TransKr[1][1] = tmpcos;
	TransKr[1][2] = -tmpsin;
	TransKr[2][0] = 0;
	TransKr[2][1] = tmpsin;
	TransKr[2][2] = tmpcos;
	tmpcos = cos(DeltaFai);
	tmpsin = sin(DeltaFai);
	TransDf[0][0] = tmpcos;
	TransDf[0][1] = -tmpsin;
	TransDf[0][2] = 0;
	TransDf[1][0] = tmpsin;
	TransDf[1][1] = tmpcos;
	TransDf[1][2] = 0;
	TransDf[2][0] = 0;
	TransDf[2][1] = 0;
	TransDf[2][2] = 1;
	mtxMtp((double *)tmpTrans, (double *)TransKr, 3, 3, (double *)TransDf, 3, 3);
	mtxMtp((double *)TransY, (double *)tmpTrans, 3, 3, (double *)TransKm, 3, 3);
	/* 从瞬时平赤道坐标系到真赤道坐标的转换关系 */
	tmp = 2 * PI / 24;
	tmpx = 18.697374558 * tmp;
	tmpy = 24.06570982441908 * tmp;
	tmpz = tmpx + tmpy * Td;
	tmp = tmpz + DeltaFai * cos(KesiReal);
	tmpcos = cos(tmp);
	tmpsin = sin(tmp);
	TransZ[0][0] = tmpcos;
	TransZ[0][1] = tmpsin;
	TransZ[0][2] = 0;
	TransZ[1][0] = -tmpsin;
	TransZ[1][1] = tmpcos;
	TransZ[1][2] = 0;
	TransZ[2][0] = 0;
	TransZ[2][1] = 0;
	TransZ[2][2] = 1;
	/* 从瞬时真赤道坐标系到准地固系的转换关系 */
	mtxMtp((double *)tmpTrans, (double *)TransZ, 3, 3, (double *)TransY, 3, 3);
	mtxMtp((double *)JtoWGS, (double *)tmpTrans, 3, 3, (double *)TransX, 3, 3);

	return;
}

/**
* [FunctionName] [获得地固系下的位置矢量]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [PosInWGS] [卫星在地固系下的位置矢量]
* [input]  [JToWGS] [惯性系到地固系的转移矩阵]
* [input]  [orbit] [位置速度信息]
**/
void posInWGSGetG(double PosInWGS[3], double JToWGS[3][3], double orbit[6])
{
	mtxMtp(PosInWGS, JToWGS[0], 3, 3, orbit, 3, 1);
	return;
}

/**
* [FunctionName] [获得大地系下的位置表示]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [GeoCord] [大地系下的位置表示]
* [input]  [PosInWGS] [地固系下的位置矢量]
**/
void geoInfoGetG(double GeoCord[3], double PosInWGS[3])
{
	double pi = 3.141592653589793;
	GeoCord[2] = norm(PosInWGS, 3);              /*半径*/
	GeoCord[1] = acos(PosInWGS[2] / GeoCord[2]); /*纬度*/
	if (fabs(PosInWGS[0]) < 1)
	{
		if (PosInWGS[1] > 0)
		{
			GeoCord[0] = pi / 2.0;
		}
		else
		{
			GeoCord[0] = pi * 3 / 2.0;
		}
	}
	else
	{
		GeoCord[0] = atan(PosInWGS[1] / PosInWGS[0]); /*经度*/

		if (PosInWGS[0] < 0)
		{
			GeoCord[0] = GeoCord[0] + pi;
		}
		if ((PosInWGS[0] > 0) & (PosInWGS[1] < 0))
		{
			GeoCord[0] = GeoCord[0] + 2 * pi;
		}
	}
	return;
}

/**
* [FunctionName] [计算从惯性系到本体系的转换矩阵]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [JToB] [惯性系到本体系的转移矩阵]
* [input]  [attitude] [卫星姿态信息]
**/
void calJToB(double JToB[3][3], double attitude[7])
{
	int i;
	double att[7];
	double JToBody[3][3];
	for (i = 0; i < 7; i++)
	{
		att[i] = attitude[i];
	}
	/*获得惯性系到本体系的转换矩阵A*/
	JToBody[0][0] = att[3] * att[3] + att[0] * att[0] -
		att[1] * att[1] - att[2] * att[2];
	JToBody[0][1] = 2 * (att[3] * att[2] + att[0] * att[1]);
	JToBody[0][2] = 2 * (-att[3] * att[1] + att[0] * att[2]);

	JToBody[1][0] = 2 * (-att[3] * att[2] + att[0] * att[1]);
	JToBody[1][1] = att[3] * att[3] - att[0] * att[0] +
		att[1] * att[1] - att[2] * att[2];
	JToBody[1][2] = 2 * (att[3] * att[0] + att[1] * att[2]);

	JToBody[2][0] = 2 * (att[3] * att[1] + att[0] * att[2]);
	JToBody[2][1] = 2 * (-att[3] * att[0] + att[1] * att[2]);
	JToBody[2][2] = att[3] * att[3] - att[0] * att[0] -
		att[1] * att[1] + att[2] * att[2];

	mtxCpy(JToB[0], JToBody[0], 3, 3);
	return;
}


//calTime	时间转换函数库

/**
* [FunctionName] [时间初始化]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [inSecond] [初始化需要的秒]
* [input]  [inUSecond] [初始化需要的微秒]
**/
double cal_JD_DateG(iTime_struct *newDate, iTime_struct refDate, int inSecond, int inUSecond)
{
	/*------tm结构体的定义------
	struct tm
	{
	int tm_year;	years since 1900
	int tm_mon;		months since January[0, 11]
	int tm_mday;	day of the month[1, 31]
	int tm_hour;	hours after midnight[0, 23]
	int tm_min;		minutes after the hour[0, 59]
	int tm_sec;		seconds after the minute[0, 59]
	};*/
	double JD;
	time_t t;
	__int64 ltime;
	struct tm date;
	//_time64(&ltime);
	/*先把秒数转换成日期（年月日时分秒）*/
	date.tm_year = refDate.year - 1900;
	date.tm_mon = refDate.month - 1;
	date.tm_mday = refDate.day;
	date.tm_hour = refDate.hour;
	date.tm_min = refDate.minute;
	date.tm_sec = refDate.second;
	date.tm_isdst = 0;
	ltime = mktime(&date);
	ltime += 28800;// mktime计算秒数少算8小时，原因可能为时区问题 wjq 2018年9月20日
	ltime += inSecond;

	gmtime_s(&date, &ltime);
	//(void) gmtime_r(&t, &date);

	/*把日期赋给仿真日期（年月日时分秒）*/
	(*newDate).year = date.tm_year + 1900;
	(*newDate).month = date.tm_mon + 1;
	(*newDate).day = date.tm_mday;
	(*newDate).hour = date.tm_hour;
	(*newDate).minute = date.tm_min;
	(*newDate).second = date.tm_sec;
	(*newDate).uSecond = inUSecond;
	/*转换成儒略日，存放仿真时间起点*/
	calJD(&JD, *newDate); /* 用以描述时间的参照起点 */
	return JD;
}

/**
* [FunctionName] [计算儒略日]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [返回值] [儒略日]
* [input]  [tTime] [计算后的儒略日]
* [input]  [iTime] [当前日期时间]
**/
double calJD(double *tTime, iTime_struct iTime)
{
	double tmpx, tmpy, tmpz, i, j, k, tmp;

	tmpx = ((iTime.month) + 9.0) / 12.0;
	modf(tmpx, &i);

	tmpy = (i + (iTime.year))*7.0 / 4.0;
	modf(tmpy, &j);

	tmpz = 275.0*(iTime.month) / 9.0;
	modf(tmpz, &k);

	*tTime = 367.0*(iTime.year) - j + k + (iTime.day) + 1721013.5; /*J0*/
	tmp = (iTime.hour) + (iTime.minute) / 60.0 + (iTime.second) / 3600.0 + (iTime.uSecond) / 3600000000.0; /*转换世界时为h单位*/
	*tTime = (*tTime) + tmp / 24.0; /*JD*/
	return *tTime;
}


//orbitLib	轨道函数库

/**
* [FunctionName] [轨道递推]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [newOrbit] [迭代后的位置速度]
* [input]  [lstOrbit] [迭代前的位置速度]
* [input]  [accel] [三轴加速度]
* [input]  [ha] [迭代步长]
**/
void orbProp_VR(double newOrbit[6], double lstOrbit[6], double ha, double accel[3])
{
	double GM = 3.986005e14; /*地球引力常数*/
	double J2 = 1.08264e-3;  /*J2摄动项*/
	double Re = 6378137;	 /*地球半径*/

	double a = 0;
	double tmp;				 /*J2摄动项引起的加速度*/
	double radius = 0;		 /*到地心的距离*/
	double orbTmp[6] = { 0 };  /*缓存变量*/
	double orbk[4][6] = { 0 }; /*四阶龙格库塔计算中间值*/
	double orb[6] = { 0 };	 /*积分结果差值*/

	mtxCpy(orbTmp, lstOrbit, 6, 1);
	mtxCpy(orb, lstOrbit, 1, 6);

	/* RK 第一步 */
	mtxMtp(orbk[0], &orbTmp[3], 3, 1, &ha, 1, 1);
	radius = norm(orbTmp, 3);
	orbk[0][3] = -GM * orbTmp[0] / pow(radius, 3) * (1 + 1.5 * J2 * pow(Re, 2) / pow(radius, 2) * (1 - 5 * pow(orbTmp[2], 2) / pow(radius, 2))) + accel[0];
	orbk[0][4] = -GM * orbTmp[1] / pow(radius, 3) * (1 + 1.5 * J2 * pow(Re, 2) / pow(radius, 2) * (1 - 5 * pow(orbTmp[2], 2) / pow(radius, 2))) + accel[1];
	orbk[0][5] = -GM * orbTmp[2] / pow(radius, 3) * (1 + 1.5 * J2 * pow(Re, 2) / pow(radius, 2) * (3 - 5 * pow(orbTmp[2], 2) / pow(radius, 2))) + accel[2];
	mtxMtp(&orbk[0][3], &orbk[0][3], 1, 3, &ha, 1, 1);

	/* RK 第二步 */
	a = 0.5;
	mtxMtp(orb, orbk[0], 1, 6, &a, 1, 1);
	mtxAdd(orbTmp, lstOrbit, orb, 1, 6);
	mtxMtp(orbk[1], &orbTmp[3], 3, 1, &ha, 1, 1);
	radius = norm(orbTmp, 3);
	orbk[1][3] = -GM * orbTmp[0] / pow(radius, 3) * (1 + 1.5 * J2 * pow(Re, 2) / pow(radius, 2) * (1 - 5 * pow(orbTmp[2], 2) / pow(radius, 2))) + accel[0];
	orbk[1][4] = -GM * orbTmp[1] / pow(radius, 3) * (1 + 1.5 * J2 * pow(Re, 2) / pow(radius, 2) * (1 - 5 * pow(orbTmp[2], 2) / pow(radius, 2))) + accel[1];
	orbk[1][5] = -GM * orbTmp[2] / pow(radius, 3) * (1 + 1.5 * J2 * pow(Re, 2) / pow(radius, 2) * (3 - 5 * pow(orbTmp[2], 2) / pow(radius, 2))) + accel[2];
	mtxMtp(&orbk[1][3], &orbk[1][3], 1, 3, &ha, 1, 1);

	/* RK 第三步 */
	mtxMtp(orb, orbk[1], 1, 6, &a, 1, 1);
	mtxAdd(orbTmp, lstOrbit, orb, 1, 6);
	mtxMtp(orbk[2], &orbTmp[3], 3, 1, &ha, 1, 1);
	radius = norm(orbTmp, 3);
	orbk[2][3] = -GM * orbTmp[0] / pow(radius, 3) * (1 + 1.5 * J2 * pow(Re, 2) / pow(radius, 2) * (1 - 5 * pow(orbTmp[2], 2) / pow(radius, 2))) + accel[0];
	orbk[2][4] = -GM * orbTmp[1] / pow(radius, 3) * (1 + 1.5 * J2 * pow(Re, 2) / pow(radius, 2) * (1 - 5 * pow(orbTmp[2], 2) / pow(radius, 2))) + accel[1];
	orbk[2][5] = -GM * orbTmp[2] / pow(radius, 3) * (1 + 1.5 * J2 * pow(Re, 2) / pow(radius, 2) * (3 - 5 * pow(orbTmp[2], 2) / pow(radius, 2))) + accel[2];
	mtxMtp(&orbk[2][3], &orbk[2][3], 1, 3, &ha, 1, 1);

	/* RK 第四步 */
	mtxAdd(orbTmp, lstOrbit, orbk[2], 1, 6);
	mtxMtp(orbk[3], &orbTmp[3], 3, 1, &ha, 1, 1);
	radius = norm(orbTmp, 3);
	orbk[3][3] = -GM * orbTmp[0] / pow(radius, 3) * (1 + 1.5 * J2 * pow(Re, 2) / pow(radius, 2) * (1 - 5 * pow(orbTmp[2], 2) / pow(radius, 2))) + accel[0];
	orbk[3][4] = -GM * orbTmp[1] / pow(radius, 3) * (1 + 1.5 * J2 * pow(Re, 2) / pow(radius, 2) * (1 - 5 * pow(orbTmp[2], 2) / pow(radius, 2))) + accel[1];
	orbk[3][5] = -GM * orbTmp[2] / pow(radius, 3) * (1 + 1.5 * J2 * pow(Re, 2) / pow(radius, 2) * (3 - 5 * pow(orbTmp[2], 2) / pow(radius, 2))) + accel[2];
	mtxMtp(&orbk[3][3], &orbk[3][3], 1, 3, &ha, 1, 1);

	mtxAdd(orb, orbk[0], orbk[1], 1, 6);
	mtxAdd(orb, orb, orbk[2], 1, 6);
	mtxAdd(orb, orb, orbk[3], 1, 6);
	mtxAdd(orb, orb, orbk[1], 1, 6);
	mtxAdd(orb, orb, orbk[2], 1, 6);
	a = 1.0 / 6.0;
	mtxMtp(orb, orb, 6, 1, &a, 1, 1);
	mtxAdd(orbTmp, lstOrbit, orb, 6, 1);
	mtxCpy(newOrbit, orbTmp, 6, 1);
	return;
}

/**
* [FunctionName] [轨道六根数计算位置和速度]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [orbInfo] [位置与速度]
* [input]  [kplInfo] [轨道六根数]
**/
void calOrbG(double orbInfo[6], double kplInfo[6])
{
	int i;
	double GM = 3.986005e14;
	/*t,*//*dt,*//* n,*/
	double xPVec[3], xQVec[3], xE;
	double rVec[3], rxPVec[3], rxQVec[3];
	double dxPVec[3], dxQVec[3], drVec[3];
	double drnorm;
	double tmpa, tmpb, tmpc, tmpd;
	/*xPVec[3]、xQVec[3]——近焦点坐标系中的位置矢量在地心赤道系中的投影p^、q^（单位矢量）*/

	xE = kplInfo[5];
	for (i = 0; i < 10; i++)
	{
		/*迭代计算偏近点角*/
		xE = kplInfo[5] + kplInfo[1] * sin(xE); /* 偏近点角Ei+1=M+e*sin(Ei)，Ei迭代初值为M*/
	}
	/****  以下见轨道力学p141，(4.44)  ***/
	xPVec[0] = cos(kplInfo[3])*cos(kplInfo[4])
		- sin(kplInfo[3])*sin(kplInfo[4])*cos(kplInfo[2]);
	xPVec[1] = sin(kplInfo[3])*cos(kplInfo[4])
		+ cos(kplInfo[3])*sin(kplInfo[4])*cos(kplInfo[2]);
	xPVec[2] = sin(kplInfo[4])*sin(kplInfo[2]);
	/*近焦点坐标系中的位置矢量x在地心赤道系中的三轴投影p^*/
	xQVec[0] = -cos(kplInfo[3])*sin(kplInfo[4])
		- sin(kplInfo[3])*cos(kplInfo[4])*cos(kplInfo[2]);
	xQVec[1] = -sin(kplInfo[3])*sin(kplInfo[4])
		+ cos(kplInfo[3])*cos(kplInfo[4])*cos(kplInfo[2]);
	xQVec[2] = cos(kplInfo[4])*sin(kplInfo[2]);
	/*近焦点坐标系中的位置矢量y在地心赤道系中的三轴投影q^*/
	tmpa = kplInfo[0] * (cos(xE) - kplInfo[1]); /*tmpa=rcos(f)=x^=a(cos(E)-e)，f——真近点角，p88*/
	tmpb = kplInfo[0] * sqrt(1 - kplInfo[1] * kplInfo[1])*sin(xE); /*tmpb=rsin(f)=y^=a*sqrt(1-e^2)sin(E),p88*/
	tmpc = sqrt(GM*kplInfo[0])*(-sin(xE)); /*近焦点坐标系中的速度分量dx,p61,88*/
	tmpd = sqrt(GM*kplInfo[0])*sqrt(1 - kplInfo[1] * kplInfo[1])*cos(xE); /*近焦点坐标系中的速度分量dy,p61,88*/

	mtxMtp(rxPVec, xPVec, 3, 1, &tmpa, 1, 1);
	mtxMtp(rxQVec, xQVec, 3, 1, &tmpb, 1, 1);
	mtxAdd(rVec, rxPVec, rxQVec, 3, 1);/***相对于近焦点系的位置矢量r=x^p^+y^q^（轨道力学p60、88）*/

	mtxMtp(dxPVec, xPVec, 3, 1, &tmpc, 1, 1);
	mtxMtp(dxQVec, xQVec, 3, 1, &tmpd, 1, 1);
	mtxAdd(drVec, dxPVec, dxQVec, 3, 1); /*相对于近焦点系的速度矢量dr=dx^p^+dy^q^*/

	drnorm = 1.0 / norm(rVec, 3);
	mtxMtp(drVec, drVec, 3, 1, &drnorm, 1, 1);

	orbInfo[0] = rVec[0];
	orbInfo[1] = rVec[1];
	orbInfo[2] = rVec[2];
	orbInfo[3] = drVec[0];
	orbInfo[4] = drVec[1];
	orbInfo[5] = drVec[2];

	return;
}

/**
* [FunctionName] [位置和速度计算轨道六根数]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [kpl] [解算后的六根数]
* [input]  [orbit] [位置和速度]
**/
void calKPLG(double kpl[6], double orbit[6])
{
	int i;
	double pi = 3.141592653589793; 	 /*圆周率*/
	double GM = 3.986005e14; 		   /*地球引力常数*/
	double vel[3] = { 0 }, pos[3] = { 0 }; /*速度，位置*/
	double r, v, vr; 				   /*位置，速度，切向速度*/
	double M; 						   /*平近点角*/
	double arrayH[3] = { 0 }, arrayN[3] = { 0 }, arrayE[3] = { 0 };
	double array1[3] = { 0 }, array2[3] = { 0 }, value1, value2;
	/*-------------------------------------------------*/
	double ikpl[6] = { 0 }; /*局部变量，ikpl轨道六根数,h,i,omega,e,w,theta*/
							/*
							--------注意：ikpl与全局变量的kpl定义不一样-------
							ikpl[0]:h,比角动量的模
							ikpl[1]:i,轨道倾角
							ikpl[2]:omega,升交点赤经
							ikpl[3]:e,偏心率
							ikpl[4]:w,近地点幅角
							ikpl[5]:theta,真近点角
							------------------------------------------------
							*/

							/*--------------------获取速度和位置------------------*/
							/*获取速度*/
	for (i = 0; i < 3; i++) {
		vel[i] = orbit[3 + i];
	}
	/*获取位置*/
	for (i = 0; i < 3; i++) {
		pos[i] = orbit[i];
	}
	/*------------由位置和速度计算出新的轨道六根数------------*/
	r = norm(pos, 3);
	v = norm(vel, 3);
	mtxMtp(&vr, pos, 1, 3, vel, 3, 1);
	vr = vr / r;
	VecCross(arrayH, pos, vel);
	/*=============比角动量=============*/
	ikpl[0] = norm(arrayH, 3);
	/*=============轨道倾角=============*/
	ikpl[1] = acos(arrayH[2] / ikpl[0]);
	/*=============升交点赤经=============*/
	array1[0] = 0.0;
	array1[1] = 0.0;
	array1[2] = 1.0;
	VecCross(arrayN, array1, arrayH);
	value1 = norm(arrayN, 3);
	if (arrayN[1] >= 0.0) {
		ikpl[2] = acos(arrayN[0] / value1);
	}
	else if (arrayN[1] < 0.0) {
		ikpl[2] = 2 * pi - acos(arrayN[0] / value1);
	}
	/*=============偏心率=============*/
	value1 = v * v - GM / r;
	for (i = 0; i < 3; i++) {
		array1[i] = pos[i] * value1;
	}
	for (i = 0; i < 3; i++) {
		array2[i] = vel[i] * vr * r;
	}
	mtxSub(arrayE, array1, array2, 1, 3);
	value1 = 1 / GM;
	for (i = 0; i < 3; i++) {
		arrayE[i] = arrayE[i] * value1;
	}
	ikpl[3] = norm(arrayE, 3);
	/*=============近地点幅角=============*/
	mtxMtp(&value1, arrayN, 1, 3, arrayE, 3, 1);
	value2 = norm(arrayN, 3) * ikpl[3];
	/*奇异点判断*/
	if (arrayE[2] >= 0.0) {
		ikpl[4] = acos(value1 / value2);
		if (value1 / value2 >= 1.0) {
			ikpl[4] = 0.0;
		}
		else if (value1 / value2 < -1.0) {
			ikpl[4] = pi;
		}
	}
	else if (arrayE[2] < 0.0) {
		ikpl[4] = 2 * pi - acos(value1 / value2);
		if (value1 / value2 >= 1.0) {
			ikpl[4] = 0.0;
		}
		else if (value1 / value2 < -1.0) {
			ikpl[4] = pi;
		}
	}
	/*=============真近点角=============*/
	/*mtxMtp(&value1, arrayE, 1, 3, pos, 3, 1);*/
	value1 = 0;
	for (i = 0; i < 3; i++) {
		value1 += arrayE[i] * pos[i];
	}
	value2 = ikpl[3] * r;
	if (vr >= 0.0) {
		if (value1 / value2 >= 1.0) {
			ikpl[5] = 0.0;
		}
		else if (value1 / value2 < -1.0) {
			ikpl[5] = pi;
		}
		else {
			ikpl[5] = acos(value1 / value2);
		}
	}
	else if (vr < 0.0) {
		if (value1 / value2 >= 1.0) {
			ikpl[5] = 0.0;
		}
		else if (value1 / value2 < -1.0) {
			ikpl[5] = pi;
		}
		else {
			ikpl[5] = 2 * pi - acos(value1 / value2);
		}
	}
	/*============平近点角============*/
	value1 = (1 - ikpl[3]) / (1 + ikpl[3]);
	value2 = 2 * atan(sqrt(value1) * tan(0.5 * ikpl[5]));
	M = value2 - ikpl[3] * sin(value2);
	ikpl[5] = M;
	/*-------------------更新轨道六根数-------------------*/
	value1 = norm(arrayH, 3);
	value1 = value1 * value1;
	kpl[0] = (value1 / GM) / (1 - ikpl[3] * ikpl[3]);
	kpl[1] = ikpl[3];
	kpl[2] = ikpl[1];
	kpl[3] = ikpl[2];
	kpl[4] = ikpl[4];
	kpl[5] = ikpl[5];
	if (kpl[5] >= 2 * pi)
	{
		kpl[5] -= 2 * pi;
	}
	else if (kpl[5] < 0)
	{
		kpl[5] += 2 * pi;
	}

	return;
}




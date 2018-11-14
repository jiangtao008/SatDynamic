#include "SatEnvDll.h"

#define PI 3.141592653589793/*Բ����*/
#define RAD2DEG 57.295779513082321/*����ת����������*/
#define DEG2RAD 0.01745329251994329/*����ת�����ȵ���*/

//slfmath	��ѧ����

/**
* [FunctionName] [������]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

pMtx1 = pMtx2

* [output] [pMtx1] [�������]
* [input]  [pMtx2] [�������]
* [input]  [row] [����]
* [input]  [col] [����]
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
* [FunctionName] [����˷�]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

pMtx0 = pMtx1 �� pMtx2

* [output] [pMtx0] [�������]
* [input]  [pMtx1] [�������1]
* [input]  [row1] [����1]
* [input]  [col1] [����1]
* [input]  [pMtx2] [�������2]
* [input]  [row2] [����2]
* [input]  [col2] [����2]
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
* [FunctionName] [������ģ]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

����ֵ = norm(pMtx0)

* [output] [����ֵ] [������ģ]
* [input]  [pMtx0] [��������]
* [input]  [cnt] [����ά��]
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
* [FunctionName] [����ӷ�]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

pMtx0 = pMtx1 + pMtx2

* [output] [pMtx0] [�������]
* [input]  [pMtx1] [�������1]
* [input]  [pMtx2] [�������2]
* [input]  [row] [����]
* [input]  [col] [����]
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
* [FunctionName] [��ά�������]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

pMtx0 = pMtx1 �� pMtx2

* [output] [pMtx0] [�������]
* [input]  [pMtx1] [��������1]
* [input]  [pMtx2] [��������2]
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
* [FunctionName] [����ת��]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

pMtx0 = pMtx1��

* [output] [pMtx0] [�������]
* [input]  [pMtx1] [�������]
* [input]  [row] [����]
* [input]  [col] [����]
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
* [FunctionName] [��˹���������]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

���ɾ�ֵΪ0������Ϊ1�ĸ�˹�����

* [output] [����ֵ] [��˹�����]
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
* [FunctionName] [�������������]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

����0~1���������

* [output] [����ֵ] [���������]
* [input]  [param] [definition]
**/
double rand01()
{
	double result;

	result = fmod(rand(), 360) / 359.0;

	return result;
}

/**
* [FunctionName] [�������]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

pMtx0 = pMtx1 - pMtx2

* [output] [pMtx0] [�������]
* [input]  [pMtx1] [�������1]
* [input]  [pMtx2] [�������2]
* [input]  [row] [����]
* [input]  [col] [����]
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
* [FunctionName] [���׾�������]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

a = a^-1

* [output] [a] [�������]
* [input]  [a] [�������]
* [input]  [n] [����]
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

//attitudeLib	��̬������

/**
* [FunctionName] [��̬��ʼ��]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

��ŷ���Ǻ�ŷ�����ٶȳ�ʼ����̬��Ԫ��

* [output] [attitude] [��Ԫ�� + ���ٶ�rad/s]
* [input]  [orbit] [����ϵ�� λ��m �ٶ�m/s]
* [input]  [angle] [����ŷ����rad]
* [input]  [angleRate] [������ٶ�rad/s]
**/
void attInit(double attitude[7], double orbit[6], double angle[3], double angleRate[3])
{
	double MtxJtoO[3][3], MtxOtoB[3][3], MtxJtoB[3][3];

	attitude[4] = angleRate[0];
	attitude[5] = angleRate[1];
	attitude[6] = angleRate[2];
	/* �������ǵĽ����� */
	MtxJtoOGetG(MtxJtoO, orbit); /*���ù����Ϣ��ù��ϵ����ڹ���ϵ����̬����*/
	MtxOtoBGet(MtxOtoB, angle); /*���㱾��ϵ����ڹ��ϵ����̬����*/
	mtxMtp(MtxJtoB[0], MtxOtoB[0], 3, 3, MtxJtoO[0], 3, 3);
	/* ������ˣ���ñ���ϵ��Թ���ϵ����̬���� */
	QuatGetG(attitude, MtxJtoB); /*����̬�����л����̬��Ԫ��*/

	return;
}

/**
* [FunctionName] [��ȡ������Թ���ϵ��̬��Ԫ��]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

����̬ת�������ȡ��̬��Ԫ��

* [output] [attitude] [��Ԫ�� + ���ٶ�rad/s]
* [input]  [MtxJtoB] [����ϵ������ϵ��̬ת������]
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
* [FunctionName] [��̬����]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [attitude] [��Ԫ�� + ���ٶ�rad/s]
* [input]  [JM] [ת����������]
* [input]  [JMI] [ת�������������]
* [input]  [deltaMom] [����Ƕ���]
* [input]  [Tq] [�����������]
* [input]  [ha] [��������]
**/
void attProp(double attitude[7], double JM[3][3], double JMI[3][3], double deltaMom[3], double Tq[3], double ha)
{
	int i;
	double A11[4][4], TransA[4][4], attTmp[7], tmp;
	double FVec[3], FV[3], attk[4][3];

	/*����Ԫ������̬�ǵĹ�ϵ������̬�˶����̵���*/
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
	/*dq=0.5All.*q(��Ԫ������̬�ǵĹ�ϵ��*/
	tmp = 0.5*ha;
	mtxMtp(TransA[0], A11[0], 4, 4, &tmp, 1, 1);/* TransA=0.05All */
	for (i = 0; i < 4; i++)
	{
		TransA[i][i] += 1;
	}
	mtxMtp(attTmp, TransA[0], 4, 4, attitude, 4, 1);
	tmp = 1.0 / norm(attTmp, 4);
	mtxMtp(attitude, attTmp, 4, 1, &tmp, 1, 1);

	/*-----------------------��Ԫ��������һ��RK�㷨��-----------------------*/
	/*-----------------------��̬����ѧ����-----------------------*/
	//���ٶ�΢�ַ��̣�w' = I^-1*(T - w �� (Iw + h) - h')
	//������ת�ٲ��䣬h' == 0

	// RK1	
	AngRateDiff(FV, &attTmp[4], JM, JMI, deltaMom, Tq); //���� w'
	mtxMtp(attk[0], FV, 3, 1, &ha, 1, 1);		//����K1

												// RK2
	tmp = 0.5;
	mtxMtp(FVec, attk[0], 1, 3, &tmp, 1, 1);			//���� 0.5 * K1
	mtxAdd(&attTmp[4], &attitude[4], FVec, 1, 3);		//���� w + 0.5 * K1
	AngRateDiff(FV, &attTmp[4], JM, JMI, deltaMom, Tq);	//���� w'
	mtxMtp(attk[1], FV, 3, 1, &ha, 1, 1);		//����K2

												// RK3
	tmp = 0.5;
	mtxMtp(FVec, attk[1], 1, 3, &tmp, 1, 1);			//���� 0.5 * K2
	mtxAdd(&attTmp[4], &attitude[4], FVec, 1, 3);		//���� w + 0.5 * K2
	AngRateDiff(FV, &attTmp[4], JM, JMI, deltaMom, Tq);	//���� w'
	mtxMtp(attk[2], FV, 3, 1, &ha, 1, 1);		//����K3

	mtxAdd(&attTmp[4], &attitude[4], attk[2], 1, 3);	//���� w + K3
	AngRateDiff(FV, &attTmp[4], JM, JMI, deltaMom, Tq);	//���� w'
	mtxMtp(attk[3], FV, 3, 1, &ha, 1, 1);		//����K4

	mtxAdd(FV, attk[0], attk[1], 1, 3);
	mtxAdd(FV, FV, attk[2], 1, 3);
	mtxAdd(FV, FV, attk[3], 1, 3);
	mtxAdd(FV, FV, attk[1], 1, 3);
	mtxAdd(FV, FV, attk[2], 1, 3);
	/*�������м���h(K1+2K2+2K3+K4)*/
	tmp = 1.0 / 6.0;
	mtxMtp(FVec, FV, 1, 3, &tmp, 1, 1);
	mtxAdd((attitude + 4), (attitude + 4), FVec, 3, 1);/*y(i+1)=y(i)+1/6h(k1+2k2+2k3+k4)*/

													   /*��Ԫ����һ��*/
	tmp = 1.0 / norm(attitude, 4);
	mtxMtp(attitude, attitude, 4, 1, &tmp, 1, 1);

	return;
}

/**
* [FunctionName] [��Ԫ��תŷ����]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [Eulera] [ŷ���� [roll pitch yaw] ]
* [input]  [Quat] [��Ԫ�� [q1 q2 q3 q0] ]
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
* [FunctionName] [��ȡ��ת����]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [RotMtx] [��ת����3��3]
* [input]  [Angle] [��ת�Ƕ�rad]
* [input]  [Axis] [��ת�� X/Y/Z ]
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
* [FunctionName] [����ϵ����Ԫ��ת��Ϊ���ϵ����Ԫ��]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [QuatO] [���ϵ����Ԫ��]
* [input]  [QuatJ] [����ϵ����Ԫ��]
* [input]  [Orb] [���λ�á��ٶ���Ϣ]
**/
void QuatJ2O(double QuatO[4], double QuatJ[4], double Orb[6])
{
	int i;
	double MtxJtoO[3][3];
	double tmpM[4][4];
	double q[4];
	double tmp;

	/*��ù��ϵ����ڹ���ϵ����̬����*/
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
	/*�������ڹ��ϵ����̬��Ԫ��*/
	mtxMtp(q, tmpM[0], 4, 4, QuatJ, 4, 1);
	tmp = 1.0 / norm(q, 4);
	mtxMtp(QuatO, q, 4, 1, &tmp, 1, 1);
}

/**
* [FunctionName] [���㵱ǰ״̬���ٶȵ�΢��]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [dw] [���ٶ�΢��]
* [input]  [w] [��ǰ���ٶ�]
* [input]  [J] [ת����������]
* [input]  [JI] [ת������������]
* [input]  [h] [�����ֽǶ���]
* [input]  [T] [��������]
**/
void AngRateDiff(double dw[3], double w[3], double J[3][3], double JI[3][3], double h[3], double T[3])
{
	double tmp = -1.0;
	double Vec[3], Vec1[3];

	mtxMtp(Vec, (double*)J, 3, 3, w, 3, 1);//���� Iw
	mtxAdd(Vec, Vec, h, 3, 1);				//���� (Iw + h)
	VecCross(Vec1, w, Vec);			//���� w �� (Iw + h)
	tmp = -1;
	mtxMtp(Vec1, Vec1, 3, 1, &tmp, 1, 1);		//���� - w �� (Iw + h)
	mtxAdd(Vec, Vec1, T, 3, 1);				//���� T - w �� (Iw + h)
	mtxMtp(dw, (double*)JI, 3, 3, Vec, 3, 1);		//���� w'
}

//calTransMtx	ת���������


/**
* [FunctionName] [���ݹ����Ϣ��ù���ϵ�����ϵ����̬����]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [MtxJtoO] [����ϵ�����ϵ����̬����]
* [input]  [orbit] [�����Ϣ]
**/
void MtxJtoOGetG(double MtxJtoO[3][3], double orbit[6])
{
	double XVec[3], YVec[3], ZVec[3], VVec[3];
	double tmp, tmpM[3][3];

	tmp = -1.0 / norm(orbit, 3);           /*normΪ������ģ����*/
	mtxMtp(ZVec, orbit, 3, 1, &tmp, 1, 1); /*����˷�����*/
										   /* ��ù��ϵZ�ڹ���ϵ�ĵ�λʸ�� */
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
	/* ��ù��ϵY�ڹ���ϵ�ĵ�λʸ�� */
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
	/* ��ù��ϵX�ڹ���ϵ�ĵ�λʸ�� */
	mtxCpy(MtxJtoO[0], XVec, 1, 3);
	mtxCpy(MtxJtoO[1], YVec, 1, 3);
	mtxCpy(MtxJtoO[2], ZVec, 1, 3);
	/* ����XYZʸ��������̬���� */
	return;
}

/**
* [FunctionName] [�����趨��321ŷ���ǣ�����ӹ��ϵ������ϵ����̬����]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [MtxOtoB] [�����Ĺ��ϵ������ϵ��ת������]
* [input]  [angle] [ŷ����rad]
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
	/* ��Z����ת���� */
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
	/* ��Y����ת���� */
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
	/* ��x����ת���� */
	mtxMtp((double *)XYMtx, (double *)XMtx, 3, 3, (double *)YMtx, 3, 3);
	mtxMtp((double *)MtxOtoB, (double *)XYMtx, 3, 3, (double *)ZMtx, 3, 3);
	/* ����3����ת������̬���� */
	return;
}

/**
* [FunctionName] [ʵʱ����ӹ��ϵ������ϵ����̬����]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [MtxOtoB] [�ɹ��ϵ������ϵ��ת������]
* [input]  [angle] [��̬��]
**/
void MtxOtoBGet_realtime(double MtxOtoB[3][3], double angle[3])
{
	double ZMtx[3][3], YMtx[3][3], XMtx[3][3], XYMtx[3][3];
	double heading, pitch, roll, tmpsin, tmpcos;
	/*��ʼ��̬�ǽǶȵ����ȵ�ת��*/
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
	/* ��Z����ת���� */
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
	/* ��Y����ת���� */
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
	/* ��x����ת���� */
	mtxMtp(XYMtx[0], XMtx[0], 3, 3, YMtx[0], 3, 3);
	mtxMtp(MtxOtoB[0], XYMtx[0], 3, 3, ZMtx[0], 3, 3);
	/* ����3����ת������̬���� */
	return;
}

/**
* [FunctionName] [���J2000��WGS84�����ת������]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [JtoWGS] [J2000��WGS84�����ת������]
* [input]  [JD] [������]
* [input]  [date] [��ǰ����]
**/
void cordMtxJToWGSGetG(double JtoWGS[3][3], double JD, iTime_struct date)
{
	double Tu, worldTime; /*worldTimeΪ����ʱUT����h��λʱ��ֵ*/
	double tmpx, tmpy, tmpz, tmpsin, tmp, tmpcos;
	double i;
	/*������������=========(t - 2451545.0) / 36525.0*/
	Tu = JD - 2451545.0;
	Tu = Tu / 36525.0;

	tmpx = Tu;
	tmpy = Tu * Tu;
	tmpz = tmpy * Tu;

	/*�������������ʱ��ʱ�ĺ���ʱ
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
	/*����������������ʱ(UT)ʱ�̵ĸ������κ���ʱ*/
	tmp += 360.98564724 * worldTime / 24.0; /*��Ϊ��������ʱUT�ĸ������κ���ʱ����λΪ��*/
											/*�ٴ��޶���ΧΪ0-360��*/
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
* [FunctionName] [���J2000��84�����ת������]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [JtoWGS] [J2000��84�����ת������]
* [input]  [JToFJ] [J2000��˲ʱƽ�������ϵ��ת������]
* [input]  [tTime] [������]
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
	/* ��J2000��˲ʱƽ�������ϵ��ת����ϵ */
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
	/* ��˲ʱƽ�������ϵ�����������ת����ϵ */
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
	/* ��˲ʱ��������ϵ��׼�ع�ϵ��ת����ϵ */
	mtxMtp((double *)tmpTrans, (double *)TransZ, 3, 3, (double *)TransY, 3, 3);
	mtxMtp((double *)JtoWGS, (double *)tmpTrans, 3, 3, (double *)TransX, 3, 3);

	return;
}

/**
* [FunctionName] [��õع�ϵ�µ�λ��ʸ��]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [PosInWGS] [�����ڵع�ϵ�µ�λ��ʸ��]
* [input]  [JToWGS] [����ϵ���ع�ϵ��ת�ƾ���]
* [input]  [orbit] [λ���ٶ���Ϣ]
**/
void posInWGSGetG(double PosInWGS[3], double JToWGS[3][3], double orbit[6])
{
	mtxMtp(PosInWGS, JToWGS[0], 3, 3, orbit, 3, 1);
	return;
}

/**
* [FunctionName] [��ô��ϵ�µ�λ�ñ�ʾ]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [GeoCord] [���ϵ�µ�λ�ñ�ʾ]
* [input]  [PosInWGS] [�ع�ϵ�µ�λ��ʸ��]
**/
void geoInfoGetG(double GeoCord[3], double PosInWGS[3])
{
	double pi = 3.141592653589793;
	GeoCord[2] = norm(PosInWGS, 3);              /*�뾶*/
	GeoCord[1] = acos(PosInWGS[2] / GeoCord[2]); /*γ��*/
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
		GeoCord[0] = atan(PosInWGS[1] / PosInWGS[0]); /*����*/

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
* [FunctionName] [����ӹ���ϵ������ϵ��ת������]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [JToB] [����ϵ������ϵ��ת�ƾ���]
* [input]  [attitude] [������̬��Ϣ]
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
	/*��ù���ϵ������ϵ��ת������A*/
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


//calTime	ʱ��ת��������

/**
* [FunctionName] [ʱ���ʼ��]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [inSecond] [��ʼ����Ҫ����]
* [input]  [inUSecond] [��ʼ����Ҫ��΢��]
**/
double cal_JD_DateG(iTime_struct *newDate, iTime_struct refDate, int inSecond, int inUSecond)
{
	/*------tm�ṹ��Ķ���------
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
	/*�Ȱ�����ת�������ڣ�������ʱ���룩*/
	date.tm_year = refDate.year - 1900;
	date.tm_mon = refDate.month - 1;
	date.tm_mday = refDate.day;
	date.tm_hour = refDate.hour;
	date.tm_min = refDate.minute;
	date.tm_sec = refDate.second;
	date.tm_isdst = 0;
	ltime = mktime(&date);
	ltime += 28800;// mktime������������8Сʱ��ԭ�����Ϊʱ������ wjq 2018��9��20��
	ltime += inSecond;

	gmtime_s(&date, &ltime);
	//(void) gmtime_r(&t, &date);

	/*�����ڸ����������ڣ�������ʱ���룩*/
	(*newDate).year = date.tm_year + 1900;
	(*newDate).month = date.tm_mon + 1;
	(*newDate).day = date.tm_mday;
	(*newDate).hour = date.tm_hour;
	(*newDate).minute = date.tm_min;
	(*newDate).second = date.tm_sec;
	(*newDate).uSecond = inUSecond;
	/*ת���������գ���ŷ���ʱ�����*/
	calJD(&JD, *newDate); /* ��������ʱ��Ĳ������ */
	return JD;
}

/**
* [FunctionName] [����������]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [����ֵ] [������]
* [input]  [tTime] [������������]
* [input]  [iTime] [��ǰ����ʱ��]
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
	tmp = (iTime.hour) + (iTime.minute) / 60.0 + (iTime.second) / 3600.0 + (iTime.uSecond) / 3600000000.0; /*ת������ʱΪh��λ*/
	*tTime = (*tTime) + tmp / 24.0; /*JD*/
	return *tTime;
}


//orbitLib	���������

/**
* [FunctionName] [�������]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [newOrbit] [�������λ���ٶ�]
* [input]  [lstOrbit] [����ǰ��λ���ٶ�]
* [input]  [accel] [������ٶ�]
* [input]  [ha] [��������]
**/
void orbProp_VR(double newOrbit[6], double lstOrbit[6], double ha, double accel[3])
{
	double GM = 3.986005e14; /*������������*/
	double J2 = 1.08264e-3;  /*J2�㶯��*/
	double Re = 6378137;	 /*����뾶*/

	double a = 0;
	double tmp;				 /*J2�㶯������ļ��ٶ�*/
	double radius = 0;		 /*�����ĵľ���*/
	double orbTmp[6] = { 0 };  /*�������*/
	double orbk[4][6] = { 0 }; /*�Ľ�������������м�ֵ*/
	double orb[6] = { 0 };	 /*���ֽ����ֵ*/

	mtxCpy(orbTmp, lstOrbit, 6, 1);
	mtxCpy(orb, lstOrbit, 1, 6);

	/* RK ��һ�� */
	mtxMtp(orbk[0], &orbTmp[3], 3, 1, &ha, 1, 1);
	radius = norm(orbTmp, 3);
	orbk[0][3] = -GM * orbTmp[0] / pow(radius, 3) * (1 + 1.5 * J2 * pow(Re, 2) / pow(radius, 2) * (1 - 5 * pow(orbTmp[2], 2) / pow(radius, 2))) + accel[0];
	orbk[0][4] = -GM * orbTmp[1] / pow(radius, 3) * (1 + 1.5 * J2 * pow(Re, 2) / pow(radius, 2) * (1 - 5 * pow(orbTmp[2], 2) / pow(radius, 2))) + accel[1];
	orbk[0][5] = -GM * orbTmp[2] / pow(radius, 3) * (1 + 1.5 * J2 * pow(Re, 2) / pow(radius, 2) * (3 - 5 * pow(orbTmp[2], 2) / pow(radius, 2))) + accel[2];
	mtxMtp(&orbk[0][3], &orbk[0][3], 1, 3, &ha, 1, 1);

	/* RK �ڶ��� */
	a = 0.5;
	mtxMtp(orb, orbk[0], 1, 6, &a, 1, 1);
	mtxAdd(orbTmp, lstOrbit, orb, 1, 6);
	mtxMtp(orbk[1], &orbTmp[3], 3, 1, &ha, 1, 1);
	radius = norm(orbTmp, 3);
	orbk[1][3] = -GM * orbTmp[0] / pow(radius, 3) * (1 + 1.5 * J2 * pow(Re, 2) / pow(radius, 2) * (1 - 5 * pow(orbTmp[2], 2) / pow(radius, 2))) + accel[0];
	orbk[1][4] = -GM * orbTmp[1] / pow(radius, 3) * (1 + 1.5 * J2 * pow(Re, 2) / pow(radius, 2) * (1 - 5 * pow(orbTmp[2], 2) / pow(radius, 2))) + accel[1];
	orbk[1][5] = -GM * orbTmp[2] / pow(radius, 3) * (1 + 1.5 * J2 * pow(Re, 2) / pow(radius, 2) * (3 - 5 * pow(orbTmp[2], 2) / pow(radius, 2))) + accel[2];
	mtxMtp(&orbk[1][3], &orbk[1][3], 1, 3, &ha, 1, 1);

	/* RK ������ */
	mtxMtp(orb, orbk[1], 1, 6, &a, 1, 1);
	mtxAdd(orbTmp, lstOrbit, orb, 1, 6);
	mtxMtp(orbk[2], &orbTmp[3], 3, 1, &ha, 1, 1);
	radius = norm(orbTmp, 3);
	orbk[2][3] = -GM * orbTmp[0] / pow(radius, 3) * (1 + 1.5 * J2 * pow(Re, 2) / pow(radius, 2) * (1 - 5 * pow(orbTmp[2], 2) / pow(radius, 2))) + accel[0];
	orbk[2][4] = -GM * orbTmp[1] / pow(radius, 3) * (1 + 1.5 * J2 * pow(Re, 2) / pow(radius, 2) * (1 - 5 * pow(orbTmp[2], 2) / pow(radius, 2))) + accel[1];
	orbk[2][5] = -GM * orbTmp[2] / pow(radius, 3) * (1 + 1.5 * J2 * pow(Re, 2) / pow(radius, 2) * (3 - 5 * pow(orbTmp[2], 2) / pow(radius, 2))) + accel[2];
	mtxMtp(&orbk[2][3], &orbk[2][3], 1, 3, &ha, 1, 1);

	/* RK ���Ĳ� */
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
* [FunctionName] [�������������λ�ú��ٶ�]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [orbInfo] [λ�����ٶ�]
* [input]  [kplInfo] [���������]
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
	/*xPVec[3]��xQVec[3]��������������ϵ�е�λ��ʸ���ڵ��ĳ��ϵ�е�ͶӰp^��q^����λʸ����*/

	xE = kplInfo[5];
	for (i = 0; i < 10; i++)
	{
		/*��������ƫ�����*/
		xE = kplInfo[5] + kplInfo[1] * sin(xE); /* ƫ�����Ei+1=M+e*sin(Ei)��Ei������ֵΪM*/
	}
	/****  ���¼������ѧp141��(4.44)  ***/
	xPVec[0] = cos(kplInfo[3])*cos(kplInfo[4])
		- sin(kplInfo[3])*sin(kplInfo[4])*cos(kplInfo[2]);
	xPVec[1] = sin(kplInfo[3])*cos(kplInfo[4])
		+ cos(kplInfo[3])*sin(kplInfo[4])*cos(kplInfo[2]);
	xPVec[2] = sin(kplInfo[4])*sin(kplInfo[2]);
	/*����������ϵ�е�λ��ʸ��x�ڵ��ĳ��ϵ�е�����ͶӰp^*/
	xQVec[0] = -cos(kplInfo[3])*sin(kplInfo[4])
		- sin(kplInfo[3])*cos(kplInfo[4])*cos(kplInfo[2]);
	xQVec[1] = -sin(kplInfo[3])*sin(kplInfo[4])
		+ cos(kplInfo[3])*cos(kplInfo[4])*cos(kplInfo[2]);
	xQVec[2] = cos(kplInfo[4])*sin(kplInfo[2]);
	/*����������ϵ�е�λ��ʸ��y�ڵ��ĳ��ϵ�е�����ͶӰq^*/
	tmpa = kplInfo[0] * (cos(xE) - kplInfo[1]); /*tmpa=rcos(f)=x^=a(cos(E)-e)��f���������ǣ�p88*/
	tmpb = kplInfo[0] * sqrt(1 - kplInfo[1] * kplInfo[1])*sin(xE); /*tmpb=rsin(f)=y^=a*sqrt(1-e^2)sin(E),p88*/
	tmpc = sqrt(GM*kplInfo[0])*(-sin(xE)); /*����������ϵ�е��ٶȷ���dx,p61,88*/
	tmpd = sqrt(GM*kplInfo[0])*sqrt(1 - kplInfo[1] * kplInfo[1])*cos(xE); /*����������ϵ�е��ٶȷ���dy,p61,88*/

	mtxMtp(rxPVec, xPVec, 3, 1, &tmpa, 1, 1);
	mtxMtp(rxQVec, xQVec, 3, 1, &tmpb, 1, 1);
	mtxAdd(rVec, rxPVec, rxQVec, 3, 1);/***����ڽ�����ϵ��λ��ʸ��r=x^p^+y^q^�������ѧp60��88��*/

	mtxMtp(dxPVec, xPVec, 3, 1, &tmpc, 1, 1);
	mtxMtp(dxQVec, xQVec, 3, 1, &tmpd, 1, 1);
	mtxAdd(drVec, dxPVec, dxQVec, 3, 1); /*����ڽ�����ϵ���ٶ�ʸ��dr=dx^p^+dy^q^*/

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
* [FunctionName] [λ�ú��ٶȼ�����������]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [kpl] [������������]
* [input]  [orbit] [λ�ú��ٶ�]
**/
void calKPLG(double kpl[6], double orbit[6])
{
	int i;
	double pi = 3.141592653589793; 	 /*Բ����*/
	double GM = 3.986005e14; 		   /*������������*/
	double vel[3] = { 0 }, pos[3] = { 0 }; /*�ٶȣ�λ��*/
	double r, v, vr; 				   /*λ�ã��ٶȣ������ٶ�*/
	double M; 						   /*ƽ�����*/
	double arrayH[3] = { 0 }, arrayN[3] = { 0 }, arrayE[3] = { 0 };
	double array1[3] = { 0 }, array2[3] = { 0 }, value1, value2;
	/*-------------------------------------------------*/
	double ikpl[6] = { 0 }; /*�ֲ�������ikpl���������,h,i,omega,e,w,theta*/
							/*
							--------ע�⣺ikpl��ȫ�ֱ�����kpl���岻һ��-------
							ikpl[0]:h,�ȽǶ�����ģ
							ikpl[1]:i,������
							ikpl[2]:omega,������ྭ
							ikpl[3]:e,ƫ����
							ikpl[4]:w,���ص����
							ikpl[5]:theta,������
							------------------------------------------------
							*/

							/*--------------------��ȡ�ٶȺ�λ��------------------*/
							/*��ȡ�ٶ�*/
	for (i = 0; i < 3; i++) {
		vel[i] = orbit[3 + i];
	}
	/*��ȡλ��*/
	for (i = 0; i < 3; i++) {
		pos[i] = orbit[i];
	}
	/*------------��λ�ú��ٶȼ�����µĹ��������------------*/
	r = norm(pos, 3);
	v = norm(vel, 3);
	mtxMtp(&vr, pos, 1, 3, vel, 3, 1);
	vr = vr / r;
	VecCross(arrayH, pos, vel);
	/*=============�ȽǶ���=============*/
	ikpl[0] = norm(arrayH, 3);
	/*=============������=============*/
	ikpl[1] = acos(arrayH[2] / ikpl[0]);
	/*=============������ྭ=============*/
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
	/*=============ƫ����=============*/
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
	/*=============���ص����=============*/
	mtxMtp(&value1, arrayN, 1, 3, arrayE, 3, 1);
	value2 = norm(arrayN, 3) * ikpl[3];
	/*������ж�*/
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
	/*=============������=============*/
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
	/*============ƽ�����============*/
	value1 = (1 - ikpl[3]) / (1 + ikpl[3]);
	value2 = 2 * atan(sqrt(value1) * tan(0.5 * ikpl[5]));
	M = value2 - ikpl[3] * sin(value2);
	ikpl[5] = M;
	/*-------------------���¹��������-------------------*/
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




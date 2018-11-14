#pragma once
#include <math.h>
#include <stdlib.h>
#include <time.h>

/**
* [FunctionName] [ʱ��ṹ��]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]
**/
typedef struct
{
	int year;
	int month;
	int day;
	int hour;
	int minute;
	int second;
	int uSecond;
} iTime_struct;

//slfmath:������ѧ���㺯��
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
int mtxCpy(double *pMtx1, double *pMtx2, int row, int col);

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
int mtxMtp(double *pMtx0, double *pMtx1, int row1, int col1, double *pMtx2, int row2, int col2);

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
double norm(double *pMtx0, int cnt);

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
void mtxAdd(double *pMtx0, double *pMtx1, double *pMtx2, int row, int col);

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
void VecCross(double *pMtx0, double *pMtx1, double *pMtx2);

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
void mtxT(double *pMtx0, double *pMtx1, int row, int col);

/**
* [FunctionName] [��˹���������]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

���ɾ�ֵΪ0������Ϊ1�ĸ�˹�����

* [output] [����ֵ] [��˹�����]
* [input]  [param] [definition]
**/
double randn();

/**
* [FunctionName] [�������������]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

����0~1���������

* [output] [����ֵ] [���������]
* [input]  [param] [definition]
**/
double rand01();

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
void mtxSub(double *pMtx0, double *pMtx1, double *pMtx2, int row, int col);

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
int mtxInv(double *a, int n);


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
void attInit(double attitude[7], double orbit[6], double angle[3], double angleRate[3]);

/**
* [FunctionName] [��ȡ������Թ���ϵ��̬��Ԫ��]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

����̬ת�������ȡ��̬��Ԫ��

* [output] [attitude] [��Ԫ�� + ���ٶ�rad/s]
* [input]  [MtxJtoB] [����ϵ������ϵ��̬ת������]
**/
void QuatGetG(double attitude[7], double MtxJtoB[3][3]);

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
void attProp(double attitude[7], double JM[3][3], double JMI[3][3], double deltaMom[3], double Tq[3], double ha);

/**
* [FunctionName] [��Ԫ��תŷ����]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [Eulera] [ŷ���� [roll pitch yaw] ]
* [input]  [Quat] [��Ԫ�� [q1 q2 q3 q0] ]
**/
void Quat2Eulera(double Quat[4], double Eulera[3]);

/**
* [FunctionName] [��ȡ��ת����]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [RotMtx] [��ת����3��3]
* [input]  [Angle] [��ת�Ƕ�rad]
* [input]  [Axis] [��ת�� X/Y/Z ]
**/
void RotMtx(double RotMtx[3][3], double Angle, char Axis);

/**
* [FunctionName] [����ϵ����Ԫ��ת��Ϊ���ϵ����Ԫ��]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [QuatO] [���ϵ����Ԫ��]
* [input]  [QuatJ] [����ϵ����Ԫ��]
* [input]  [Orb] [���λ�á��ٶ���Ϣ]
**/
void QuatJ2O(double QuatO[4], double QuatJ[4], double Orb[6]);

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
void AngRateDiff(double dw[3], double w[3], double J[3][3], double JI[3][3], double h[3], double T[3]);


//calTransMtx	ת���������

/**
* [FunctionName] [���ݹ����Ϣ��ù���ϵ�����ϵ����̬����]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [MtxJtoO] [����ϵ�����ϵ����̬����]
* [input]  [orbit] [�����Ϣ]
**/
void MtxJtoOGetG(double MtxJtoO[3][3], double orbit[6]);

/**
* [FunctionName] [�����趨��321ŷ���ǣ�����ӹ��ϵ������ϵ����̬����]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [MtxOtoB] [�����Ĺ��ϵ������ϵ��ת������]
* [input]  [angle] [ŷ����rad]
**/
void MtxOtoBGet(double MtxOtoB[3][3], double angle[3]);

/**
* [FunctionName] [ʵʱ����ӹ��ϵ������ϵ����̬����]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [MtxOtoB] [�ɹ��ϵ������ϵ��ת������]
* [input]  [angle] [��̬��]
**/
void MtxOtoBGet_realtime(double MtxOtoB[3][3], double angle[3]);

/**
* [FunctionName] [���J2000��WGS84�����ת������]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [JtoWGS] [J2000��WGS84�����ת������]
* [input]  [JD] [������]
* [input]  [date] [��ǰ����]
**/
void cordMtxJToWGSGetG(double JtoWGS[3][3], double JD, iTime_struct date);

/**
* [FunctionName] [���J2000��84�����ת������]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [JtoWGS] [J2000��84�����ת������]
* [input]  [JToFJ] [J2000��˲ʱƽ�������ϵ��ת������]
* [input]  [tTime] [������]
**/
void cordMtxJToWGSGet_Jin(double JtoWGS[3][3], double JToFJ[3][3], double tTime);

/**
* [FunctionName] [��õع�ϵ�µ�λ��ʸ��]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [PosInWGS] [�����ڵع�ϵ�µ�λ��ʸ��]
* [input]  [JToWGS] [����ϵ���ع�ϵ��ת�ƾ���]
* [input]  [orbit] [λ���ٶ���Ϣ]
**/
void posInWGSGetG(double PosInWGS[3], double JToWGS[3][3], double orbit[6]);

/**
* [FunctionName] [��ô��ϵ�µ�λ�ñ�ʾ]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [GeoCord] [���ϵ�µ�λ�ñ�ʾ]
* [input]  [PosInWGS] [�ع�ϵ�µ�λ��ʸ��]
**/
void geoInfoGetG(double GeoCord[3], double PosInWGS[3]);

/**
* [FunctionName] [����ӹ���ϵ������ϵ��ת������]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [JToB] [����ϵ������ϵ��ת�ƾ���]
* [input]  [attitude] [������̬��Ϣ]
**/
void calJToB(double JToB[3][3], double attitude[7]);

//calTime	ʱ��ת��������

/**
* [FunctionName] [ʱ���ʼ��]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [inSecond] [��ʼ����Ҫ����]
* [input]  [inUSecond] [��ʼ����Ҫ��΢��]
**/
double cal_JD_DateG(iTime_struct *newDate, iTime_struct refDate, int inSecond, int inUSecond);

/**
* [FunctionName] [����������]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [����ֵ] [������]
* [input]  [tTime] [������������]
* [input]  [iTime] [��ǰ����ʱ��]
**/
double calJD(double *tTime, iTime_struct iTime);


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
void calOrbG(double orbInfo[6], double kplInfo[6]);

/**
* [FunctionName] [�������������λ�ú��ٶ�]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [orbInfo] [λ�����ٶ�]
* [input]  [kplInfo] [���������]
**/
void orbProp_VR(double newOrbit[6], double lstOrbit[6], double ha, double accel[3]);

/**
* [FunctionName] [λ�ú��ٶȼ�����������]
* [Author] [�����]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018��11��14��]

* [output] [kpl] [������������]
* [input]  [orbit] [λ�ú��ٶ�]
**/
void calKPLG(double kpl[6], double orbit[6]);



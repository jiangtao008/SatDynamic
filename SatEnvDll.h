#pragma once
#include <math.h>
#include <stdlib.h>
#include <time.h>

/**
* [FunctionName] [时间结构体]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]
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

//slfmath:矩阵数学计算函数
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
int mtxCpy(double *pMtx1, double *pMtx2, int row, int col);

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
int mtxMtp(double *pMtx0, double *pMtx1, int row1, int col1, double *pMtx2, int row2, int col2);

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
double norm(double *pMtx0, int cnt);

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
void mtxAdd(double *pMtx0, double *pMtx1, double *pMtx2, int row, int col);

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
void VecCross(double *pMtx0, double *pMtx1, double *pMtx2);

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
void mtxT(double *pMtx0, double *pMtx1, int row, int col);

/**
* [FunctionName] [高斯随机数生成]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

生成均值为0，方差为1的高斯随机数

* [output] [返回值] [高斯随机数]
* [input]  [param] [definition]
**/
double randn();

/**
* [FunctionName] [均匀随机数生成]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

生成0~1均匀随机数

* [output] [返回值] [均匀随机数]
* [input]  [param] [definition]
**/
double rand01();

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
void mtxSub(double *pMtx0, double *pMtx1, double *pMtx2, int row, int col);

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
int mtxInv(double *a, int n);


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
void attInit(double attitude[7], double orbit[6], double angle[3], double angleRate[3]);

/**
* [FunctionName] [获取本体相对惯性系姿态四元数]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

从姿态转换矩阵获取姿态四元数

* [output] [attitude] [四元数 + 角速度rad/s]
* [input]  [MtxJtoB] [惯性系到本体系姿态转换矩阵]
**/
void QuatGetG(double attitude[7], double MtxJtoB[3][3]);

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
void attProp(double attitude[7], double JM[3][3], double JMI[3][3], double deltaMom[3], double Tq[3], double ha);

/**
* [FunctionName] [四元数转欧拉角]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [Eulera] [欧拉角 [roll pitch yaw] ]
* [input]  [Quat] [四元数 [q1 q2 q3 q0] ]
**/
void Quat2Eulera(double Quat[4], double Eulera[3]);

/**
* [FunctionName] [获取旋转矩阵]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [RotMtx] [旋转矩阵3×3]
* [input]  [Angle] [旋转角度rad]
* [input]  [Axis] [旋转轴 X/Y/Z ]
**/
void RotMtx(double RotMtx[3][3], double Angle, char Axis);

/**
* [FunctionName] [惯性系下四元数转换为轨道系下四元数]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [QuatO] [轨道系下四元数]
* [input]  [QuatJ] [惯性系下四元数]
* [input]  [Orb] [轨道位置、速度信息]
**/
void QuatJ2O(double QuatO[4], double QuatJ[4], double Orb[6]);

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
void AngRateDiff(double dw[3], double w[3], double J[3][3], double JI[3][3], double h[3], double T[3]);


//calTransMtx	转换矩阵计算

/**
* [FunctionName] [根据轨道信息获得惯性系到轨道系的姿态矩阵]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [MtxJtoO] [惯性系到轨道系的姿态矩阵]
* [input]  [orbit] [轨道信息]
**/
void MtxJtoOGetG(double MtxJtoO[3][3], double orbit[6]);

/**
* [FunctionName] [根据设定的321欧拉角，计算从轨道系到本体系的姿态矩阵]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [MtxOtoB] [解算后的轨道系到本体系的转换矩阵]
* [input]  [angle] [欧拉角rad]
**/
void MtxOtoBGet(double MtxOtoB[3][3], double angle[3]);

/**
* [FunctionName] [实时计算从轨道系到本体系的姿态矩阵]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [MtxOtoB] [由轨道系到本体系的转换矩阵]
* [input]  [angle] [姿态角]
**/
void MtxOtoBGet_realtime(double MtxOtoB[3][3], double angle[3]);

/**
* [FunctionName] [获得J2000到WGS84坐标的转换矩阵]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [JtoWGS] [J2000到WGS84坐标的转换矩阵]
* [input]  [JD] [儒略日]
* [input]  [date] [当前日期]
**/
void cordMtxJToWGSGetG(double JtoWGS[3][3], double JD, iTime_struct date);

/**
* [FunctionName] [获得J2000到84坐标的转换矩阵]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [JtoWGS] [J2000到84坐标的转换矩阵]
* [input]  [JToFJ] [J2000到瞬时平赤道坐标系的转换矩阵]
* [input]  [tTime] [儒略日]
**/
void cordMtxJToWGSGet_Jin(double JtoWGS[3][3], double JToFJ[3][3], double tTime);

/**
* [FunctionName] [获得地固系下的位置矢量]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [PosInWGS] [卫星在地固系下的位置矢量]
* [input]  [JToWGS] [惯性系到地固系的转移矩阵]
* [input]  [orbit] [位置速度信息]
**/
void posInWGSGetG(double PosInWGS[3], double JToWGS[3][3], double orbit[6]);

/**
* [FunctionName] [获得大地系下的位置表示]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [GeoCord] [大地系下的位置表示]
* [input]  [PosInWGS] [地固系下的位置矢量]
**/
void geoInfoGetG(double GeoCord[3], double PosInWGS[3]);

/**
* [FunctionName] [计算从惯性系到本体系的转换矩阵]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [JToB] [惯性系到本体系的转移矩阵]
* [input]  [attitude] [卫星姿态信息]
**/
void calJToB(double JToB[3][3], double attitude[7]);

//calTime	时间转换函数库

/**
* [FunctionName] [时间初始化]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [inSecond] [初始化需要的秒]
* [input]  [inUSecond] [初始化需要的微秒]
**/
double cal_JD_DateG(iTime_struct *newDate, iTime_struct refDate, int inSecond, int inUSecond);

/**
* [FunctionName] [计算儒略日]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [返回值] [儒略日]
* [input]  [tTime] [计算后的儒略日]
* [input]  [iTime] [当前日期时间]
**/
double calJD(double *tTime, iTime_struct iTime);


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
void calOrbG(double orbInfo[6], double kplInfo[6]);

/**
* [FunctionName] [轨道六根数计算位置和速度]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [orbInfo] [位置与速度]
* [input]  [kplInfo] [轨道六根数]
**/
void orbProp_VR(double newOrbit[6], double lstOrbit[6], double ha, double accel[3]);

/**
* [FunctionName] [位置和速度计算轨道六根数]
* [Author] [吴佳奇]
* [E-mail] [wjqjump@163.com]
* [Data]   [2018年11月14日]

* [output] [kpl] [解算后的六根数]
* [input]  [orbit] [位置和速度]
**/
void calKPLG(double kpl[6], double orbit[6]);



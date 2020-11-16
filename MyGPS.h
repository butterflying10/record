#pragma once


#include <conio.h>
#include <stdlib.h>


#include "stdio.h"
#include "math.h"
#include "string.h"

//GPS网平差的类
class CCoGPS
{
public:

	CCoGPS();//构造函数

	virtual ~CCoGPS();//析构函数

	char *ResultFile;//输出文件路径

	


	char *Netname;//网名

	int Anumber;//同步区总数

	int Pnumber;//总点数

	int Vnumber;//基线向量总数

	int* dir0;//各同步区首向量在观测值中的序号
	int* dir1, *dir2;//向量的端点号
	int* AeraNumber;//各向量所属的向量组的编号


	char** Pname;//点名指针数组

	int* IsKnown;//标识点数组




	double *L;//全区观测的基线向量
	double** P;//权矩阵数组

	double* ATPA, * ATPL;

	double* XYZ;//现在用来存储读取坐标文件获取的坐标
	//存储秩亏自由网平差之后的所有点的XYZ坐标，但这个坐标并不具备实际的意义，它只是以网的重心当作
	//坐标系的原点

	double* V;//存储秩亏自由网平差之后残差

	double Sigma;//验后单位权中误差




public :

	bool InputData(char* DataFile);//输入原始数据,观测值向量

	bool InputP_Data(char* PointDataFile);//输入坐标文件，只是坐标的近似文件，没有基准点



	void PrintData();//打印原始数据

	void Free();//秩亏自由网平差计算

	void PrintXYZ();

	void PrintLV();


	



private:
	int GetStationNumber(char* buff);//保存点名，返回点号


	void CaATPA();//组成法方程
	
	void CaATPAi(int nk, int* dir1, int* dir2, double* Pk,double *Lk);
	/*void CaL(double V[]);*/

	double Ca_dX(bool isPrint);//计算未知参数的改正数

	void CaLV(double V[]);

	double CaVTPV();//计算VTPV，然后用它算单位权中误差

};

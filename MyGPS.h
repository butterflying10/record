#pragma once


#include <conio.h>
#include <stdlib.h>


#include "stdio.h"
#include "math.h"
#include "string.h"

//GPS��ƽ�����
class CCoGPS
{
public:

	CCoGPS();//���캯��

	virtual ~CCoGPS();//��������

	char *ResultFile;//����ļ�·��

	


	char *Netname;//����

	int Anumber;//ͬ��������

	int Pnumber;//�ܵ���

	int Vnumber;//������������

	int* dir0;//��ͬ�����������ڹ۲�ֵ�е����
	int* dir1, *dir2;//�����Ķ˵��
	int* AeraNumber;//������������������ı��


	char** Pname;//����ָ������

	int* IsKnown;//��ʶ������




	double *L;//ȫ���۲�Ļ�������
	double** P;//Ȩ��������

	double* ATPA, * ATPL;

	double* XYZ;//���������洢��ȡ�����ļ���ȡ������
	//�洢�ȿ�������ƽ��֮������е��XYZ���꣬��������겢���߱�ʵ�ʵ����壬��ֻ�����������ĵ���
	//����ϵ��ԭ��

	double* V;//�洢�ȿ�������ƽ��֮��в�

	double Sigma;//���λȨ�����




public :

	bool InputData(char* DataFile);//����ԭʼ����,�۲�ֵ����

	bool InputP_Data(char* PointDataFile);//���������ļ���ֻ������Ľ����ļ���û�л�׼��



	void PrintData();//��ӡԭʼ����

	void Free();//�ȿ�������ƽ�����

	void PrintXYZ();

	void PrintLV();


	



private:
	int GetStationNumber(char* buff);//������������ص��


	void CaATPA();//��ɷ�����
	
	void CaATPAi(int nk, int* dir1, int* dir2, double* Pk,double *Lk);
	/*void CaL(double V[]);*/

	double Ca_dX(bool isPrint);//����δ֪�����ĸ�����

	void CaLV(double V[]);

	double CaVTPV();//����VTPV��Ȼ�������㵥λȨ�����

};

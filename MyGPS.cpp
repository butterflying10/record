

#include "MyGPS.h"
#include "public.h"
#include <iostream>
#include<fstream>
#include <string>
#include < iomanip >//����С�����λ��
//#include <stdlib.h>
//#include <stdio.h>

using namespace std;

//���캯��
CCoGPS::CCoGPS()
{
	
}

//��������

CCoGPS :: ~CCoGPS()
{

}

//�������뺯��
//����ԭʼ����
//����ȫ���Ĺ۲�����

bool CCoGPS::InputData(char* file)
{
	ifstream infile(file,ios::in|ios::_Nocreate);


	if (!infile)
	{
		cerr << "��ԭʼ����ʧ��" << endl;
		return 0;
	}
	//��ȡ����
	char netname[20];
	infile >> netname;
	Netname = netname;
	cout << Netname << endl;
	//��ȡ��ע
	string note;
	infile >> note;
	//��ȡͬ��������
	infile >> Anumber;
	cout << Anumber << endl;
	infile >> note;
	//��ȡ��������
	infile >> Vnumber;
	cout << Vnumber << endl;
	infile >> note;
	//��ȡȫ���ܵ���
	infile >> Pnumber;
	cout << Pnumber << endl;
	infile >> note;

	int t = 3 * Pnumber;
	int tt = t * (t + 1) / 2;
	
	ATPA = new double[tt]; //������ϵ�����󣨻������   �������ǽ��д洢   t*t�ľ���
	ATPL = new double[t];  //����������������


	XYZ = new double[t];//��������
	


	Pname = new char* [Pnumber];//������ַ����
	for (int i = 0; i < Pnumber; i++)Pname[i] = NULL;



	dir0 = new int[Anumber + 1];//��ͬ�����������ڹ۲�ֵ�е����
	dir1 = new int[Vnumber];//��������ʼ���
	dir2 = new int[Vnumber];//�����Ľ������
	AeraNumber = new int[Vnumber];//������������������ı��

	int n=Vnumber*3;//�۲�ֵ������

	P = new double* [Anumber]; //�۲�ֵȨ�����ַ

	L = new double[n];//�����Ĺ۲�ֵ

	V = new double[n];//�����۲�ֵ�Ĳв�

	dir0[0] = 0;//��һ��ͬ�������������ڹ۲�ֵ�е����Ϊ0
	for (int i = 0; i <= Anumber - 1; i++)
	{
		int d;//ͬ�����ı��
		infile >> d;
		cout << "ͬ�����ı��" << d << endl;
		infile >> note;//    ��ע"//ͬ�������"
		if (d != i + 1)
		{
			cout << "ͬ�����ı�Ų�����" << endl;
			return 0;
		}

		int jd;//ͬ�����ڵ�������
		infile >> jd;
		infile >> note;//    ��ע"//������������"
		/*ע�⣺
		infile >> note   �����������������/n  ������ʽ�е�getline(infile, note);�Ľ���������/n
		������Ҫ�ٶ�ȡһ�У���ĩβ��/n�������ٶ�ȡ��һ��
		*/
		getline(infile, note);

		dir0[i + 1] = dir0[i] + jd;//��һ��ͬ�������������ڹ۲�ֵ�е����

		getline(infile, note);// ��ע���������� D:\��ҵ���\ofile2002\1\oxisma.001  //   ����Դ�ļ�
		cout <<"fileinformation"<< note << endl;

		getline(infile, note);//��ע����������2002   1   1   0    0    0.000   //   �۲�ʱ��
		cout <<"observeTime"<< note << endl;
		getline(infile, note);
		getline(infile, note);
		getline(infile, note);//   ��ע��������//     �۲�����
		cout<<"observe Vertor" << note << endl;
	
		//��ȡ���ͬ�����ڵ������۲�ֵ��ÿ����������ʼ�����ֹ��
		for (int j = dir0[i]; j < dir0[i + 1]; j++)
		{
			char P1[20], P2[20];
			infile >> P1;
			infile >> P2;

			//cout << P1 <<"   "<<P2<< endl;
			dir1[j] = GetStationNumber(P1);//��j�������۲�ֵ������ʼ����P1
			dir2[j] = GetStationNumber(P2);//��j�������۲�ֵ������ֹ����P2

			if (dir1[j] < 0 || dir2[j] < 0)  
			{
				cout << "��ȡ�����۲�ֵ������Ŵ����ܵ���";
				return 0;
			}
			//��ȡ�����۲�ֵ

			infile >> L[3 * j];

			infile>>L[3 * j + 1];
			infile >> L[3 * j + 2];

			cout << P1 << "   " << P2<<"  "<< L[3 * j] <<"  "<< L[3 * j + 1] <<"  "<< L[3 * j + 2]<< endl;

			AeraNumber[j] = i;//j�����۲�ֵ���ڵ�һ��ͬ����
			getline(infile, note);
		}
		//----------------------------------------------------------
		getline(infile, note);
		getline(infile, note);
		getline(infile, note);
		//��ȡ����Э�������

		int ni = 3 * jd;

		double *Pi = new double[(ni + 1) * ni / 2];

		P[i] = Pi;//��i��ͬ���������۲�ֵ�ķ���Э������
		int s = 0;
		for (int j = 0; j < ni; j++)
		{
			int d;
			infile >> d;; //�к�
			for (int k = 0; k <d; k++)
			{
				
				infile >> Pi[s];//�����ķ���Э����
				//cout<< Pi[s] <<endl;
				s++;
			}
		}


	}
	infile.close();
	return 1;

}

bool CCoGPS::InputP_Data(char* file)
{
	ifstream infile(file, ios::in | ios::_Nocreate);

	if (!infile)
	{
		cerr << "��ԭʼ����ʧ��" << endl;
		return 0;
	}
	//��ȡ����
	/*char netname[20];
	infile >> netname;
	Netname = netname;
	cout << Netname << endl;*/
	//��ȡȫ���ܵ���,Ӧ�ú������۲�ֵ�ļ��е�ȫ���ܵ���һ�²Ŷ�

	int pnumber;
	infile >> pnumber;
	if (pnumber != Pnumber)
	{
		cout << "�����ļ�����������۲�ֵ�ļ���ȫ���ܵ�����һ��" << endl;
		//Pnumber = pnumber;

	}
	IsKnown = new int[Pnumber];

	string  note;//��ע
	infile >> note;//  " //ȫ���ܵ���"
	getline(infile, note);//  ��ȡ��һ��ĩβ�Ļ��з�"/n"
	infile >> note;//      "XYZ"
	getline(infile, note);

	for (int k = 0; k < Pnumber; k++)
	{
		cout << Pname[k] << "  ";
	}


	for (int i = 0; i < pnumber; i++)    //pnumber��34��  0-33��34����
	{
		int xuhao;
		infile >> xuhao;
		char dianming[20];
		infile >> dianming;
		int flag = 0;
		for (int j = 0; j < Pnumber; j++)
		{
		
			if ( strcmp(dianming , Pname[j])==0)//�����ַ����Ƚϴ�С����ȷ���0
			{
				flag = 1;
				infile >> XYZ[j * 3];
				infile >> XYZ[j * 3 + 1];
				infile >> XYZ[j * 3 + 2];
	
				int isknown;//��ʶ   ��Ϊ0ʱ˵�����������ǽ������꣬���߱���Ϊ��׼��
				//��Ϊ1ʱ��˵����������Ǿ��Ⱥܸߵ����꣬������Ϊ��׼
				infile >> isknown;
				IsKnown[j] = isknown;

				cout << Pname[j] <<"         "<<j<< "  " << XYZ[j * 3] << "  " << XYZ[j * 3 + 1] << "  " << XYZ[j * 3 + 2] << endl;

			}
			if (j==Pnumber-1&&flag==0)
			{
				cout << dianming << "�˵��ڻ��������ļ���δ�õ�" << endl;
				getline(infile, note);//��ȡ��һ�в�Ҫ��
				
			}
		}
	
	}

	infile.close();

	return 1;
}

void CCoGPS::Free()
{
//���۲�ֵ�ķ���Э�������棬��ΪȨ����
	for (int i = 0; i < Anumber; i++)
	{
	
		int ni = 3 * (dir0[i + 1] - dir0[i]);//���ͬ�����۲�ֵ������
		//����Э���������ʵ��һ���ֿ����
		//�ֿ�������棿����
		//�ֿ����Ϊ������������
		inverse(P[i], ni);
	}

	//����ATPA
	CaATPA();
	//�����ʵ���Ǽ�����SST��           ��STS�������ǵ���*��λ�󣬳��Ե����γɵ�λ�󣬣�����ô
	//Ӧ����S�е�Ԫ�س���Sqrt(����),,�������������STS��Ϊ��λ���Ǵ�ʱSST�еķ���Ԫ��Ϊ
	//1/��������
	for (int i = 0; i < 3 * Pnumber; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			//����ط�Ҳ�����i��j,��3�����������Ϊ0�Ļ���x��λ��
			if (i % 3 == j % 3) ATPA[ij(i, j)] += 1.0 / Pnumber;
		}
	}

	Ca_dX(true);//����������Ѿ�����������

	//����в�
	CaLV(V);

	//����VTPV
	double VTPV=CaVTPV();
	//�������λȨ�����
	//n+3-t,n�۲�ֵ������tδ֪��������
	Sigma = sqrt(VTPV / (3 * Vnumber + 3 - 3 * Pnumber));

	PrintXYZ();
	PrintLV();


	ofstream outfile(ResultFile, ios::app);
	outfile << "\n\n ������ƽ��:\n[pvv]=" << VTPV;
	outfile << "\n���λȨ����� :��=��" << Sigma;
	outfile.close();


}
double CCoGPS::Ca_dX(bool isPrint)
{
	ofstream outfile(ResultFile, ios::out);
	int t = 3 * Pnumber;//δ֪�����ĸ���

	if (!inverse(ATPA, t))//��ATPA����
	{
		cout << "�����̲�����!" << endl;
		//�ر�����ļ�
		outfile.close();
		exit(0);
	}
	if (isPrint)
	{
		outfile << "\n           =================== ��������� ===================\n";

	}
	for (int i = 0; i < t; i++)
	{

		double xi = 0.0;
		for (int j = 0; j < t; j++)
			xi -= ATPA[ij(i, j)] * ATPL[j];  //ATPA�����Ѿ���ATPA�������ˣ�����x=ATPA_ni * ATPL
											//һЩ������x=-ATPA_ni * ATPL �����ǹ���ĺ���ģ�ͳ����ͬ�����ؾ���
											//t*t���󣬣�ATPL��t*1�ľ���������ˣ�����Ҫ������Ͳ���

		XYZ[i] += xi;
		if (isPrint)
		{
			if (i % 3 == 0)
			{
				outfile << "\n"<<Pname[i / 3]<<setw(10);
			}
			outfile << xi<<setw(20);
		}
	}
	outfile.close();
}

/// <summary>
/// �������ƽ��ֵ�������
/// </summary>
void CCoGPS::PrintXYZ()
{
	ofstream outfile(ResultFile, ios::app);
	//outfile.precision(6);

	outfile << "\n\nSTATION       X                         Y                         Z                       RMS_X        RMS_Y        RMS_Z";
	outfile << "\n              (m)                         (m)                      (m)                         (cm)         (cm)         (cm)\n";
	outfile << "---------------------------------------------------------------------------------\n";

	//XYZ�����ľ�������ƽ��ֵ
	for (int i = 0; i < Pnumber; i++)
	{
		//��������ţ�����
		outfile << i + 1 << "   " << Pname[i]<<"    ";
		double X = XYZ[3 * i ];
		double Y = XYZ[3 * i + 1];
		double Z = XYZ[3 * i + 2];



		outfile <<right<<setw(20)<< setiosflags(ios::fixed)<< X << setw(20) << Y << setw(20) << Z  << endl;
	
	}
	outfile.close();


}

void CCoGPS::PrintLV()
{
	ofstream outfile(ResultFile, ios::app);

	outfile << "\n\n        	         ===== V =====\n\n";
	outfile << "    ������  �� ��     �� ��      Vx         Vy         Vz\n";
	//Vnumber  �ܹ��������ĸ���
	for (int i = 0; i < Vnumber; i++)
	{
		int ni = AeraNumber[i];//�ж��������λ����һ��ͬ����
		outfile << ni + 1 << setw(10);

		outfile << Pname[dir1[i]] << setw(15) << Pname[dir2[i]] << setw(20);

		outfile << V[3 * i] << setw(20) << V[3 * i + 1] << setw(20) << V[3 * i + 2]<<endl;
	}

	
	
	outfile.close();


}





/// <summary>
/// ���㳣�����ƽ��֮��Ĳв�
/// </summary>
/// <param name="V"></param>
void CCoGPS::CaLV(double V[])
{
	for (int i = 0; i < Vnumber; i++)
	{
	// dir1 ��������ʼ�������   Vnumber
	//dir2  �����Ľ����������   Vnumber
		int k1 = dir1[i];
		int k2 = dir2[i];
		///cout << k1 << "     " << k2 << endl;


		V[3 * i] = XYZ[3 * k2] - XYZ[3 * k1] - L[i * 3];
		V[3 * i+1] = XYZ[3 * k2+1] - XYZ[3 * k1+1] - L[i * 3+1];
		V[3 * i+2] = XYZ[3 * k2+2] - XYZ[3 * k1+2] - L[i * 3+2];

		//cout << V[3 * i] << "  " << V[3 * i + 1] << "  " << V[3 * i + 2] << "		`��������   "<<i<< endl;

	}

}



void CCoGPS::CaATPA()
{

	int t = 3 * Pnumber;
	int tt = t * (t + 1) / 2;
	for (int i = 0; i < t; i++) ATPL[i] = 0.0;
	for (int i = 0; i < tt; i++) ATPA[i] = 0.0;

	double* l = new double[3 * Vnumber];//������

	CaLV(l); //���̳�����ļ���

	for (int i = 0; i < Anumber; i++) //����������ɷ�����
	{
		int s = dir0[i];//��i��ͬ�����������ڹ۲�ֵ�е����
		int n = dir0[i + 1] - s;//��i��ͬ�������ж��ٸ�����
		CaATPAi(3 * n, dir1 + s, dir2 + s, P[i], l + 3 * s);
	}

	delete[]l;
}


void CCoGPS::CaATPAi(int nk, //����������ڹ۲�ֵ������   3*Vnumber
	int* dir1, //����������ڸ���������ʼ���
	int* dir2,//����������ڸ���������ֹ���
	double* Pk,//����������Ȩ����
	double* Lk) // ���������ĳ�����
{
	for (int s1 = 0; s1 < nk; s1++)  //s1�����s1���۲�ֵ
	{
		int i1 = 3 * dir1[s1 / 3] + s1 % 3; //����Vk��ϵ��Ϊ-1��δ֪�����
		int i2 = 3 * dir2[s1 / 3] + s1 % 3; //����Vk��ϵ��Ϊ+1��δ֪�����

		for (int s2 = 0; s2 < nk; s2++)//s2�����s2���۲�ֵ
		{
			int j1 = 3 * dir1[s2 / 3] + s2 % 3;
			int j2 = 3 * dir2[s2 / 3] + s2 % 3;

			double p12 = Pk[ij(s1, s2)];
			//cout << s1 << "  " << s2 << "  " << p12 << endl;

			double l = Lk[s2];

			ATPL[i1] -= p12 * l;
			ATPL[i2] += p12 * l;

			//cout << "l:" << l << "   ATPL[  ]   " << i1 << "  " << ATPL[i1] << endl;



			if (i1 >= j1)ATPA[ij(i1, j1)] += p12;
			if (i1 >= j2)ATPA[ij(i1, j2)] -= p12;
			if (i2 >= j1)ATPA[ij(i2, j1)] -= p12;
			if (i2 >= j2)ATPA[ij(i2, j2)] += p12;
		}
	}
}
/// <summary>
/// ����VTPV
/// </summary>
/// <returns></returns>
double CCoGPS::CaVTPV()
{
	double VTPV = 0.0;
	
	for (int i = 0; i< Anumber; i++)
	{
		double* Pi = P[i]; //i����Ȩ����
		int k1 = 3 * dir0[i];//i���������ڹ۲�ֵ�е����
		int k2 = 3 * dir0[i + 1];//i+1���������ڹ۲�ֵ�е����

		

		//˼�����
		//V1TP1V1+V2TP2V2+V3TP2V3
		double* Vi = V + k1;//��ʾ  ��i��  �۲�ֵ��������ʼ��ַ
		int  ni = k2 - k1;//��ʾ  ��i��  �۲�ֵ�����ĳ���
		//�ȳ�PV������ǰ���   VT
		double *PiVi = new double[ni];

		double vtpvi = 0.0;//��i����VTPV���
		for (int j = 0; j < k2 - k1; j++)
		{
			double pvj = 0.0;//P�ĵ�j�г���V�Ľ��
			for (int k = 0; k < k2 - k1; k++)
			{
				pvj += Pi[ij(j, k)] * Vi[k];
			
			}
			PiVi[j] = pvj;
			//����ǰ�����VT
			vtpvi += Vi[j] * pvj;
		}
		VTPV += vtpvi;
	}

	return VTPV;
}





int CCoGPS::GetStationNumber(char *buff)
{
	for (int i = 0; i < Pnumber; i++)
	{
		if (Pname[i] != NULL)//���ж��ǲ��ǿյģ��������λ�û�û����
		{
			//���Ѿ����������ַ����ĵ����Ƚ�
			if (strcmp(buff, Pname[i]) == 0)return i;
		}
		else
		{
			//���µ����浽�ڴ棬��ַ�ŵ�Pname������
			int len = strlen(buff);
			Pname[i] = new char[len + 1];
			strcpy_s(Pname[i],strlen(buff)+1, buff);
			return i;
		}
	}
	return -1;
}






#include "MyGPS.h"
#include "public.h"
#include <iostream>
#include<fstream>
#include <string>
#include < iomanip >//保留小数点后几位的
//#include <stdlib.h>
//#include <stdio.h>

using namespace std;

//构造函数
CCoGPS::CCoGPS()
{
	
}

//析构函数

CCoGPS :: ~CCoGPS()
{

}

//数据输入函数
//输入原始数据
//输入全区的观测向量

bool CCoGPS::InputData(char* file)
{
	ifstream infile(file,ios::in|ios::_Nocreate);


	if (!infile)
	{
		cerr << "打开原始数据失败" << endl;
		return 0;
	}
	//读取网名
	char netname[20];
	infile >> netname;
	Netname = netname;
	cout << Netname << endl;
	//读取备注
	string note;
	infile >> note;
	//读取同步区总数
	infile >> Anumber;
	cout << Anumber << endl;
	infile >> note;
	//读取向量总数
	infile >> Vnumber;
	cout << Vnumber << endl;
	infile >> note;
	//读取全区总点数
	infile >> Pnumber;
	cout << Pnumber << endl;
	infile >> note;

	int t = 3 * Pnumber;
	int tt = t * (t + 1) / 2;
	
	ATPA = new double[tt]; //法方程系数矩阵（或逆矩阵）   用下三角进行存储   t*t的矩阵
	ATPL = new double[t];  //法方程自由项向量


	XYZ = new double[t];//坐标数组
	


	Pname = new char* [Pnumber];//点名地址数据
	for (int i = 0; i < Pnumber; i++)Pname[i] = NULL;



	dir0 = new int[Anumber + 1];//各同步区首向量在观测值中的序号
	dir1 = new int[Vnumber];//向量的起始点号
	dir2 = new int[Vnumber];//向量的结束点号
	AeraNumber = new int[Vnumber];//各向量所属的向量组的编号

	int n=Vnumber*3;//观测值的总数

	P = new double* [Anumber]; //观测值权矩阵地址

	L = new double[n];//向量的观测值

	V = new double[n];//向量观测值的残差

	dir0[0] = 0;//第一个同步区的首向量在观测值中的序号为0
	for (int i = 0; i <= Anumber - 1; i++)
	{
		int d;//同步区的编号
		infile >> d;
		cout << "同步区的编号" << d << endl;
		infile >> note;//    备注"//同步区序号"
		if (d != i + 1)
		{
			cout << "同步区的编号不连续" << endl;
			return 0;
		}

		int jd;//同步区内的向量数
		infile >> jd;
		infile >> note;//    备注"//区内向量总数"
		/*注意：
		infile >> note   这个结束条件并不是/n  ，而下式中的getline(infile, note);的结束条件是/n
		所以需要再读取一行，把末尾的/n消掉，再读取下一行
		*/
		getline(infile, note);

		dir0[i + 1] = dir0[i] + jd;//下一个同步区的首向量在观测值中的序号

		getline(infile, note);// 备注：：：：： D:\毕业设计\ofile2002\1\oxisma.001  //   数据源文件
		cout <<"fileinformation"<< note << endl;

		getline(infile, note);//备注：：：：：2002   1   1   0    0    0.000   //   观测时间
		cout <<"observeTime"<< note << endl;
		getline(infile, note);
		getline(infile, note);
		getline(infile, note);//   备注：：：：//     观测向量
		cout<<"observe Vertor" << note << endl;
	
		//读取这个同步区内的向量观测值，每个向量的起始点和终止点
		for (int j = dir0[i]; j < dir0[i + 1]; j++)
		{
			char P1[20], P2[20];
			infile >> P1;
			infile >> P2;

			//cout << P1 <<"   "<<P2<< endl;
			dir1[j] = GetStationNumber(P1);//第j个向量观测值它的起始点是P1
			dir2[j] = GetStationNumber(P2);//第j个向量观测值它的终止点是P2

			if (dir1[j] < 0 || dir2[j] < 0)  
			{
				cout << "读取向量观测值出错：点号大于总点数";
				return 0;
			}
			//读取向量观测值

			infile >> L[3 * j];

			infile>>L[3 * j + 1];
			infile >> L[3 * j + 2];

			cout << P1 << "   " << P2<<"  "<< L[3 * j] <<"  "<< L[3 * j + 1] <<"  "<< L[3 * j + 2]<< endl;

			AeraNumber[j] = i;//j向量观测值属于第一个同步区
			getline(infile, note);
		}
		//----------------------------------------------------------
		getline(infile, note);
		getline(infile, note);
		getline(infile, note);
		//读取方差协方差矩阵

		int ni = 3 * jd;

		double *Pi = new double[(ni + 1) * ni / 2];

		P[i] = Pi;//第i个同步区向量观测值的方差协方差阵
		int s = 0;
		for (int j = 0; j < ni; j++)
		{
			int d;
			infile >> d;; //行号
			for (int k = 0; k <d; k++)
			{
				
				infile >> Pi[s];//向量的方差协方差
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
		cerr << "打开原始数据失败" << endl;
		return 0;
	}
	//读取网名
	/*char netname[20];
	infile >> netname;
	Netname = netname;
	cout << Netname << endl;*/
	//读取全区总点数,应该和向量观测值文件中的全区总点数一致才对

	int pnumber;
	infile >> pnumber;
	if (pnumber != Pnumber)
	{
		cout << "坐标文件与基线向量观测值文件的全区总点数不一致" << endl;
		//Pnumber = pnumber;

	}
	IsKnown = new int[Pnumber];

	string  note;//备注
	infile >> note;//  " //全区总点数"
	getline(infile, note);//  读取第一行末尾的换行符"/n"
	infile >> note;//      "XYZ"
	getline(infile, note);

	for (int k = 0; k < Pnumber; k++)
	{
		cout << Pname[k] << "  ";
	}


	for (int i = 0; i < pnumber; i++)    //pnumber是34个  0-33共34个点
	{
		int xuhao;
		infile >> xuhao;
		char dianming[20];
		infile >> dianming;
		int flag = 0;
		for (int j = 0; j < Pnumber; j++)
		{
		
			if ( strcmp(dianming , Pname[j])==0)//两个字符串比较大小，相等返回0
			{
				flag = 1;
				infile >> XYZ[j * 3];
				infile >> XYZ[j * 3 + 1];
				infile >> XYZ[j * 3 + 2];
	
				int isknown;//标识   当为0时说明这个坐标就是近似坐标，不具备作为基准，
				//当为1时，说明这个坐标是精度很高的坐标，可以作为基准
				infile >> isknown;
				IsKnown[j] = isknown;

				cout << Pname[j] <<"         "<<j<< "  " << XYZ[j * 3] << "  " << XYZ[j * 3 + 1] << "  " << XYZ[j * 3 + 2] << endl;

			}
			if (j==Pnumber-1&&flag==0)
			{
				cout << dianming << "此点在基线向量文件中未用到" << endl;
				getline(infile, note);//读取这一行不要的
				
			}
		}
	
	}

	infile.close();

	return 1;
}

void CCoGPS::Free()
{
//将观测值的方差协方差求逆，化为权矩阵
	for (int i = 0; i < Anumber; i++)
	{
	
		int ni = 3 * (dir0[i + 1] - dir0[i]);//这个同步区观测值的总数
		//方差协方差矩阵其实是一个分块矩阵
		//分块矩阵求逆？？？
		//分块矩阵为方阵，所以满足
		inverse(P[i], ni);
	}

	//计算ATPA
	CaATPA();
	//这个其实就是加上了SST，           而STS出来就是点数*单位阵，除以点数形成单位阵，，，那么
	//应该让S中的元素除以Sqrt(点数),,这样构造出来的STS才为单位阵。那此时SST中的非零元素为
	//1/（点数）
	for (int i = 0; i < 3 * Pnumber; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			//这个地方也很巧妙，i，j,对3的余数相等且为0的话在x的位置
			if (i % 3 == j % 3) ATPA[ij(i, j)] += 1.0 / Pnumber;
		}
	}

	Ca_dX(true);//这里坐标就已经被改正过了

	//计算残差
	CaLV(V);

	//计算VTPV
	double VTPV=CaVTPV();
	//计算验后单位权中误差
	//n+3-t,n观测值数量，t未知参数数量
	Sigma = sqrt(VTPV / (3 * Vnumber + 3 - 3 * Pnumber));

	PrintXYZ();
	PrintLV();


	ofstream outfile(ResultFile, ios::app);
	outfile << "\n\n 自由网平差:\n[pvv]=" << VTPV;
	outfile << "\n验后单位权中误差 :σ=±" << Sigma;
	outfile.close();


}
double CCoGPS::Ca_dX(bool isPrint)
{
	ofstream outfile(ResultFile, ios::out);
	int t = 3 * Pnumber;//未知参数的个数

	if (!inverse(ATPA, t))//对ATPA求逆
	{
		cout << "法方程不满秩!" << endl;
		//关闭输出文件
		outfile.close();
		exit(0);
	}
	if (isPrint)
	{
		outfile << "\n           =================== 坐标改正数 ===================\n";

	}
	for (int i = 0; i < t; i++)
	{

		double xi = 0.0;
		for (int j = 0; j < t; j++)
			xi -= ATPA[ij(i, j)] * ATPL[j];  //ATPA现在已经是ATPA的逆阵了，这里x=ATPA_ni * ATPL
											//一些书上是x=-ATPA_ni * ATPL ，这是构造的函数模型常数项不同，不必纠结
											//t*t矩阵，，ATPL是t*1的矩阵，两者相乘，还需要进行求和操作

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
/// 输出坐标平差值和中误差
/// </summary>
void CCoGPS::PrintXYZ()
{
	ofstream outfile(ResultFile, ios::app);
	//outfile.precision(6);

	outfile << "\n\nSTATION       X                         Y                         Z                       RMS_X        RMS_Y        RMS_Z";
	outfile << "\n              (m)                         (m)                      (m)                         (cm)         (cm)         (cm)\n";
	outfile << "---------------------------------------------------------------------------------\n";

	//XYZ里面存的就是坐标平差值
	for (int i = 0; i < Pnumber; i++)
	{
		//输出点的序号，点名
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
	outfile << "    向量组  点 名     点 名      Vx         Vy         Vz\n";
	//Vnumber  总共的向量的个数
	for (int i = 0; i < Vnumber; i++)
	{
		int ni = AeraNumber[i];//判断这个向量位于哪一个同步区
		outfile << ni + 1 << setw(10);

		outfile << Pname[dir1[i]] << setw(15) << Pname[dir2[i]] << setw(20);

		outfile << V[3 * i] << setw(20) << V[3 * i + 1] << setw(20) << V[3 * i + 2]<<endl;
	}

	
	
	outfile.close();


}





/// <summary>
/// 计算常数项还有平差之后的残差
/// </summary>
/// <param name="V"></param>
void CCoGPS::CaLV(double V[])
{
	for (int i = 0; i < Vnumber; i++)
	{
	// dir1 向量的起始点号数组   Vnumber
	//dir2  向量的结束点号数组   Vnumber
		int k1 = dir1[i];
		int k2 = dir2[i];
		///cout << k1 << "     " << k2 << endl;


		V[3 * i] = XYZ[3 * k2] - XYZ[3 * k1] - L[i * 3];
		V[3 * i+1] = XYZ[3 * k2+1] - XYZ[3 * k1+1] - L[i * 3+1];
		V[3 * i+2] = XYZ[3 * k2+2] - XYZ[3 * k1+2] - L[i * 3+2];

		//cout << V[3 * i] << "  " << V[3 * i + 1] << "  " << V[3 * i + 2] << "		`这个点号是   "<<i<< endl;

	}

}



void CCoGPS::CaATPA()
{

	int t = 3 * Pnumber;
	int tt = t * (t + 1) / 2;
	for (int i = 0; i < t; i++) ATPL[i] = 0.0;
	for (int i = 0; i < tt; i++) ATPA[i] = 0.0;

	double* l = new double[3 * Vnumber];//常数项

	CaLV(l); //误差方程常数项的计算

	for (int i = 0; i < Anumber; i++) //逐向量组组成法方程
	{
		int s = dir0[i];//第i个同步区首向量在观测值中的序号
		int n = dir0[i + 1] - s;//第i个同步区共有多少个向量
		CaATPAi(3 * n, dir1 + s, dir2 + s, P[i], l + 3 * s);
	}

	delete[]l;
}


void CCoGPS::CaATPAi(int nk, //这个向量组内观测值的总数   3*Vnumber
	int* dir1, //这个向量组内各向量的起始点号
	int* dir2,//这个向量组内各向量的终止点号
	double* Pk,//这个向量组的权矩阵
	double* Lk) // 这个向量组的常数项
{
	for (int s1 = 0; s1 < nk; s1++)  //s1代表第s1个观测值
	{
		int i1 = 3 * dir1[s1 / 3] + s1 % 3; //误差方程Vk中系数为-1的未知数编号
		int i2 = 3 * dir2[s1 / 3] + s1 % 3; //误差方程Vk中系数为+1的未知数编号

		for (int s2 = 0; s2 < nk; s2++)//s2代表第s2个观测值
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
/// 计算VTPV
/// </summary>
/// <returns></returns>
double CCoGPS::CaVTPV()
{
	double VTPV = 0.0;
	
	for (int i = 0; i< Anumber; i++)
	{
		double* Pi = P[i]; //i区的权矩阵
		int k1 = 3 * dir0[i];//i区首向量在观测值中的序号
		int k2 = 3 * dir0[i + 1];//i+1区首向量在观测值中的序号

		

		//思想就是
		//V1TP1V1+V2TP2V2+V3TP2V3
		double* Vi = V + k1;//表示  第i区  观测值向量的起始地址
		int  ni = k2 - k1;//表示  第i区  观测值向量的长度
		//先乘PV，再在前面乘   VT
		double *PiVi = new double[ni];

		double vtpvi = 0.0;//第i区的VTPV相乘
		for (int j = 0; j < k2 - k1; j++)
		{
			double pvj = 0.0;//P的第j行乘以V的结果
			for (int k = 0; k < k2 - k1; k++)
			{
				pvj += Pi[ij(j, k)] * Vi[k];
			
			}
			PiVi[j] = pvj;
			//再在前面乘以VT
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
		if (Pname[i] != NULL)//先判断是不是空的，就是这个位置还没存点号
		{
			//与已经存入点名地址数组的点名比较
			if (strcmp(buff, Pname[i]) == 0)return i;
		}
		else
		{
			//将新点名存到内存，地址放到Pname数组中
			int len = strlen(buff);
			Pname[i] = new char[len + 1];
			strcpy_s(Pname[i],strlen(buff)+1, buff);
			return i;
		}
	}
	return -1;
}




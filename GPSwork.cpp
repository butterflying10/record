// GPSwork.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#include <iostream>
#include "MyGPS.h"
#include "public.h"
#include "conio.h"
#include < iomanip >//保留小数点后几位的
using namespace std;



//自由网平差主函数
void main()
{
    CCoGPS coGPS;

    char* inputfile = "C:\\白腾飞\\研一课程学习\\平差\\基线向量52.txt";

    char* inputPfile= "C:\\白腾飞\\研一课程学习\\平差\\坐标文件.txt";

    char* resultfile = "C:\\白腾飞\\研一课程学习\\平差\\result52.txt";

    coGPS.ResultFile = resultfile;
    



    
    coGPS.InputData(inputfile);
    coGPS.InputP_Data(inputPfile);
   

    coGPS.Free();
}

//void main()
//{
//	double x12 = -807801.685499;
//	double y12 = 831845.667811;
//	double z12 = -2051577.244583;
//
//	double x1 = -1512848.531482;
//	double y1 = 4816556.788763;
//	double z1 = 3886841.239701;
//
//	
//	double x23 = 256110.677662;
//	double y23 = -809301.154800;
//	double z23 = 1759175.794680;
//
//	double x3 = -2064539.539397;
//	double y3 = 4839101.304744;
//	double z3 = 3594439.797728;
//
//
//	double x2 = (x1 + x12 + x3 - x23) / 2.0;
//	double y2 = (y1 + y3 + y12 - y23) / 2.0;
//	double z2 = (z1 + z3 + z12 - z23) / 2.0;
//
//	cout << setiosflags(ios::fixed) << x2 << "  " << y2 << "  " << z2 << endl;
//
//
//}
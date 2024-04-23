#pragma once

#include <string>
#include <vector>

using namespace std;

struct curline
{
	char line[256];
	static std::vector<double> GetLineByNumbers(const char* s);
};

class doserateTable
{
public:
	void Load(std::vector<std::vector<double>> data);
	void Clear();

	double getValue(double EFS, double depth);
	std::vector<std::vector<double>> Values;
};

class wedgesTable
{
public:
	void Load(std::vector<std::vector<double>> data);
	void Clear();

	double getValue(double EFS, double wangle);
	std::vector<std::vector<double>> Values;
};

class QATable
{
public:
	void Load(const char* fname);
	void Clear();
	void Dump();

	double getDouble(string line);

	doserateTable DRT;

	wedgesTable WT;

	int CalcPointNumber, secondpointnumber;

	//Данные пациента
	string Name, ID, Machine;
	
	//Точка расчета
	std::vector<string> Location;

	//Шапка таблицы
	std::vector<string> TableHeader;

	//Параметры полей
	std::vector<string> Field;
	std::vector<int> WAngle;
	std::vector<int> Energy;
	std::vector<double> fsx;
	std::vector<double> fsy;
	std::vector<double> efs;
	std::vector<double> depth;
	std::vector<double> dr;
	std::vector<double> edfw;
	std::vector<double> weight;
	std::vector<double> MU;
	std::vector<double> Dose;
	std::vector<string> comment;

	//Анализ
	std::vector<double> Sum;
	std::vector<double> TPS;
	std::vector<double> delta;
};
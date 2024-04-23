#include "dataReader.h"
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>  

using namespace std;

std::vector<double> curline::GetLineByNumbers(const char* s)
{
	string s1 = s;
	string num;
	std::vector<double> result;
	int length = s1.length();
	for (int i = 0; i < length; i++)
	{
		if (isdigit(s[i]) || s[i] == '.')
			num += s[i];
		else if (num != "")
		{
			result.push_back(stof(num));
			num.clear();
		}
		if (i == length - 1)
		{
			result.push_back(stof(num));
			num.clear();
		}
	}
	return result;
}

void doserateTable::Load(std::vector<std::vector<double>> data)
{
	Values = data;
}

void doserateTable::Clear()
{
	Values.clear();
}

double doserateTable::getValue(double EFS, double depth)
{
	double result = 0.0;
	int fsID = 0;
	for (fsID = 0; fsID < Values[0].size(); fsID++)
	{
		if (Values[0][fsID] > EFS)
			break;
	}
	int depthID = 1;
	for (depthID; depthID < Values.size(); depthID++)
	{
		if (Values[depthID][0] > depth)
			break;
	}
	if (fsID == Values[0].size() && depthID == Values.size())
		result = Values[depthID - 1][fsID];
	else if (fsID == Values[0].size())
	{
		result = Values[depthID - 1][fsID] + (Values[depthID][0] - depth) / (Values[depthID][0] - Values[depthID - 1][0]) * (Values[depthID][fsID] -
			Values[depthID - 1][fsID]);
	}
	else if (depthID == Values.size())
	{
		result = Values[depthID - 1][fsID] + (Values[0][fsID] - EFS) / (Values[0][fsID] - Values[0][fsID - 1]) * (Values[depthID - 1][fsID + 1] -
			Values[depthID - 1][fsID]);
	}
	else
	{
		result = (Values[0][fsID] - EFS) * (Values[depthID][0] - depth) / (Values[0][fsID] - Values[0][fsID - 1]) /
			(Values[depthID][0] - Values[depthID - 1][0]) * Values[depthID - 1][fsID] +
			(EFS - Values[0][fsID - 1]) * (Values[depthID][0] - depth) / (Values[0][fsID] - Values[0][fsID - 1]) /
			(Values[depthID][0] - Values[depthID - 1][0]) * Values[depthID - 1][fsID + 1] +
			(Values[0][fsID] - EFS) * (depth - Values[depthID - 1][0]) / (Values[0][fsID] - Values[0][fsID - 1]) /
			(Values[depthID][0] - Values[depthID - 1][0]) * Values[depthID][fsID] +
			(EFS - Values[0][fsID - 1]) * (depth - Values[depthID - 1][0]) / (Values[0][fsID] - Values[0][fsID - 1]) /
			(Values[depthID][0] - Values[depthID - 1][0]) * Values[depthID][fsID + 1];
	}
	return result;
}

void QATable::Load(const char* fname)
{
	CalcPointNumber = 0;
	int linenumber = 0;
	bool isPatientLoaded = false;
	bool isFirstPoint = false, isSecondPoint = false;
	ifstream in(fname);
	if (in.fail())
		throw exception("Can't open QA_Table file!");
	string line, s1;
	while (!in.fail())
	{
		getline(in, line);
		if (!isPatientLoaded)
		{
			line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
			line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
			line.erase(line.find("Patient"), 7);
			Name = line.substr(0, line.find("ID"));
			int i = 0;
			line.erase(0, line.find("ID") + 2);
			ID = line.substr(0, line.find("Machine"));
			line.erase(0, line.find("Machine") + 7);
			if (ID[0] == '0')
				ID = "_" + ID;
			Machine = line;
			isPatientLoaded = true;
		}
		else if (line.find("Calculation") != std::string::npos)
		{
			CalcPointNumber++;
			Location.push_back(line);
			for (int i = 0; i < 2; i++)
			{
				getline(in, line);
				Location.push_back(line);
			}
		}
		else if (line.find("Field") != std::string::npos)
		{
			if (isFirstPoint)
			{
				isSecondPoint = true;
				secondpointnumber = Field.size();
			}
			isFirstPoint = true;
			TableHeader.push_back(line);
			getline(in, line);
			while (line.find("T") != std::string::npos || line.find("V") != std::string::npos || line.find("Neck") != std::string::npos || line.find("L") != std::string::npos)
			{
				Field.push_back(line.substr(0, line.find("\t")));
				line.erase(0, line.find("\t") + 1);
				Energy.push_back(stoi(line.substr(0, line.find("\t"))));
				line.erase(0, line.find("\t") + 1);
				fsx.push_back(QATable::getDouble(line.substr(0, line.find("\t"))));
				line.erase(0, line.find("\t") + 1);
				fsy.push_back(QATable::getDouble(line.substr(0, line.find("\t"))));
				line.erase(0, line.find("\t") + 2);
				depth.push_back(10 * QATable::getDouble(line.substr(0, line.find("\t"))));
				line.erase(0, line.find("\t") + 3);
				weight.push_back(QATable::getDouble(line.substr(0, line.find("\t"))));
				line.erase(0, line.find("\t") + 1);
				MU.push_back(QATable::getDouble(line.substr(0, line.find("\t"))));
				getline(in, line);
			}
		}
		else if (line.find("Sum") != std::string::npos)
		{
			getline(in, line);
			line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
			line.erase(0, 3);
			TPS.push_back(QATable::getDouble(line));
		}
	}
	if (CalcPointNumber == 1 && isFirstPoint)
	{
		for (int i = 0; i < Field.size(); i++)
		{
			if (Field[i].find("W") != std::string::npos)
			{
				WAngle.push_back(stoi(Field[i].substr(Field[i].find("W") + 1, 2)));
			}
			else WAngle.push_back(0);
			efs.push_back(20 * fsx[i] * fsy[i] / (fsx[i] + fsy[i]));
			dr.push_back(DRT.getValue(efs[i], depth[i]));
			if (Energy[i] == 10)
				dr[i] += 0.05;
			else if (Energy[i] == 15)
				dr[i] += 0.1;
			edfw.push_back(WT.getValue(efs[i], WAngle[i]));
			Dose.push_back(dr[i] * edfw[i] * weight[i] * MU[i]);
			if (Dose[i] == 0)
				comment.push_back("Calculation point has no weight");
			else comment.push_back("");
		}
		Sum.push_back(0.0);
		for (int i = 0; i < Field.size(); i++)
			Sum[0] += Dose[i];
		delta.push_back((1 - Sum[0] / TPS[0]) * 100);
	}
	else if (CalcPointNumber == 2 && isSecondPoint)
	{
		for (int i = 0; i < Field.size(); i++)
		{
			if (Field[i].find("W") != std::string::npos)
			{
				WAngle.push_back(stoi(Field[i].substr(Field[i].find("W") + 1, 2)));
			}
			else WAngle.push_back(0);
			efs.push_back(20 * fsx[i] * fsy[i] / (fsx[i] + fsy[i]));
			dr.push_back(DRT.getValue(efs[i], depth[i]));
			if (Energy[i] == 10)
				dr[i] += 0.05;
			else if (Energy[i] == 15)
				dr[i] += 0.1;
			edfw.push_back(WT.getValue(efs[i], WAngle[i]));
			Dose.push_back(dr[i] * edfw[i] * weight[i] * MU[i]);
			if (Dose[i] == 0)
				comment.push_back("Calculation point has no weight");
			else comment.push_back("");
		}
		Sum.push_back(0.0);
		for (int i = 0; i < Field.size() - secondpointnumber; i++)
			Sum[0] += Dose[i];
		delta.push_back((1 - Sum[0] / TPS[0]) * 100);
		Sum.push_back(0.0);
		for (int i = Field.size() - secondpointnumber; i < Field.size(); i++)
			Sum[1] += Dose[i];
		delta.push_back((1 - Sum[1] / TPS[1]) * 100);
	}
	else throw exception("Program doesn't support more than 2 calculation points!");
	linenumber++;
}

void QATable::Clear()
{
	Name.clear();
	ID.clear();
	Machine.clear();
	Field.clear();
	Energy.clear();
	fsx.clear();
	fsy.clear();
	efs.clear();
	depth.clear();
	dr.clear();
	edfw.clear();
	weight.clear();
	MU.clear();
	Dose.clear();
	comment.clear();
	DRT.Clear();
	WT.Clear();

	Sum.clear();
	TPS.clear();
	delta.clear();
}

void QATable::Dump()
{
	ofstream out;
	string outpath = "Reports\\Report_" + Name + ".txt";
	out.open(outpath);
	if (out.is_open() && CalcPointNumber == 1)
	{
		out << "Patient\t" << Name << "\tID\t" << ID << "\tMachine\t" << Machine << "\n\n";
		out << Location[0] << endl << Location[1] << endl << Location[2] << endl << endl;
		out << TableHeader[0] << "\n";
		for (int i = 0; i < Field.size(); i++)
		{
			out << Field[i] << "\t" << Energy[i] << "\t" << fsx[i] << "\t" << fsy[i] << "\t" << efs[i] / 10 << "\t" << depth[i] / 10<< "\t" <<
				dr[i] << "\t" << edfw[i] << "\t" << weight[i] << "\t" << MU[i] << "\t" << Dose[i] << "\t" << comment[i];
			out << "\n";
		}
		out << "\n\n";
		out << "\t\t\t\t\t\t\t\tSum\t" << Sum[0] << "\n";
		out << "\t\t\t\t\t\t\t\tTPS\t" << TPS[0] << "\n";
		out << "\t\t\t\t\t\t\t\tdelta, %\t" << delta[0] << "\n";
		out.close();
	}
	else
	{
		out << "Patient\t" << Name << "\tID\t" << ID << "\tMachine\t" << Machine << "\n\n";
		out << Location[0] << endl << Location[1] << endl << Location[2] << endl << endl;
		out << TableHeader[0] << "\n";
		for (int i = 0; i < Field.size() - secondpointnumber; i++)
		{
			out << Field[i] << "\t" << Energy[i] << "\t" << fsx[i] << "\t" << fsy[i] << "\t" << efs[i] / 10 << "\t" << depth[i] / 10 << "\t" <<
				dr[i] << "\t" << edfw[i] << "\t" << weight[i] << "\t" << MU[i] << "\t" << Dose[i] << "\t" << comment[i];
			out << "\n";
		}
		out << "\n\n";
		out << "\t\t\t\t\t\t\t\tSum\t" << Sum[0] << "\n";
		out << "\t\t\t\t\t\t\t\tTPS\t" << TPS[0] << "\n";
		out << "\t\t\t\t\t\t\t\tdelta, %\t" << delta[0] << "\n";
		out << Location[3] << endl << Location[4] << endl << Location[5] << endl << endl;
		out << "\n" << TableHeader[1] << "\n";
		for (int i = Field.size() - secondpointnumber; i < Field.size(); i++)
		{
			out << Field[i] << "\t" << Energy[i] << "\t" << fsx[i] << "\t" << fsy[i] << "\t" << efs[i] / 10 << "\t" << depth[i] / 10 << "\t" <<
				dr[i] << "\t" << edfw[i] << "\t" << weight[i] << "\t" << MU[i] << "\t" << Dose[i] << "\t" << comment[i];
			out << "\n";
		}
		out << "\n\n";
		out << "\t\t\t\t\t\t\t\tSum\t" << Sum[1] << "\n";
		out << "\t\t\t\t\t\t\t\tTPS\t" << TPS[1] << "\n";
		out << "\t\t\t\t\t\t\t\tdelta, %\t" << delta[1] << "\n";
		out.close();
	}
}

double QATable::getDouble(string line)
{
	if (line.find(",") != std::string::npos)
		line.replace(line.find(","), 1, ".");
	double result = stof(line);
	return result;
}

void wedgesTable::Load(std::vector<std::vector<double>> data)
{
	Values = data;
}

void wedgesTable::Clear()
{
	Values.clear();
}

double wedgesTable::getValue(double EFS, double wangle)
{
	EFS /= 10;
	double result = 0.0;
	int fsID = 0;
	for (fsID = 0; fsID < Values[0].size(); fsID++)
	{
		if (Values[0][fsID] > EFS)
			break;
	}
	int wangleID = 1;
	for (wangleID; wangleID < Values.size(); wangleID++)
	{
		if (Values[wangleID][0] > wangle)
			break;
	}
	if (wangle != 0)
	{
		if (fsID == Values[0].size() && wangleID == Values.size())
			result = Values[wangleID - 1][fsID];
		else if (fsID == Values[0].size())
		{
			result = Values[wangleID - 1][fsID] + (Values[wangleID][0] - wangle) / (Values[wangleID][0] - Values[wangleID - 1][0]) * (Values[wangleID][fsID] -
				Values[wangleID - 1][fsID]);
		}
		else if (wangleID == Values.size())
		{
			result = Values[wangleID - 1][fsID] + (Values[0][fsID] - EFS) / (Values[0][fsID] - Values[0][fsID - 1]) * (Values[wangleID - 1][fsID + 1] -
				Values[wangleID - 1][fsID]);
		}
		else
		{
			result = (Values[0][fsID] - EFS) * (Values[wangleID][0] - wangle) / (Values[0][fsID] - Values[0][fsID - 1]) /
				(Values[wangleID][0] - Values[wangleID - 1][0]) * Values[wangleID - 1][fsID] +
				(EFS - Values[0][fsID - 1]) * (Values[wangleID][0] - wangle) / (Values[0][fsID] - Values[0][fsID - 1]) /
				(Values[wangleID][0] - Values[wangleID - 1][0]) * Values[wangleID - 1][fsID + 1] +
				(Values[0][fsID] - EFS) * (wangle - Values[wangleID - 1][0]) / (Values[0][fsID] - Values[0][fsID - 1]) /
				(Values[wangleID][0] - Values[wangleID - 1][0]) * Values[wangleID][fsID] +
				(EFS - Values[0][fsID - 1]) * (wangle - Values[wangleID - 1][0]) / (Values[0][fsID] - Values[0][fsID - 1]) /
				(Values[wangleID][0] - Values[wangleID - 1][0]) * Values[wangleID][fsID + 1];
		}
	}
	else result = 1.0;
	return result;
}

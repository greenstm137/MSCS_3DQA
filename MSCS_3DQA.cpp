// MSCS_3DQA.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <fstream>
#include <string>
#include "dataReader.h"
#include <filesystem>
#include <ctype.h>

namespace fs = std::filesystem;

using namespace std;

int main()
{
    curline curline_;
    string line;
    ifstream in("6X.dat");
    std::vector<std::vector<double>> data6x;
    while (!in.eof())
    {
        getline(in, line);
        ::memcpy(&curline_, line.c_str(), 256);
        data6x.push_back(curline::GetLineByNumbers(curline_.line));
    }
    in.close();
    ifstream wedges("EDWF.dat");
    std::vector<std::vector<double>> wedges_;
    while (!wedges.eof())
    {
        getline(wedges, line);
        ::memcpy(&curline_, line.c_str(), 256);
        wedges_.push_back(curline::GetLineByNumbers(curline_.line));
    }
    wedges.close();
	for (const auto& entry : fs::directory_iterator("QA_calc"))
	{
		if (!fs::path(entry.path()).has_stem() || !fs::path(entry.path()).has_extension())
			continue;

		string fname = fs::path(entry.path()).stem().string();
		string ext = fs::path(entry.path()).extension().string();
        QATable patient;
        if (ext == ".txt")
        {
            patient.DRT.Load(data6x);
            patient.WT.Load(wedges_);
            patient.Load(fs::path(entry.path()).string().c_str());
            patient.Dump();
            patient.Clear();
        }
	}
}
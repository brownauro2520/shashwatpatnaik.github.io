#ifndef READSAVEFILES_H_INCLUDED
#define READSAVEFILES_H_INCLUDED
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>
#include <fstream>
#include <sstream>

using namespace std;

void savefilev(vector<double> vec, string file_name) {
    ofstream outfile(file_name);
    int max = std::numeric_limits<double>::digits10;
    for (int i = 0; i < vec.size(); i++) {
        outfile << setprecision(max) << vec[i] << endl;
    }
    outfile.close();
    cout << "Vector saved to file: " << file_name << endl;
}

void savefilem(vector<vector<double>>& matrix, string fileName) {
    ofstream file;
    file.open(fileName);
    int max = std::numeric_limits<double>::digits10;
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            file << setprecision(max) <<  matrix[i][j] << " ";
        }
        file << endl;
    }
    file.close();
    cout << "Matrix saved to file: " << fileName << endl;
}


vector<vector<double>> readmatrix(string fileName) {
    vector<vector<double>> data;
    ifstream file(fileName);
    if (file.is_open()) {
        string line;
        while (getline(file, line)) {
            vector<double> row;
            stringstream ss(line);
            double value;
            while (ss >> value) {
                row.push_back(value);
            }
            data.push_back(row);
        }
        file.close();
    } else {
        cout << "Unable to open file." << fileName << endl;
    }
    return data;
}

vector<double> readvector(const string& fileName) {
    vector<double> values;
    ifstream file(fileName);
    if (file.is_open()) {
        double value;
        while (file >> value) {
            values.push_back(value);
        }
        file.close();
    } else {
        cerr << "Failed to open file " << fileName << endl;
    }
    return values;
}


#endif // READSAVEFILES_H_INCLUDED

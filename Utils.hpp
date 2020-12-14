//
// Created by Raffaele Montella on 13/12/20.
//

#ifndef WACOMMPLUSPLUS_UTILS_HPP
#define WACOMMPLUSPLUS_UTILS_HPP


#include <string>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

class Utils {
public:
    Utils();
    ~Utils();

    static string trim(const string &str, const string &whitespace);
    static string reduce(const string &str, const string &fill, const string &whitespace);
    static void tokenize(string &str, char delim, vector<string> &out);

};


#endif //WACOMMPLUSPLUS_UTILS_HPP

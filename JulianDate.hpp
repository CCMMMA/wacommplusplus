//
// Created by Raffaele Montella on 31/12/20.
//

#ifndef WACOMMPLUSPLUS_JULIANDATE_HPP
#define WACOMMPLUSPLUS_JULIANDATE_HPP

#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>

using namespace std;

class Calendar {
public:
    constexpr static const int YEAR = 0;
    constexpr static const int MONTH = 1;
    constexpr static const int DAY_OF_MONTH = 2;
    constexpr static const int HOUR_OF_DAY = 3;

    Calendar();
    Calendar(string value);
    Calendar(string format, string value);
    Calendar(int year, int month, int day, int hour);
    ~Calendar();

    int get(int idx);
    void set(int idx, int value);

    string asNCEPdate();
    string format(string format);
    void parse(string format, string value);

private:
    int _data[4];

};

class JulianDate {
    /**
* Returns the Julian day number that begins at noon of
* this day, Positive year signifies A.D., negative year B.C.
* Remember that the year after 1 B.C. was 1 A.D.
*
* ref :
*  Numerical Recipes in C, 2nd ed., Cambridge University Press 1992
*/
    // Gregorian Calendar adopted Oct. 15, 1582 (2299161)

public:
    constexpr static const int JGREG= 15 + 31*(10+12*1582);
    constexpr static const double HALFSECOND = 0.5;

    static double toModJulian(Calendar cal);
    static void fromModJulian(double modJulian, Calendar &cal);

    static double toModJulian(double julianRef, Calendar cal);
    static void fromModJulian(double julianRef, double modJulian, Calendar &cal);

    static double get19680523();
    static double get19700101();

    static double toJulian(Calendar cal);
    static void fromJulian(double injulian, Calendar &cal);

    static double toJulian(int y, int m, int d, int h);
    static void fromJulian(double injulian, int &y, int &m, int &d, int &h);

    static void playground();

};


#endif //WACOMMPLUSPLUS_JULIANDATE_HPP

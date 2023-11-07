//
// Created by Raffaele Montella on 31/12/20.
//

#include "JulianDate.hpp"
#include <math.h>
#include <chrono>
#include <ctime>

double JulianDate::toJulian(Calendar cal) {
    return toJulian(cal.get(Calendar::YEAR),
                    cal.get(Calendar::MONTH)+1,
                    cal.get(Calendar::DAY_OF_MONTH),
                    cal.get(Calendar::HOUR_OF_DAY),
                    cal.get(Calendar::MINUTES),
                    cal.get(Calendar::SECONDS));
}

double JulianDate::toJulian(int year, int month, int day, int hour, int minute, int second) {
    int julianYear = year;
    if (year < 0) julianYear++;
    int julianMonth = month;
    if (month > 2) {
        julianMonth++;
    }
    else {
        julianYear--;
        julianMonth += 13;
    }

    double julian = (floor(365.25 * julianYear)
                     + floor(30.6001*julianMonth) + day + 1720995.0);
    if (day + 31 * (month + 12 * year) >= JGREG) {
        // change over to Gregorian calendar
        int ja = (int)(0.01 * julianYear);
        julian += 2 - ja + (0.25 * ja);
    }
    //julian=floor(julian)+hour/24.+minute/1440.+second/86400.;
    julian=floor(julian)+hour/24.;

    return julian;
}

void JulianDate::fromJulian(double injulian, Calendar &cal) {
    int year, month, day, hour, minute, second;
    fromJulian(injulian,year, month, day, hour, minute, second);
    cal.set(Calendar::YEAR, year);
    cal.set(Calendar::MONTH, month-1);
    cal.set(Calendar::DAY_OF_MONTH, day);
    cal.set(Calendar::HOUR_OF_DAY, hour);
    cal.set(Calendar::MINUTES, minute);
    cal.set(Calendar::SECONDS, second);
}

/**
 * Converts a Julian day to a calendar date
 * ref :
 * Numerical Recipes in C, 2nd ed., Cambridge University Press 1992
 */
void JulianDate::fromJulian(double injulian, int &year, int &month, int &day, int &hour, int &minute, int &second) {
    int jalpha,ja,jb,jc,jd,je;
    double hf=injulian-floor(injulian);
    hour=(int)round(24*hf);
    double mf=24*hf - hour;
    minute=(int)round(60*mf);
    double sf=60*mf-minute;
    second=(int)round(60*sf);
    double julian = injulian + HALFSECOND / 86400.0;
    ja = (int) julian;
    if (ja>= JGREG) {
        jalpha = (int) (((ja - 1867216) - 0.25) / 36524.25);
        ja = ja + 1 + jalpha - jalpha / 4;
    }

    jb = ja + 1524;
    jc = (int) (6680.0 + ((jb - 2439870) - 122.1) / 365.25);
    jd = 365 * jc + jc / 4;
    je = (int) ((jb - jd) / 30.6001);
    day = jb - jd - (int) (30.6001 * je);
    month = je - 1;
    if (month > 12) month = month - 12;
    year = jc - 4715;
    if (month > 2) year--;
    if (year <= 0) year--;

}

double JulianDate::toModJulian(double julianRef, Calendar cal) {
    return toJulian(cal)-julianRef;
}

void JulianDate::fromModJulian(double julianRef, double modJulian, Calendar &cal) {
    double julian=modJulian+julianRef;
    fromJulian(julian, cal);
}

double JulianDate::toModJulian(Calendar cal) {
    return toModJulian(get19680523(), cal);
}

void JulianDate::fromModJulian(double modJulian, Calendar &cal) {
    fromModJulian(get19680523(),modJulian,cal);
}

double JulianDate::get19680523() {
    return toJulian(1968, 5, 23,0, 0,0 );
}

double JulianDate::get19700101() {
    return toJulian(1970, 1, 1,0,0,0 );
}

void JulianDate::playground() {
    // FIRST TEST reference point
    double jd=toJulian( 1968, 5, 23, 0,0,0 );
    cout << "Julian date for May 23, 1968 : " << jd << endl;
    // output : 2440000
    int year, month, day, hour, minute, second;
    jd=toJulian(1968, 5, 23,0,0,0 );
    fromJulian(jd, year, month, day, hour, minute, second);
    cout << "... back to calendar : " << year << " " << month << " " << day << " " << hour << " " << minute << " " << second << endl;


    // SECOND TEST today
    Calendar today;
    double todayJulian = toJulian (
            today.get(Calendar::YEAR),
            today.get(Calendar::MONTH)+1,
            today.get(Calendar::DAY_OF_MONTH),
            today.get(Calendar::HOUR_OF_DAY),
            today.get(Calendar::MINUTES),
            today.get(Calendar::SECONDS)
    );

    cout << "Julian date for today : " << todayJulian << endl;
    fromJulian(todayJulian, year,month, day, hour, minute, second);
    cout << "... back to calendar : " << year << " " << month << " " << day << " " << hour << endl;

    // THIRD TEST
    double date1 = toJulian(2005,1,1,0,0,0);
    double date2 = toJulian(2005,1,31,0,0,0);
    cout << "Between 2005-01-01 and 2005-01-31 : " <<
         (date2 - date1) << " days" << endl;

    /*
       expected output :
          Julian date for May 23, 1968 : 2440000.0
          ... back to calendar 1968 5 23
          Julian date for today : 2453487.0
          ... back to calendar 2005 4 26
          Between 2005-01-01 and 2005-01-31 : 30.0 days
    */
}

string Calendar::format(string format) {

    tm my_tm;
    my_tm.tm_year = _data[YEAR]-1900;
    my_tm.tm_mon = _data[MONTH];
    my_tm.tm_mday = _data[DAY_OF_MONTH];
    my_tm.tm_hour = _data[HOUR_OF_DAY];
    my_tm.tm_min = _data[MINUTES];
    my_tm.tm_sec = _data[SECONDS];
    char s[256];
    strftime(s,256,format.c_str(),&my_tm);
    return string(s);
}

void Calendar::parse(string format, string value) {
    tm my_tm;
    strptime(value.c_str(), format.c_str(), &my_tm);
    _data[YEAR]=my_tm.tm_year+1900;
    _data[MONTH] = my_tm.tm_mon;
    _data[DAY_OF_MONTH] = my_tm.tm_mday;
    _data[HOUR_OF_DAY] = my_tm.tm_hour;
    _data[MINUTES] = my_tm.tm_min;
    //_data[SECONDS] = my_tm.tm_sec;
    _data[SECONDS] = 0;
}

Calendar::Calendar() {
    auto now = std::chrono::system_clock::now();
    std::time_t dt = std::chrono::system_clock::to_time_t(now);
    auto dateTime = std::gmtime(&dt);

    _data[YEAR]=dateTime->tm_year+1900;
    _data[MONTH]=dateTime->tm_mon;
    _data[DAY_OF_MONTH]=dateTime->tm_mday;
    _data[HOUR_OF_DAY]=dateTime->tm_hour;
    _data[MINUTES]=dateTime->tm_min;
    _data[SECONDS]=dateTime->tm_sec;
}

Calendar::Calendar(string value) {
    parse("%Y%m%dZ%H", value);

}

Calendar::Calendar(string format, string value) {
    parse(format, value);
}

Calendar::Calendar(int year, int month, int day, int hour, int minute, int second) {
    _data[YEAR]=year;
    _data[MONTH]=month-1;
    _data[DAY_OF_MONTH]=day;
    _data[HOUR_OF_DAY]=hour;
    _data[MINUTES]=minute;
    _data[SECONDS]=second;
}

Calendar::~Calendar() = default;

int Calendar::get(int idx) { return _data[idx];}
void Calendar::set(int idx, int value) { _data[idx] = value;}

string Calendar::asNCEPdate() {
    return format("%Y%m%dZ%H");
}


//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#ifndef AGM3D_UTIL_H
#define AGM3D_UTIL_H

#include <cstdio>
#include <string>
#include <cstring>
#include <iostream>
#include <array>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <limits>
#include <sstream>
#include <algorithm>
#include <memory>
#include <iomanip>
#include <numeric>
#include <functional>
#include <chrono>
#include <ctime>
#include <omp.h>
#include <sys/resource.h>

constexpr double UNITVALUE{1.0000000000000000E0};
constexpr double HALFVALUE{5.0000000000000000E-1};
constexpr double ZEROVALUE{0.0000000000000000E0};
constexpr double NEARZERO{1.0000000000000000E-10};
constexpr int NT{10};

//#ifndef UNITVALUE
//#define UNITVALUE 1.0000000000000000E0
//#endif
//#ifndef HALFVALUE
//#define HALFVALUE 5.0000000000000000E-1
//#endif
//#ifndef ZEROVALUE
//#define ZEROVALUE 0.0000000000000000E0
//#endif
//#ifndef NEARZERO
//#define NEARZERO 1.0000000000000000E-10
//#endif
//#ifndef NT
//#define NT 10
//#endif

namespace AGM {
    auto isclose(double x, double y, double eps = NEARZERO) -> bool;

    void printError(const std::string &function_name);

    void printError(const char *function_name, const char *fmt, ...);

    auto iszero(double x, double eps = NEARZERO) -> bool;

    auto sgn(double d) -> double;

    auto ispositive(double d) -> bool;

    auto isnegative(double d) -> bool;

    auto ispositivezero(double d) -> bool;

    auto isnegativezero(double d) -> bool;

    enum STENCIL {
        R, L, F, B, U, D, RF, RB, RU, RD, LF, LB, LU, LD, FR, FL, FU, FD, BR, BL, BU, BD, UR, UL, UF, UB, DR, DL, DF, DB
//      0, 1, 2, 3, 4, 5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29
    };

    static std::unordered_map<std::string, int> valueMap{{"sol", 0},
                                                         {"phi", 1},
                                                         {"psi", 2},
                                                         {"rhs", 3},
                                                         {"bdv", 4},
                                                         {"dx",  5},
                                                         {"dy",  6},
                                                         {"dz",  7},
                                                         {"dxx", 8},
                                                         {"dyy", 9},
                                                         {"dzz", 10}};
}


#endif //AGM3D_UTIL_H

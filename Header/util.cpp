//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#include "util.h"

bool AGM::isclose(double x, double y, double eps) {
    return std::fabs(x - y) < eps;
}

void AGM::printError(const std::string &function_name) {
    std::cout << std::endl << function_name << std::endl;
    exit(1);
}

void AGM::printError(const char *function_name, const char *fmt, ...) {
    char buf[256] = {0,};
    va_list ap;

    printf("\n");
    printf("Fatal error has occur in %s\n", function_name);
    sprintf(buf, "Massage: ");

    va_start(ap, fmt);
    vsprintf(buf + strlen(buf), fmt, ap);
    va_end(ap);

    puts(buf);
    exit(1);
}

bool AGM::iszero(double x, double eps) {
    return std::fabs(x) < eps;
}

double AGM::sgn(double d) {
    return (d > ZEROVALUE) - (d < ZEROVALUE);
}

bool AGM::ispositive(double d) {
    return d > NEARZERO;
}

bool AGM::isnegative(double d) {
    return d < -NEARZERO;
}

bool AGM::ispositivezero(double d) {
    return ispositive(d) || iszero(d);
}

bool AGM::isnegativezero(double d) {
    return isnegative(d) || iszero(d);
}

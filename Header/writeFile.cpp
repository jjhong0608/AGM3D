//
// Created by NIMS-JUNHONG on 2022/01/07.
//

#include "writeFile.h"

template<typename T>
AGM::writeFile<T>::writeFile() = default;

template<typename T>
AGM::writeFile<T>::writeFile(std::vector<T> *pts) : pts(pts) {}

template<typename T>
auto AGM::writeFile<T>::getPt() const -> T * {
    return pt;
}

template<typename T>
void AGM::writeFile<T>::setPt(T *t) {
    writeFile::pt = t;
}

template<typename T>
auto AGM::writeFile<T>::getPts() const -> std::vector<T> * {
    return pts;
}

template<typename T>
void AGM::writeFile<T>::setPts(std::vector<T> *vector) {
    writeFile::pts = vector;
}

template<typename T>
auto AGM::writeFile<T>::calculateErrorAtPoint(const std::string &string) {
    auto printErrorToDouble = [&](const std::string &string1) -> double {
        std::cout << "\nidx = " << pt->getIdx() << ", (" << (*pt)[0] << ", " << (*pt)[1] << ", " << (*pt)[2] << ")\n";
        printError("AGM::writeFile<T>::calculateErrorAtPoint()", string1.c_str());
        return ZEROVALUE;
    };
    double xm = (*pt)[L] ? (*pt)[L]->getXyz()[0] : (*pt)[LF] ? (*pt)[LF]->getXyz()[0] : (*pt)[LU]
                                                                                        ? (*pt)[LU]->getXyz()[0]
                                                                                        : printErrorToDouble("xm");
    double xp = (*pt)[R] ? (*pt)[R]->getXyz()[0] : (*pt)[RF] ? (*pt)[RF]->getXyz()[0] : (*pt)[RU]
                                                                                        ? (*pt)[RU]->getXyz()[0]
                                                                                        : printErrorToDouble("xp");
    double ym = (*pt)[B] ? (*pt)[B]->getXyz()[1] : (*pt)[BR] ? (*pt)[BR]->getXyz()[1] : (*pt)[BU]
                                                                                        ? (*pt)[BU]->getXyz()[1]
                                                                                        : printErrorToDouble("ym");
    double yp = (*pt)[F] ? (*pt)[F]->getXyz()[1] : (*pt)[FR] ? (*pt)[FR]->getXyz()[1] : (*pt)[FU]
                                                                                        ? (*pt)[FU]->getXyz()[1]
                                                                                        : printErrorToDouble("yp");
    double zm = (*pt)[D] ? (*pt)[D]->getXyz()[2] : (*pt)[DR] ? (*pt)[DR]->getXyz()[2] : (*pt)[DF]
                                                                                        ? (*pt)[DF]->getXyz()[2]
                                                                                        : printErrorToDouble("zm");
    double zp = (*pt)[U] ? (*pt)[U]->getXyz()[2] : (*pt)[UR] ? (*pt)[UR]->getXyz()[2] : (*pt)[UF]
                                                                                        ? (*pt)[UF]->getXyz()[2]
                                                                                        : printErrorToDouble("zp");
    auto g{AGM::heatFunction()};
    pointHeat temp{};
    auto f = [&](const point &point) -> double {
        temp.point::operator=(point);
        return g.u(pointHeat::getTime(), temp);
    };
    double value =
            string == "grad" ? std::sqrt(std::pow((*pt)["dx"], 2) + std::pow((*pt)["dy"], 2) + std::pow((*pt)["dz"], 2))
                             : (*pt)[string];
    double numerator{std::pow(value - f(*pt), 2) * (xp - xm) * (yp - ym) * (zp - zm) * 1.25E-1};
    double denominator{std::pow(f(*pt), 2) * (xp - xm) * (yp - ym) * (zp - zm) * 1.25E-1};
    return std::make_pair(numerator, denominator);
}

template<typename T>
auto AGM::writeFile<T>::calculateError(const std::string &string) -> double {
    double numerator{};
    double denominator{};
    auto error = std::pair<double, double>{};
    for (auto &item: *pts) {
        pt = &item;
        error = calculateErrorAtPoint(string);
        numerator += error.first;
        denominator += error.second;
    }
    return std::sqrt(numerator / denominator);
}

template<typename T>
void AGM::writeFile<T>::writeResult(const std::string &string) {
    std::ofstream f(string);
    f.precision(16);
    for (const auto &item: *pts) {
        f << std::scientific;
        f << item.getIdx() << "\t";
        f << item[0] << "\t";
        f << item[1] << "\t";
        f << item[2] << "\t";
        f << item["sol"] << "\t";
        f << item["phi"] << "\t";
        f << item["psi"] << "\t";
        f << item["dx"] << "\t";
        f << item["dy"] << "\t";
        f << item["dz"] << "\n";
    }
    f.close();
}

template
class AGM::writeFile<AGM::point>;

template
class AGM::writeFile<AGM::pointHeat>;




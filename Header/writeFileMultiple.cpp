//
// Created by NIMS-JUNHONG on 2022/03/04.
//

#include "writeFileMultiple.h"
#include "StdVector"

template<typename T0, typename T1, typename T2, typename T3>
AGM::writeFileMultiple<T0, T1, T2, T3>::writeFileMultiple() = default;

template<typename T0, typename T1, typename T2, typename T3>
AGM::writeFileMultiple<T0, T1, T2, T3>::writeFileMultiple(std::vector<T0> *pts0, std::vector<T1> *pts1,
                                                          std::vector<T2> *pts2, std::vector<T3> *pts3):pts0(pts0),
                                                                                                        pts1(pts1),
                                                                                                        pts2(pts2),
                                                                                                        pts3(pts3) {}

template<typename T0, typename T1, typename T2, typename T3>
void AGM::writeFileMultiple<T0, T1, T2, T3>::writeResult(const std::string &string) {
    int bc{};
    std::ofstream f(string);
    f.precision(16);
    for (int i = 0; i < point::getNPts(); ++i) {
        bc = pts0->at(i).getCondition() == 'C' ? 0 : pts0->at(i).getCondition() == 'D' ? 1 :
                                                     pts0->at(i).getCondition() == 'N' ? 2 :
                                                     pts0->at(i).getCondition() == 'I' ? 3 : 4;
        f << std::scientific;
        f << i << "\t";                    // idx
        f << pts0->at(i)[0] << "\t";       // x
        f << pts0->at(i)[1] << "\t";       // y
        f << pts0->at(i)[2] << "\t";       // z
        f << pts0->at(i)["sol"] << "\t";   // u
        f << pts1->at(i)["sol"] << "\t";   // v
        f << pts2->at(i)["sol"] << "\t";   // w
        f << pts3->at(i)["sol"] << "\t";   // p
        f << pts0->at(i)["dx"] << "\t";    // ux
        f << pts0->at(i)["dy"] << "\t";    // uy
        f << pts0->at(i)["dz"] << "\t";    // uz
        f << pts1->at(i)["dx"] << "\t";    // vx
        f << pts1->at(i)["dy"] << "\t";    // vy
        f << pts1->at(i)["dz"] << "\t";    // vz
        f << pts2->at(i)["dx"] << "\t";    // wx
        f << pts2->at(i)["dy"] << "\t";    // wy
        f << pts2->at(i)["dz"] << "\t";    // wz
        f << pts0->at(i)["phi"] << "\t";   // phi of u
        f << pts0->at(i)["psi"] << "\t";   // psi of u
        f << pts1->at(i)["phi"] << "\t";   // phi of v
        f << pts1->at(i)["psi"] << "\t";   // psi of v
        f << pts2->at(i)["phi"] << "\t";   // phi of w
        f << pts2->at(i)["psi"] << "\t";   // psi of w
        f << bc << "\n";
    }
    f.close();
}

template<typename T0, typename T1, typename T2, typename T3>
AGM::writeFileMultiple<T0, T1, T2, T3>::~writeFileMultiple() = default;

template
class AGM::writeFileMultiple<AGM::pointHeat, AGM::pointHeat, AGM::pointHeat, AGM::point>;

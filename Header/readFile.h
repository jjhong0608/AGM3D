//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#ifndef AGM3D_READFILE_H
#define AGM3D_READFILE_H

#include "matrixNormal.h"

namespace AGM {
    class readFile {
    public:
        static void
        loadData(const std::string &filename, std::vector<AGM::point> *pts, std::vector<AGM::axialLine> *alineX,
                 std::vector<AGM::axialLine> *alineY, std::vector<AGM::axialLine> *alineZ,
                 std::array<std::vector<AGM::plane>, 2> *planeXY,
                 std::array<std::vector<AGM::plane>, 2> *planeYZ, std::array<std::vector<AGM::plane>, 2> *planeXZ);
    };
}


#endif //AGM3D_READFILE_H

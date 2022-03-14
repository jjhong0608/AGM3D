//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#include "readFile.h"

void
AGM::readFile::loadData(const std::string &filename, std::vector<AGM::point> *pts, std::vector<AGM::axialLine> *alineX,
                        std::vector<AGM::axialLine> *alineY, std::vector<AGM::axialLine> *alineZ,
                        std::array<std::vector<AGM::plane>, 2> *planeXY, std::array<std::vector<AGM::plane>, 2> *planeYZ,
                        std::array<std::vector<AGM::plane>, 2> *planeXZ) {
    std::ifstream f(filename);
    std::string line{}, tempString{};
    int nRegion{}, idx{}, index{}, tempInteger{}, nAxial{}, nPlane{};
    int nCross{}, nBound{}, nXAxial{}, nYAxial{}, nZAxial{}, nXYPlane{}, nYZPlane{}, nXZPlane{};
    coordinate normal{};
    double mp{}, bv{};
    char bc{};
    point pt{};
    std::vector<std::vector<int>> xline{}, yline{}, zline{};
    std::vector<std::vector<int>> xyXplane{}, xyYplane{}, yzYplane{}, yzZplane{}, xzXplane{}, xzZplane{};
    std::vector<int> aline{}, pln{};

    if (!f.is_open()) {
        printError("AGM::readFile::loadData", "No Axial data file: \"%s\"\nPlease check file name", filename.c_str());
    }

    std::cout << "Axial file: \"" << filename << "\" open\n";

    while (!f.eof()) {
        std::getline(f, line);
        if (line.empty()) std::getline(f, line);
        if (line.find("ENDREGION") != std::string::npos) {
            idx = 0;
            std::getline(f, line);
        }
        if (line.find("REGION") != std::string::npos) {
            ++nRegion;
            std::cout << "REGION " << nRegion << " reading...\n";
            std::getline(f, line);
            f >> tempString >> tempString >> tempString >> mp;
        }
        if (line.find('#') != std::string::npos) {
            ++idx;
            switch (idx) {
                case 1:
                    nCross = std::stoi(line.substr(line.find('=') + 2, line.size()));
                    for (int i = 0; i < nCross; ++i) {
                        f >> tempInteger >> pt[0] >> pt[1] >> pt[2];
                        pt.setMp(mp);
                        pt.setCondition('C');
                        pt.setIdx(index++);
                        pts->emplace_back(pt);
                    }
                    break;
                case 2:
                    nBound = std::stoi(line.substr(line.find('=') + 2, line.size()));
                    for (int i = 0; i < nBound; ++i) {
                        f >> tempInteger >> pt[0] >> pt[1] >> pt[2] >> bc >> bv >> normal[0] >> normal[1] >> normal[2];
                        pt.setMp(mp);
                        if (bc != 'D' && bc != 'N' && bc != 'I' && bc != 'F') {
                            printError("AGM::ReadFile::loadAxialData", "boundary condition (which is %s) is wrong", bc);
                        }
                        pt.setCondition(bc);
                        pt["bdv"] = bv;
                        pt.setNormal(normal);
                        pt.setIdx(index++);
                        pts->emplace_back(pt);
                    }
                    break;
                case 3:
                    nXAxial = std::stoi(line.substr(line.find('=') + 2, line.size()));
                    nAxial = 0;
                    while (nAxial != nXAxial) {
                        f >> tempString;
                        while (tempString != "/") {
                            aline.emplace_back(std::stoi(tempString));
                            f >> tempString;
                        }
                        xline.emplace_back(aline);
                        aline.clear();
                        ++nAxial;
                    }
                    break;
                case 4:
                    nYAxial = std::stoi(line.substr(line.find('=') + 2, line.size()));
                    nAxial = 0;
                    while (nAxial != nYAxial) {
                        f >> tempString;
                        while (tempString != "/") {
                            aline.emplace_back(std::stoi(tempString));
                            f >> tempString;
                        }
                        yline.emplace_back(aline);
                        aline.clear();
                        ++nAxial;
                    }
                    break;
                case 5:
                    nZAxial = std::stoi(line.substr(line.find('=') + 2, line.size()));
                    nAxial = 0;
                    while (nAxial != nZAxial) {
                        f >> tempString;
                        while (tempString != "/") {
                            aline.emplace_back(std::stoi(tempString));
                            f >> tempString;
                        }
                        zline.emplace_back(aline);
                        aline.clear();
                        ++nAxial;
                    }
                    break;
                case 6:
                    nXYPlane = std::stoi(line.substr(line.find('=') + 2, line.size()));
                    nPlane = 0;
                    while (nPlane != nXYPlane) {
                        f >> tempString;
                        while (tempString != "/") {
                            pln.emplace_back(std::stoi(tempString));
                            f >> tempString;
                        }
                        xyXplane.emplace_back(pln);
                        pln.clear();
                        ++nPlane;
                    }
                    break;
                case 7:
                    nXYPlane = std::stoi(line.substr(line.find('=') + 2, line.size()));
                    nPlane = 0;
                    while (nPlane != nXYPlane) {
                        f >> tempString;
                        while (tempString != "/") {
                            pln.emplace_back(std::stoi(tempString));
                            f >> tempString;
                        }
                        xyYplane.emplace_back(pln);
                        pln.clear();
                        ++nPlane;
                    }
                    break;
                case 8:
                    nYZPlane = std::stoi(line.substr(line.find('=') + 2, line.size()));
                    nPlane = 0;
                    while (nPlane != nYZPlane) {
                        f >> tempString;
                        while (tempString != "/") {
                            pln.emplace_back(std::stoi(tempString));
                            f >> tempString;
                        }
                        yzYplane.emplace_back(pln);
                        pln.clear();
                        ++nPlane;
                    }
                    break;
                case 9:
                    nYZPlane = std::stoi(line.substr(line.find('=') + 2, line.size()));
                    nPlane = 0;
                    while (nPlane != nYZPlane) {
                        f >> tempString;
                        while (tempString != "/") {
                            pln.emplace_back(std::stoi(tempString));
                            f >> tempString;
                        }
                        yzZplane.emplace_back(pln);
                        pln.clear();
                        ++nPlane;
                    }
                    break;
                case 10:
                    nXZPlane = std::stoi(line.substr(line.find('=') + 2, line.size()));
                    nPlane = 0;
                    while (nPlane != nXZPlane) {
                        f >> tempString;
                        while (tempString != "/") {
                            pln.emplace_back(std::stoi(tempString));
                            f >> tempString;
                        }
                        xzXplane.emplace_back(pln);
                        pln.clear();
                        ++nPlane;
                    }
                    break;
                case 11:
                    nXZPlane = std::stoi(line.substr(line.find('=') + 2, line.size()));
                    nPlane = 0;
                    while (nPlane != nXZPlane) {
                        f >> tempString;
                        while (tempString != "/") {
                            pln.emplace_back(std::stoi(tempString));
                            f >> tempString;
                        }
                        xzZplane.emplace_back(pln);
                        pln.clear();
                        ++nPlane;
                    }
                    break;
                default:
                    printError("AGM::ReadFile::loadAxialData", "idx (which is %s) is wrong", idx);
            }
        }
    }
    f.close();

    auto xAxialLine = axialLine('x');
    auto yAxialLine = axialLine('y');
    auto zAxialLine = axialLine('z');

    for (const auto &item: xline) {
        xAxialLine[0] = pts->at(item.front())[0];
        xAxialLine[1] = pts->at(item.back())[0];
        xAxialLine[2] = pts->at(item.front())[1];
        xAxialLine[3] = pts->at(item.back())[1];
        xAxialLine[4] = pts->at(item.front())[2];
        xAxialLine[5] = pts->at(item.back())[2];
        for (const auto &item0: item) {
            xAxialLine.emplace_back(&(pts->at(item0)));
        }
        alineX->emplace_back(xAxialLine);
        xAxialLine.clear();
    }
    for (const auto &item: yline) {
        yAxialLine[0] = pts->at(item.front())[0];
        yAxialLine[1] = pts->at(item.back())[0];
        yAxialLine[2] = pts->at(item.front())[1];
        yAxialLine[3] = pts->at(item.back())[1];
        yAxialLine[4] = pts->at(item.front())[2];
        yAxialLine[5] = pts->at(item.back())[2];
        for (const auto &item0: item) {
            yAxialLine.emplace_back(&(pts->at(item0)));
        }
        alineY->emplace_back(yAxialLine);
        yAxialLine.clear();
    }
    for (const auto &item: zline) {
        zAxialLine[0] = pts->at(item.front())[0];
        zAxialLine[1] = pts->at(item.back())[0];
        zAxialLine[2] = pts->at(item.front())[1];
        zAxialLine[3] = pts->at(item.back())[1];
        zAxialLine[4] = pts->at(item.front())[2];
        zAxialLine[5] = pts->at(item.back())[2];
        for (const auto &item0: item) {
            zAxialLine.emplace_back(&(pts->at(item0)));
        }
        alineZ->emplace_back(zAxialLine);
        zAxialLine.clear();
    }

    for (auto &item: *alineX) {
        for (const auto &item0: item) {
            item0->setAxialLine(&item, 'x');
        }
    }
    for (auto &item: *alineY) {
        for (const auto &item0: item) {
            item0->setAxialLine(&item, 'y');
        }
    }
    for (auto &item: *alineZ) {
        for (const auto &item0: item) {
            item0->setAxialLine(&item, 'z');
        }
    }

    point *prev{};

    for (auto &item: *alineX) {
        for (auto &item0: item) {
            if (item0 == item.front()) {
                if (item0->getCondition() != 'I') {
                    (*item0)[L] = item0;
                    if (ispositive(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[F] = item0;
                    if (isnegative(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[B] = item0;
                    if (iszero(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[F] = (*item0)[B] = item0;
                    if (ispositive(item0->getNormal()[2]) && !item0->getAxialLine('z')) (*item0)[U] = item0;
                    if (isnegative(item0->getNormal()[2]) && !item0->getAxialLine('z')) (*item0)[D] = item0;
                    if (iszero(item0->getNormal()[2]) && !item0->getAxialLine('z')) (*item0)[U] = (*item0)[D] = item0;
                }
                prev = item0;
            } else if (item0 == item.back()) {
                (*prev)[R] = item0;
                if (item0->getCondition() != 'I') {
                    (*item0)[R] = item0;
                    if (ispositive(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[F] = item0;
                    if (isnegative(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[B] = item0;
                    if (iszero(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[F] = (*item0)[B] = item0;
                    if (ispositive(item0->getNormal()[2]) && !item0->getAxialLine('z')) (*item0)[U] = item0;
                    if (isnegative(item0->getNormal()[2]) && !item0->getAxialLine('z')) (*item0)[D] = item0;
                    if (iszero(item0->getNormal()[2]) && !item0->getAxialLine('z')) (*item0)[U] = (*item0)[D] = item0;
                }
                (*item0)[L] = prev;
            } else {
                (*prev)[R] = item0;
                (*item0)[L] = prev;
                prev = item0;
            }
        }
    }
    for (auto &item: *alineY) {
        for (auto &item0: item) {
            if (item0 == item.front()) {
                if (item0->getCondition() != 'I') {
                    (*item0)[B] = item0;
                    if (ispositive(item0->getNormal()[0]) && !item0->getAxialLine('x')) (*item0)[R] = item0;
                    if (isnegative(item0->getNormal()[0]) && !item0->getAxialLine('x')) (*item0)[L] = item0;
                    if (iszero(item0->getNormal()[0]) && !item0->getAxialLine('x')) (*item0)[R] = (*item0)[L] = item0;
                    if (ispositive(item0->getNormal()[2]) && !item0->getAxialLine('z')) (*item0)[U] = item0;
                    if (isnegative(item0->getNormal()[2]) && !item0->getAxialLine('z')) (*item0)[D] = item0;
                    if (iszero(item0->getNormal()[2]) && !item0->getAxialLine('z')) (*item0)[U] = (*item0)[D] = item0;
                }
                prev = item0;
            } else if (item0 == item.back()) {
                (*prev)[F] = item0;
                if (item0->getCondition() != 'I') {
                    (*item0)[F] = item0;
                    if (ispositive(item0->getNormal()[0]) && !item0->getAxialLine('x')) (*item0)[R] = item0;
                    if (isnegative(item0->getNormal()[0]) && !item0->getAxialLine('x')) (*item0)[L] = item0;
                    if (iszero(item0->getNormal()[0]) && !item0->getAxialLine('x')) (*item0)[R] = (*item0)[L] = item0;
                    if (ispositive(item0->getNormal()[2]) && !item0->getAxialLine('z')) (*item0)[U] = item0;
                    if (isnegative(item0->getNormal()[2]) && !item0->getAxialLine('z')) (*item0)[D] = item0;
                    if (iszero(item0->getNormal()[2]) && !item0->getAxialLine('z')) (*item0)[U] = (*item0)[D] = item0;
                }
                (*item0)[B] = prev;
            } else {
                (*prev)[F] = item0;
                (*item0)[B] = prev;
                prev = item0;
            }
        }
    }
    for (auto &item: *alineZ) {
        for (auto &item0: item) {
            if (item0 == item.front()) {
                if (item0->getCondition() != 'I') {
                    (*item0)[D] = item0;
                    if (ispositive(item0->getNormal()[0]) && !item0->getAxialLine('x')) (*item0)[R] = item0;
                    if (isnegative(item0->getNormal()[0]) && !item0->getAxialLine('x')) (*item0)[L] = item0;
                    if (iszero(item0->getNormal()[0]) && !item0->getAxialLine('x')) (*item0)[R] = (*item0)[L] = item0;
                    if (ispositive(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[F] = item0;
                    if (isnegative(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[B] = item0;
                    if (iszero(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[F] = (*item0)[B] = item0;
                }
                prev = item0;
            } else if (item0 == item.back()) {
                (*prev)[U] = item0;
                if (item0->getCondition() != 'I') {
                    (*item0)[U] = item0;
                    if (ispositive(item0->getNormal()[0]) && !item0->getAxialLine('x')) (*item0)[R] = item0;
                    if (isnegative(item0->getNormal()[0]) && !item0->getAxialLine('x')) (*item0)[L] = item0;
                    if (iszero(item0->getNormal()[0]) && !item0->getAxialLine('x')) (*item0)[R] = (*item0)[L] = item0;
                    if (ispositive(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[F] = item0;
                    if (isnegative(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[B] = item0;
                    if (iszero(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[F] = (*item0)[B] = item0;
                }
                (*item0)[D] = prev;
            } else {
                (*prev)[U] = item0;
                (*item0)[D] = prev;
                prev = item0;
            }
        }
    }

    auto xyxPlane = plane("xy-x");
    auto xyyPlane = plane("xy-y");
    auto yzyPlane = plane("yz-y");
    auto yzzPlane = plane("yz-z");
    auto xzxPlane = plane("xz-x");
    auto xzzPlane = plane("xz-z");

    for (const auto &item: xyXplane) {
        xyxPlane.setCoordinate(alineX->at(item.front())[4]);
        for (const auto &item0: item) {
            xyxPlane.emplace_back(&(alineX->at(item0)));
        }
        planeXY->at(0).emplace_back(xyxPlane);
        xyxPlane.clear();
    }
    for (const auto &item: xyYplane) {
        xyyPlane.setCoordinate(alineY->at(item.front())[4]);
        for (const auto &item0: item) {
            xyyPlane.emplace_back(&(alineY->at(item0)));
        }
        planeXY->at(1).emplace_back(xyyPlane);
        xyyPlane.clear();
    }
    for (const auto &item: yzYplane) {
        yzyPlane.setCoordinate(alineY->at(item.front())[0]);
        for (const auto &item0: item) {
            yzyPlane.emplace_back(&(alineY->at(item0)));
        }
        planeYZ->at(0).emplace_back(yzyPlane);
        yzyPlane.clear();
    }
    for (const auto &item: yzZplane) {
        yzzPlane.setCoordinate(alineZ->at(item.front())[0]);
        for (const auto &item0: item) {
            yzzPlane.emplace_back(&(alineZ->at(item0)));
        }
        planeYZ->at(1).emplace_back(yzzPlane);
        yzzPlane.clear();
    }
    for (const auto &item: xzXplane) {
        xzxPlane.setCoordinate(alineX->at(item.front())[2]);
        for (const auto &item0: item) {
            xzxPlane.emplace_back(&(alineX->at(item0)));
        }
        planeXZ->at(0).emplace_back(xzxPlane);
        xzxPlane.clear();
    }
    for (const auto &item: xzZplane) {
        xzzPlane.setCoordinate(alineZ->at(item.front())[2]);
        for (const auto &item0: item) {
            xzzPlane.emplace_back(&(alineZ->at(item0)));
        }
        planeXZ->at(1).emplace_back(xzzPlane);
        xzzPlane.clear();
    }

    for (auto &item: planeXY->at(0)) {
        for (const auto &item0: item) {
            item0->setPlane(&item, 0);
        }
    }
    for (auto &item: planeXY->at(1)) {
        for (const auto &item0: item) {
            item0->setPlane(&item, 0);
        }
    }
    for (auto &item: planeYZ->at(0)) {
        for (const auto &item0: item) {
            item0->setPlane(&item, 1);
        }
    }
    for (auto &item: planeYZ->at(1)) {
        for (const auto &item0: item) {
            item0->setPlane(&item, 0);
        }
    }
    for (auto &item: planeXZ->at(0)) {
        for (const auto &item0: item) {
            item0->setPlane(&item, 1);
        }
    }
    for (auto &item: planeXZ->at(1)) {
        for (const auto &item0: item) {
            item0->setPlane(&item, 1);
        }
    }

    auto element = AGM::axialElement{};
    for (auto &item: *pts) {
        element = item.getElement();
        if (item.getAxialLine('x')) {
            if (element[R]) element[R] = item[R]->getElement()[R];
            if (element[L]) element[L] = item[L]->getElement()[L];
        }
        if (item.getAxialLine('y')) {
            if (element[F]) element[F] = item[F]->getElement()[F];
            if (element[B]) element[B] = item[B]->getElement()[B];
        }
        if (item.getAxialLine('z')) {
            if (element[U]) element[U] = item[U]->getElement()[U];
            if (element[D]) element[D] = item[D]->getElement()[D];
        }
        item.setElement1(element);
    }
}

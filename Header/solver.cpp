//
// Created by NIMS-JUNHONG on 2022/01/06.
//

#include "solver.h"

AGM::solver::solver(std::vector<AGM::point> *pts) : pts(pts) {}

auto AGM::solver::getPts() const -> std::vector<AGM::point> * {
    return pts;
}

void AGM::solver::setPts(std::vector<AGM::point> *vector) {
    solver::pts = vector;
}

void AGM::solver::ellipticSolver() {
    auto f{AGM::ellipticFunction()};
    auto rhsX = [&](int i) -> double {
        return ZEROVALUE;
    };
    auto rhsY = [&](int i) -> double {
        return HALFVALUE * f.f(pts->at(i));
    };
    auto rhsZ = [&](int i) -> double {
        return HALFVALUE * f.f(pts->at(i));
    };
    auto rhsXp = [&](int i) -> double {
        return ZEROVALUE;
    };
    auto rhsYp = [&](int i) -> double {
        return ZEROVALUE;
    };
    auto rhsZp = [&](int i) -> double {
        return ZEROVALUE;
    };

    #pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
        f.assignBoundaryValue(*item);
    }

    point::setNPts(int(pts->size()));
    #pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
        item->calculateRepresentationFormula();
        item->makeDerivatives();
        item->updateRightHandSide(rhsX, rhsY, rhsZ);
        item->updateRightHandSidePart(rhsXp, rhsYp, rhsZp);
    }

    auto matrix = AGM::matrix<point>(pts);
    matrix.makeMatrix();
    matrix.factorizeMatrix();
    matrix.calculateMatrix();
    matrix.releaseMatrix();

    #pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
        item->calculateDerivatives(pts, rhsX, rhsY, rhsZ, rhsXp, rhsYp, rhsZp);
    }
    #pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
        item->approximateNaNDerivatives(pts);
    }
}

void AGM::solver::heatSolver() {
    auto f{AGM::heatFunction()};
    point::setNPts(int(pts->size()));
    auto ptsHeat{std::vector<pointHeat>(point::getNPts())};
    auto previousValues{std::vector<value>(point::getNPts())};
    pointHeat::setTime(UNITVALUE);
    pointHeat::setDelta(1e-2 / 4);
    auto rhsX = [&](int i) -> double {
//        return previousValues.at(i)["sol"] / (3.0 * pointHeat::getDelta());
        return previousValues.at(i)["phi"] + previousValues.at(i)["psi"] +
               2.0 / 3.0 / pointHeat::getDelta() * previousValues.at(i)["sol"];
    };
    auto rhsY = [&](int i) -> double {
//        return HALFVALUE * (ptsHeat.at(i)["rhs"]) + previousValues.at(i)["sol"] / (3.0 * pointHeat::getDelta());
        return HALFVALUE * (ptsHeat.at(i)["rhs"] + previousValues.at(i)["rhs"]) - previousValues.at(i)["phi"] +
               2.0 / 3.0 / pointHeat::getDelta() * previousValues.at(i)["sol"];
    };
    auto rhsZ = [&](int i) -> double {
//        return HALFVALUE * (ptsHeat.at(i)["rhs"]) + previousValues.at(i)["sol"] / (3.0 * pointHeat::getDelta());
        return HALFVALUE * (ptsHeat.at(i)["rhs"] + previousValues.at(i)["rhs"]) - previousValues.at(i)["psi"] +
               2.0 / 3.0 / pointHeat::getDelta() * previousValues.at(i)["sol"];
    };
    auto rhsXp = [&](int i) -> double {
//        return ZEROVALUE;
        return -ptsHeat.at(i).getMp() * previousValues.at(i)["dx"];
    };
    auto rhsYp = [&](int i) -> double {
//        return ZEROVALUE;
        return -ptsHeat.at(i).getMp() * previousValues.at(i)["dy"];
    };
    auto rhsZp = [&](int i) -> double {
//        return ZEROVALUE;
        return -ptsHeat.at(i).getMp() * previousValues.at(i)["dz"];
    };
    auto copyPointInformation = [this, &ptsHeat]() -> void {
        for (int i = 0; i < point::getNPts(); ++i) {
            ptsHeat.at(i).point::operator=(pts->at(i));
            ptsHeat.at(i).findStencil(&(pts->at(i).getElement()), &ptsHeat);
        }
    };
    auto assignInitial = [&f, &ptsHeat, &previousValues]() -> void {
        #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            f.assignPreviousValue(previousValues.at(i), ptsHeat.at(i));
            f.assignBoundaryValue(ptsHeat.at(i));
        }
    };
    auto assignBoundaryValue = [&f, &ptsHeat]() -> void {
        #pragma omp parallel for
        for (auto item = ptsHeat.begin(); item != ptsHeat.end(); ++item) {
            f.assignBoundaryValue(*item);
        }
    };
    auto makeMatrix = [&rhsX, &rhsY, &rhsZ, &rhsXp, &rhsYp, &rhsZp, &ptsHeat]() -> void {
        #pragma omp parallel for
        for (auto item = ptsHeat.begin(); item != ptsHeat.end(); ++item) {
            item->calculateRepresentationFormula();
            item->makeDerivatives();
            item->updateRightHandSide(rhsX, rhsY, rhsZ);
            item->updateRightHandSidePart(rhsXp, rhsYp, rhsZp);
        }
    };
    auto calculateDifferentitation = [&rhsX, &rhsY, &rhsZ, &rhsXp, &rhsYp, &rhsZp, &ptsHeat]() -> void {
        #pragma omp parallel for
        for (auto item = ptsHeat.begin(); item != ptsHeat.end(); ++item) {
            item->calculateDerivatives(&ptsHeat, rhsX, rhsY, rhsZ, rhsXp, rhsYp, rhsZp);
        }
        #pragma omp parallel for
        for (auto item = ptsHeat.begin(); item != ptsHeat.end(); ++item) {
            item->approximateNaNDerivatives(&ptsHeat);
        }
    };
    auto updateRHS = [&rhsX, &rhsY, &rhsZ, &rhsXp, &rhsYp, &rhsZp, &ptsHeat]() -> void {
        #pragma omp parallel for
        for (auto item = ptsHeat.begin(); item != ptsHeat.end(); ++item) {
            item->updateRightHandSide(rhsX, rhsY, rhsZ);
            item->updateRightHandSidePart(rhsXp, rhsYp, rhsZp);
        }
    };
    auto updateValue = [&ptsHeat, &previousValues]() -> void {
        #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            previousValues.at(i) = ptsHeat.at(i).getValue();
        }
    };
    auto updateTime = []() -> void {
        pointHeat::setTime(pointHeat::getTime() + pointHeat::getDelta());
        std::cout << "current time = " << pointHeat::getTime() << "\n";
    };
    std::cout << "# of the points = " << point::getNPts() << "\n";
    copyPointInformation();
    assignInitial();
    makeMatrix();
    auto matrix = AGM::matrix<pointHeat>(&ptsHeat);
    matrix.makeMatrix();
    matrix.factorizeMatrix();
    matrix.calculateMatrix();
    calculateDifferentitation();
    updateValue();
    while (pointHeat::getTime() + pointHeat::getDelta() < 1.25 - HALFVALUE * pointHeat::getDelta()) {
        updateTime();
        assignBoundaryValue();
        updateRHS();
        matrix.calculateMatrix();
        calculateDifferentitation();
        updateValue();
    }
    updateTime();
    matrix.releaseMatrix();
    for (int i = 0; i < point::getNPts(); ++i) {
        pts->at(i).setValue(ptsHeat.at(i).getValue());
    }
}

void AGM::solver::NavierStokesSolver() {
    auto f{AGM::NavierStokesFunction()};
    int fixedPointIndex{};
    int presentIter{};
    int saveIter{};
    point::setNPts(int(pts->size()));
    auto uvel{std::vector<pointHeat>(point::getNPts())};
    auto vvel{std::vector<pointHeat>(point::getNPts())};
    auto wvel{std::vector<pointHeat>(point::getNPts())};
    auto puvel{std::vector<value>(point::getNPts())};
    auto pvvel{std::vector<value>(point::getNPts())};
    auto pwvel{std::vector<value>(point::getNPts())};
    auto ppvel{std::vector<value>(point::getNPts())};
    auto tildeuvel{std::vector<value>(point::getNPts())};
    auto tildevvel{std::vector<value>(point::getNPts())};
    auto tildewvel{std::vector<value>(point::getNPts())};
    auto hatuvel{std::vector<value>(point::getNPts())};
    auto hatvvel{std::vector<value>(point::getNPts())};
    auto hatwvel{std::vector<value>(point::getNPts())};
    pointHeat::setTime(AGM::NavierStokesFunction::initialTime());
    pointHeat::setDelta(AGM::NavierStokesFunction::deltaTime());
    auto old_uRhsX = [&](int i) -> double {
        return 2.0 / 3.0 * puvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto old_uRhsY = [&](int i) -> double {
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) +
               2.0 / 3.0 * puvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto old_uRhsZ = [&](int i) -> double {
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) +
               2.0 / 3.0 * puvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto old_vRhsX = [&](int i) -> double {
        return 2.0 / 3.0 * pvvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto old_vRhsY = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) +
               2.0 / 3.0 * pvvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto old_vRhsZ = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) +
               2.0 / 3.0 * pvvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto old_wRhsX = [&](int i) -> double {
        return 2.0 / 3.0 * pwvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto old_wRhsY = [&](int i) -> double {
        return HALFVALUE * (wvel.at(i)["rhs"] + pwvel.at(i)["rhs"]) +
               2.0 / 3.0 * pwvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto old_wRhsZ = [&](int i) -> double {
        return HALFVALUE * (wvel.at(i)["rhs"] + pwvel.at(i)["rhs"]) +
               2.0 / 3.0 * pwvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto uRhsX = [&](int i) -> double {
        return 4.0 / 3.0 * puvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto uRhsY = [&](int i) -> double {
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) +
               4.0 / 3.0 * puvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto uRhsZ = [&](int i) -> double {
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) +
               4.0 / 3.0 * puvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto vRhsX = [&](int i) -> double {
        return 4.0 / 3.0 * pvvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto vRhsY = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) +
               4.0 / 3.0 * pvvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto vRhsZ = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) +
               4.0 / 3.0 * pvvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto wRhsX = [&](int i) -> double {
        return 4.0 / 3.0 * pwvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto wRhsY = [&](int i) -> double {
        return HALFVALUE * (wvel.at(i)["rhs"] + pwvel.at(i)["rhs"]) +
               4.0 / 3.0 * pwvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto wRhsZ = [&](int i) -> double {
        return HALFVALUE * (wvel.at(i)["rhs"] + pwvel.at(i)["rhs"]) +
               4.0 / 3.0 * pwvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto pRhsX = [&](int i) -> double {
        return ZEROVALUE;
    };
    auto pRhsY = [&](int i) -> double {
        return ZEROVALUE;
    };
    auto pRhsZ = [&](int i) -> double {
        return ZEROVALUE;
    };
    auto old_uRhsX1 = [&](int i) -> double {
        return old_uRhsX(i);
    };
    auto old_uRhsY1 = [&](int i) -> double {
        return old_uRhsY(i);
    };
    auto old_uRhsZ1 = [&](int i) -> double {
        return old_uRhsZ(i);
    };
    auto old_vRhsX1 = [&](int i) -> double {
        return old_vRhsX(i);
    };
    auto old_vRhsY1 = [&](int i) -> double {
        return old_vRhsY(i);
    };
    auto old_vRhsZ1 = [&](int i) -> double {
        return old_vRhsZ(i);
    };
    auto old_wRhsX1 = [&](int i) -> double {
        return old_wRhsX(i);
    };
    auto old_wRhsY1 = [&](int i) -> double {
        return old_wRhsY(i);
    };
    auto old_wRhsZ1 = [&](int i) -> double {
        return old_wRhsZ(i);
    };
    auto uRhsX1 = [&](int i) -> double {
        return uRhsX(i);
    };
    auto uRhsY1 = [&](int i) -> double {
        return uRhsY(i);
    };
    auto uRhsZ1 = [&](int i) -> double {
        return uRhsZ(i);
    };
    auto vRhsX1 = [&](int i) -> double {
        return vRhsX(i);
    };
    auto vRhsY1 = [&](int i) -> double {
        return vRhsY(i);
    };
    auto vRhsZ1 = [&](int i) -> double {
        return vRhsZ(i);
    };
    auto wRhsX1 = [&](int i) -> double {
        return wRhsX(i);
    };
    auto wRhsY1 = [&](int i) -> double {
        return wRhsY(i);
    };
    auto wRhsZ1 = [&](int i) -> double {
        return wRhsZ(i);
    };
    auto old_uRhsXp = [&](int i) -> double {
        return 2 * puvel.at(i)["sol"] * puvel.at(i)["sol"] - uvel.at(i).getMp() * puvel.at(i)["dx"];
    };
    auto old_uRhsYp = [&](int i) -> double {
        return 2 * puvel.at(i)["sol"] * pvvel.at(i)["sol"] - uvel.at(i).getMp() * puvel.at(i)["dy"];
    };
    auto old_uRhsZp = [&](int i) -> double {
        return 2 * puvel.at(i)["sol"] * pwvel.at(i)["sol"] - uvel.at(i).getMp() * puvel.at(i)["dz"];
    };
    auto old_vRhsXp = [&](int i) -> double {
        return 2 * puvel.at(i)["sol"] * pvvel.at(i)["sol"] - vvel.at(i).getMp() * pvvel.at(i)["dx"];
    };
    auto old_vRhsYp = [&](int i) -> double {
        return 2 * pvvel.at(i)["sol"] * pvvel.at(i)["sol"] - vvel.at(i).getMp() * pvvel.at(i)["dy"];
    };
    auto old_vRhsZp = [&](int i) -> double {
        return 2 * pvvel.at(i)["sol"] * pwvel.at(i)["sol"] - vvel.at(i).getMp() * pvvel.at(i)["dz"];
    };
    auto old_wRhsXp = [&](int i) -> double {
        return 2 * puvel.at(i)["sol"] * pwvel.at(i)["sol"] - wvel.at(i).getMp() * pwvel.at(i)["dx"];
    };
    auto old_wRhsYp = [&](int i) -> double {
        return 2 * pvvel.at(i)["sol"] * pwvel.at(i)["sol"] - wvel.at(i).getMp() * pwvel.at(i)["dy"];
    };
    auto old_wRhsZp = [&](int i) -> double {
        return 2 * pwvel.at(i)["sol"] * pwvel.at(i)["sol"] - wvel.at(i).getMp() * pwvel.at(i)["dz"];
    };
    auto uRhsXp = [&](int i) -> double {
        return 2 * puvel.at(i)["sol"] * puvel.at(i)["sol"];
    };
    auto uRhsYp = [&](int i) -> double {
        return 2 * puvel.at(i)["sol"] * pvvel.at(i)["sol"];
    };
    auto uRhsZp = [&](int i) -> double {
        return 2 * puvel.at(i)["sol"] * pwvel.at(i)["sol"];
    };
    auto vRhsXp = [&](int i) -> double {
        return 2 * puvel.at(i)["sol"] * pvvel.at(i)["sol"];
    };
    auto vRhsYp = [&](int i) -> double {
        return 2 * pvvel.at(i)["sol"] * pvvel.at(i)["sol"];
    };
    auto vRhsZp = [&](int i) -> double {
        return 2 * pvvel.at(i)["sol"] * pwvel.at(i)["sol"];
    };
    auto wRhsXp = [&](int i) -> double {
        return 2 * puvel.at(i)["sol"] * pwvel.at(i)["sol"];
    };
    auto wRhsYp = [&](int i) -> double {
        return 2 * pvvel.at(i)["sol"] * pwvel.at(i)["sol"];
    };
    auto wRhsZp = [&](int i) -> double {
        return 2 * pwvel.at(i)["sol"] * pwvel.at(i)["sol"];
    };
    auto pRhsXp = [&](int i) -> double {
        return uvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto pRhsYp = [&](int i) -> double {
        return vvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto pRhsZp = [&](int i) -> double {
        return wvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto old_uRhsXp1 = [&](int i) -> double {
        return hatuvel.at(i)["sol"] * hatuvel.at(i)["sol"] + puvel.at(i)["sol"] * puvel.at(i)["sol"] + pts->at(i)["sol"] +
               ppvel.at(i)["sol"] - uvel.at(i).getMp() * puvel.at(i)["dx"];
    };
    auto old_uRhsYp1 = [&](int i) -> double {
        return hatuvel.at(i)["sol"] * hatvvel.at(i)["sol"] + puvel.at(i)["sol"] * pvvel.at(i)["sol"] -
               uvel.at(i).getMp() * puvel.at(i)["dy"];
    };
    auto old_uRhsZp1 = [&](int i) -> double {
        return hatuvel.at(i)["sol"] * hatwvel.at(i)["sol"] + puvel.at(i)["sol"] * pwvel.at(i)["sol"] -
               uvel.at(i).getMp() * puvel.at(i)["dz"];
    };
    auto old_vRhsXp1 = [&](int i) -> double {
        return hatuvel.at(i)["sol"] * hatvvel.at(i)["sol"] + puvel.at(i)["sol"] * pvvel.at(i)["sol"] -
               vvel.at(i).getMp() * pvvel.at(i)["dx"];
    };
    auto old_vRhsYp1 = [&](int i) -> double {
        return hatvvel.at(i)["sol"] * hatvvel.at(i)["sol"] + pvvel.at(i)["sol"] * pvvel.at(i)["sol"] + pts->at(i)["sol"] +
               ppvel.at(i)["sol"] - vvel.at(i).getMp() * pvvel.at(i)["dy"];
    };
    auto old_vRhsZp1 = [&](int i) -> double {
        return hatvvel.at(i)["sol"] * hatwvel.at(i)["sol"] + pvvel.at(i)["sol"] * pwvel.at(i)["sol"] -
               vvel.at(i).getMp() * pvvel.at(i)["dz"];
    };
    auto old_wRhsXp1 = [&](int i) -> double {
        return hatuvel.at(i)["sol"] * hatwvel.at(i)["sol"] + puvel.at(i)["sol"] * pwvel.at(i)["sol"] -
               wvel.at(i).getMp() * pwvel.at(i)["dx"];
    };
    auto old_wRhsYp1 = [&](int i) -> double {
        return hatvvel.at(i)["sol"] * hatwvel.at(i)["sol"] + pvvel.at(i)["sol"] * pwvel.at(i)["sol"] -
               wvel.at(i).getMp() * pwvel.at(i)["dy"];
    };
    auto old_wRhsZp1 = [&](int i) -> double {
        return hatwvel.at(i)["sol"] * hatwvel.at(i)["sol"] + pwvel.at(i)["sol"] * pwvel.at(i)["sol"] + pts->at(i)["sol"] +
               ppvel.at(i)["sol"] - wvel.at(i).getMp() * pwvel.at(i)["dz"];
    };
    auto uRhsXp1 = [&](int i) -> double {
        return hatuvel.at(i)["sol"] * hatuvel.at(i)["sol"] + puvel.at(i)["sol"] * puvel.at(i)["sol"] + pts->at(i)["sol"] +
               ppvel.at(i)["sol"];
    };
    auto uRhsYp1 = [&](int i) -> double {
        return hatuvel.at(i)["sol"] * hatvvel.at(i)["sol"] + puvel.at(i)["sol"] * pvvel.at(i)["sol"];
    };
    auto uRhsZp1 = [&](int i) -> double {
        return hatuvel.at(i)["sol"] * hatwvel.at(i)["sol"] + puvel.at(i)["sol"] * pwvel.at(i)["sol"];
    };
    auto vRhsXp1 = [&](int i) -> double {
        return hatuvel.at(i)["sol"] * hatvvel.at(i)["sol"] + puvel.at(i)["sol"] * pvvel.at(i)["sol"];
    };
    auto vRhsYp1 = [&](int i) -> double {
        return hatvvel.at(i)["sol"] * hatvvel.at(i)["sol"] + pvvel.at(i)["sol"] * pvvel.at(i)["sol"] + pts->at(i)["sol"] +
               ppvel.at(i)["sol"];
    };
    auto vRhsZp1 = [&](int i) -> double {
        return hatvvel.at(i)["sol"] * hatwvel.at(i)["sol"] + pvvel.at(i)["sol"] * pwvel.at(i)["sol"];
    };
    auto wRhsXp1 = [&](int i) -> double {
        return hatuvel.at(i)["sol"] * hatwvel.at(i)["sol"] + puvel.at(i)["sol"] * pwvel.at(i)["sol"];
    };
    auto wRhsYp1 = [&](int i) -> double {
        return hatvvel.at(i)["sol"] * hatwvel.at(i)["sol"] + pvvel.at(i)["sol"] * pwvel.at(i)["sol"];
    };
    auto wRhsZp1 = [&](int i) -> double {
        return hatwvel.at(i)["sol"] * hatwvel.at(i)["sol"] + pwvel.at(i)["sol"] * pwvel.at(i)["sol"] + pts->at(i)["sol"] +
               ppvel.at(i)["sol"];
    };
    auto pRhsX1 = [&](int i) -> double {
        return -uvel.at(i)["dx"] / pointHeat::getDelta();
    };
    auto pRhsY1 = [&](int i) -> double {
        return -vvel.at(i)["dy"] / pointHeat::getDelta();
    };
    auto pRhsZ1 = [&](int i) -> double {
        return -wvel.at(i)["dz"] / pointHeat::getDelta();
    };
    auto findSaveIter = [&presentIter, &saveIter]() -> void {
        saveIter = int(std::floor(NavierStokesFunction::writeTime() / NavierStokesFunction::deltaTime() + HALFVALUE));
        presentIter = int(std::floor(
                (std::fmod(NavierStokesFunction::initialTime(), NavierStokesFunction::writeTime())) /
                NavierStokesFunction::deltaTime() + HALFVALUE));
        std::cout << "Initial iteration number = " << presentIter << ",\n" << "file will be saved every " << saveIter
                  << " iterations\n";
    };
    auto copyPointInformation = [this, &uvel, &vvel, &wvel, &fixedPointIndex]() -> void {
//        for (int i = 0; i < point::getNPts(); ++i) {
//            if (pts->at(i).getCondition() == 'D') {
//                if (isclose(pts->at(i).getNormal()[0], UNITVALUE) && (pts->at(i)[1] < 4e0)) {
//                    pts->at(i).setCondition('N');
//                }
//
//                if (pts->at(i).getNormal()[0] > 1e-5 && iszero(pts->at(i).getNormal()[2], 1e-3) && isclose(pts->at(i).getNormal()[0], pts->at(i).getNormal()[1], 1e-3)) {
//                    pts->at(i).setCondition('N');
//                }
//            }
//        }

        for (int i = 0; i < point::getNPts(); ++i) {
            uvel.at(i).point::operator=(pts->at(i));
            uvel.at(i).findStencil(&(pts->at(i).getElement()), &uvel);
            vvel.at(i).point::operator=(pts->at(i));
            vvel.at(i).findStencil(&(pts->at(i).getElement()), &vvel);
            wvel.at(i).point::operator=(pts->at(i));
            wvel.at(i).findStencil(&(pts->at(i).getElement()), &wvel);
            if (isclose(pts->at(i)[0], 0.5) && isclose(pts->at(i)[1], 0.5) && isclose(pts->at(i)[2], 0.5)) {
                fixedPointIndex = i;
                std::cout << "fixed point index = " << fixedPointIndex << "\n";
            }
        }
    };
    auto assignInitial = [this, &f, &uvel, &vvel, &wvel, &puvel, &pvvel, &pwvel, &ppvel]() -> void {
        #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            pts->at(i).setMp(UNITVALUE);
            if (pts->at(i).getCondition() == 'D') {
                pts->at(i).setCondition('N');
                pts->at(i)["bdv"] = ZEROVALUE;
            } else if (pts->at(i).getCondition() == 'N') {
                pts->at(i).setCondition('D');
                pts->at(i)["bdv"] = ZEROVALUE;
            }
            f.assignPreviousValue(puvel.at(i), pvvel.at(i), pwvel.at(i), ppvel.at(i), uvel.at(i), vvel.at(i),
                                  wvel.at(i), pts->at(i));
            f.assignBoundaryValue(uvel.at(i), vvel.at(i), wvel.at(i));
            uvel.at(i)["bdv"] *= HALFVALUE;
            vvel.at(i)["bdv"] *= HALFVALUE;
            wvel.at(i)["bdv"] *= HALFVALUE;
        }
//        f.loadPreviousValue(
//                "/home/jjhong0608/docker/AGM3D/Navier-Stokes/Blood_flow/10/AGM_Result_3.3", &puvel,
//                &pvvel, &pwvel, &ppvel);
    };
    auto assignBoundaryValue = [&f, &uvel, &vvel, &wvel]() -> void {
        #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            f.assignBoundaryValue(uvel.at(i), vvel.at(i), wvel.at(i));
        }
    };
    auto makeMatrixVelocity = [&]() -> void {
        #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            uvel.at(i).calculateRepresentationFormula();
            uvel.at(i).makeDerivatives();
            uvel.at(i).updateRightHandSide(old_uRhsX, old_uRhsY, old_uRhsZ);
            uvel.at(i).updateRightHandSidePart(old_uRhsXp, old_uRhsYp, old_uRhsZp);
            vvel.at(i).calculateRepresentationFormula();
            vvel.at(i).makeDerivatives();
            vvel.at(i).updateRightHandSide(old_vRhsX, old_vRhsY, old_vRhsZ);
            vvel.at(i).updateRightHandSidePart(old_vRhsXp, old_vRhsYp, old_vRhsZp);
            wvel.at(i).calculateRepresentationFormula();
            wvel.at(i).makeDerivatives();
            wvel.at(i).updateRightHandSide(old_wRhsX, old_wRhsY, old_wRhsZ);
            wvel.at(i).updateRightHandSidePart(old_wRhsXp, old_wRhsYp, old_wRhsZp);
        }
    };
    auto makeMatrixPressure = [&]() -> void {
        #pragma omp parallel for
        for (auto item = pts->begin(); item != pts->end(); ++item) {
            item->calculateRepresentationFormula();
            item->makeDerivatives();
            item->updateRightHandSide(pRhsX, pRhsY, pRhsZ);
            item->updateRightHandSidePart(pRhsXp, pRhsYp, pRhsZp);
        }
    };
    auto calculateDifferentiationVelocityOld = [&]() -> void {
        #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            uvel.at(i).calculateDerivatives(&uvel, old_uRhsX, old_uRhsY, old_uRhsZ, old_uRhsXp, old_uRhsYp, old_uRhsZp);
            vvel.at(i).calculateDerivatives(&vvel, old_vRhsX, old_vRhsY, old_vRhsZ, old_vRhsXp, old_vRhsYp, old_vRhsZp);
            wvel.at(i).calculateDerivatives(&wvel, old_wRhsX, old_wRhsY, old_wRhsZ, old_wRhsXp, old_wRhsYp, old_wRhsZp);
        }
        #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            uvel.at(i).approximateNaNDerivatives(&uvel);
            vvel.at(i).approximateNaNDerivatives(&vvel);
            wvel.at(i).approximateNaNDerivatives(&wvel);
        }
    };
    auto calculateDifferentiationVelocity = [&]() -> void {
        #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            uvel.at(i).calculateDerivatives(&uvel, uRhsX, uRhsY, uRhsZ, uRhsXp, uRhsYp, uRhsZp);
            vvel.at(i).calculateDerivatives(&vvel, vRhsX, vRhsY, vRhsZ, vRhsXp, vRhsYp, vRhsZp);
            wvel.at(i).calculateDerivatives(&wvel, wRhsX, wRhsY, wRhsZ, wRhsXp, wRhsYp, wRhsZp);
        }
        #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            uvel.at(i).approximateNaNDerivatives(&uvel);
            vvel.at(i).approximateNaNDerivatives(&vvel);
            wvel.at(i).approximateNaNDerivatives(&wvel);
        }
    };
    auto calculateDifferentiationVelocityOld1 = [&]() -> void {
        auto dx{ZEROVALUE}, dy{ZEROVALUE}, dz{ZEROVALUE};
        #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            uvel.at(i).calculateDerivatives(&uvel, old_uRhsX1, old_uRhsY1, old_uRhsZ1, old_uRhsXp1, old_uRhsYp1, old_uRhsZp1);
            vvel.at(i).calculateDerivatives(&vvel, old_vRhsX1, old_vRhsY1, old_vRhsZ1, old_vRhsXp1, old_vRhsYp1, old_vRhsZp1);
            wvel.at(i).calculateDerivatives(&wvel, old_wRhsX1, old_wRhsY1, old_wRhsZ1, old_wRhsXp1, old_wRhsYp1, old_wRhsZp1);
        }
        #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            uvel.at(i).approximateNaNDerivatives(&uvel);
            vvel.at(i).approximateNaNDerivatives(&vvel);
            wvel.at(i).approximateNaNDerivatives(&wvel);
        }
    };
    auto calculateDifferentiationVelocity1 = [&]() -> void {
        auto dx{ZEROVALUE}, dy{ZEROVALUE}, dz{ZEROVALUE};
        #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            uvel.at(i).calculateDerivatives(&uvel, uRhsX1, uRhsY1, uRhsZ1, uRhsXp1, uRhsYp1, uRhsZp1);
            vvel.at(i).calculateDerivatives(&vvel, vRhsX1, vRhsY1, vRhsZ1, vRhsXp1, vRhsYp1, vRhsZp1);
            wvel.at(i).calculateDerivatives(&wvel, wRhsX1, wRhsY1, wRhsZ1, wRhsXp1, wRhsYp1, wRhsZp1);
        }
        #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            uvel.at(i).approximateNaNDerivatives(&uvel);
            vvel.at(i).approximateNaNDerivatives(&vvel);
            wvel.at(i).approximateNaNDerivatives(&wvel);
        }
    };
    auto calculateDifferentiationPressure = [&]() -> void {
        #pragma omp parallel for
        for (auto item = pts->begin(); item != pts->end(); ++item) {
            item->calculateDerivatives(pts, pRhsX, pRhsY, pRhsZ, pRhsXp, pRhsYp, pRhsZp);
        }
        #pragma omp parallel for
        for (auto item = pts->begin(); item != pts->end(); ++item) {
            item->approximateNaNDerivatives(pts);
        }
        #pragma omp parallel for
        for (auto item = pts->begin(); item != pts->end(); ++item) {
            item->calculateSecondDerivatives(pRhsX1, pRhsY1, pRhsZ1);
        }
    };
    auto updateRHSVelocity = [&]() -> void {
        #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            uvel.at(i).updateRightHandSide(uRhsX, uRhsY, uRhsZ);
            uvel.at(i).updateRightHandSidePart(uRhsXp, uRhsYp, uRhsZp);
            vvel.at(i).updateRightHandSide(vRhsX, vRhsY, vRhsZ);
            vvel.at(i).updateRightHandSidePart(vRhsXp, vRhsYp, vRhsZp);
            wvel.at(i).updateRightHandSide(wRhsX, wRhsY, wRhsZ);
            wvel.at(i).updateRightHandSidePart(wRhsXp, wRhsYp, wRhsZp);
        }
    };
    auto updateRHSVelocityOld1 = [&]() -> void {
        #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            uvel.at(i).updateRightHandSide(old_uRhsX1, old_uRhsY1, old_uRhsZ1);
            uvel.at(i).updateRightHandSidePart(old_uRhsXp1, old_uRhsYp1, old_uRhsZp1);
            vvel.at(i).updateRightHandSide(old_vRhsX1, old_vRhsY1, old_vRhsZ1);
            vvel.at(i).updateRightHandSidePart(old_vRhsXp1, old_vRhsYp1, old_vRhsZp1);
            wvel.at(i).updateRightHandSide(old_wRhsX1, old_wRhsY1, old_wRhsZ1);
            wvel.at(i).updateRightHandSidePart(old_wRhsXp1, old_wRhsYp1, old_wRhsZp1);
        }
    };
    auto updateRHSVelocity1 = [&]() -> void {
        #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            uvel.at(i).updateRightHandSide(uRhsX1, uRhsY1, uRhsZ1);
            uvel.at(i).updateRightHandSidePart(uRhsXp1, uRhsYp1, uRhsZp1);
            vvel.at(i).updateRightHandSide(vRhsX1, vRhsY1, vRhsZ1);
            vvel.at(i).updateRightHandSidePart(vRhsXp1, vRhsYp1, vRhsZp1);
            wvel.at(i).updateRightHandSide(wRhsX1, wRhsY1, wRhsZ1);
            wvel.at(i).updateRightHandSidePart(wRhsXp1, wRhsYp1, wRhsZp1);
        }
    };
    auto updateRHSPressrue = [&]() -> void {
        #pragma omp parallel for
        for (auto item = pts->begin(); item != pts->end(); ++item) {
            item->updateRightHandSide(pRhsX, pRhsY, pRhsZ);
            item->updateRightHandSidePart(pRhsXp, pRhsYp, pRhsZp);
        }
    };
    auto updateValues = [this, &uvel, &vvel, &wvel, &hatuvel, &hatvvel, &hatwvel, &ppvel]() -> void {
        #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            if (pts->at(i).getCondition() != 'N') {
                uvel.at(i)["sol"] -= pts->at(i)["dx"] * pointHeat::getDelta();
                vvel.at(i)["sol"] -= pts->at(i)["dy"] * pointHeat::getDelta();
                wvel.at(i)["sol"] -= pts->at(i)["dz"] * pointHeat::getDelta();
            }
//            pts->at(i)["sol"] -=
//                    HALFVALUE * uvel.at(i).getMp() * (uvel.at(i)["dx"] + vvel.at(i)["dy"] + wvel.at(i)["dz"]);
            pts->at(i)["sol"] = 2 * pts->at(i)["sol"] -
                                uvel.at(i).getMp() * (uvel.at(i)["dx"] + vvel.at(i)["dy"] + wvel.at(i)["dz"]) -
                                ppvel.at(i)["sol"];
            uvel.at(i)["dx"] -= pts->at(i)["dxx"] * pointHeat::getDelta();
            vvel.at(i)["dy"] -= pts->at(i)["dyy"] * pointHeat::getDelta();
            wvel.at(i)["dz"] -= pts->at(i)["dzz"] * pointHeat::getDelta();
            hatuvel.at(i) = uvel.at(i).getValue();
            hatvvel.at(i) = vvel.at(i).getValue();
            hatwvel.at(i) = wvel.at(i).getValue();
        }
    };
    auto updateValues1 = [this, &uvel, &vvel, &wvel, &puvel, &pvvel, &pwvel, &ppvel]() -> void {
        #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            puvel.at(i) = uvel.at(i).getValue();
            pvvel.at(i) = vvel.at(i).getValue();
            pwvel.at(i) = wvel.at(i).getValue();
            ppvel.at(i) = pts->at(i).getValue();
        }
    };
    auto subtractPreviousVelocity = [&uvel, &vvel, &wvel, &puvel, &pvvel, &pwvel]() -> void {
    #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            uvel.at(i)["sol"] -= puvel.at(i)["sol"];
            vvel.at(i)["sol"] -= pvvel.at(i)["sol"];
            wvel.at(i)["sol"] -= pwvel.at(i)["sol"];
            uvel.at(i)["dx"] -= puvel.at(i)["dx"];
            vvel.at(i)["dx"] -= pvvel.at(i)["dx"];
            wvel.at(i)["dx"] -= pwvel.at(i)["dx"];
            uvel.at(i)["dy"] -= puvel.at(i)["dy"];
            vvel.at(i)["dy"] -= pvvel.at(i)["dy"];
            wvel.at(i)["dy"] -= pwvel.at(i)["dy"];
            uvel.at(i)["dz"] -= puvel.at(i)["dz"];
            vvel.at(i)["dz"] -= pvvel.at(i)["dz"];
            wvel.at(i)["dz"] -= pwvel.at(i)["dz"];
        }
    };
    auto stopCondition{ZEROVALUE};
    auto checkStopCondition = [&presentIter, &saveIter, &stopCondition, &uvel, &vvel, &wvel, &puvel, &pvvel, &pwvel]() -> bool {
        if (presentIter % saveIter != 0) {
            return false;
        }
        auto max_vel{ZEROVALUE}, max_err{ZEROVALUE};
        auto vel{ZEROVALUE}, pvel{ZEROVALUE};
        auto tolerance{1e-6};
        for (int i = 0; i < point::getNPts(); ++i) {
            vel = std::sqrt(uvel.at(i)["sol"] * uvel.at(i)["sol"] + vvel.at(i)["sol"] * vvel.at(i)["sol"] + wvel.at(i)["sol"] * wvel.at(i)["sol"]);
            pvel = std::sqrt(puvel.at(i)["sol"] * puvel.at(i)["sol"] + pvvel.at(i)["sol"] * pvvel.at(i)["sol"] + pwvel.at(i)["sol"] * pwvel.at(i)["sol"]);
            if (max_vel < vel) {
                max_vel = vel;
            }
            if (max_err < std::abs(vel - pvel)) {
                max_err = std::abs(vel - pvel);
            }
        }
        stopCondition = max_err / max_vel;
        std::cout << "Stop Condition: [" << stopCondition << " / " << tolerance << "]\n";
        if (stopCondition < tolerance) {
            return true;
        }
        return false;
    };
    auto wf{writeFileMultiple<pointHeat, pointHeat, pointHeat, point>(&uvel, &vvel, &wvel, pts)};
    auto updateTime = [&presentIter, &saveIter, &wf, &stopCondition]() -> void {
        ++presentIter;
        pointHeat::setTime(pointHeat::getTime() + pointHeat::getDelta());
        std::cout << presentIter << "-th iteration, " << "current time = [" << pointHeat::getTime() << " / "
                  << AGM::NavierStokesFunction::terminalTime() << "], Stopping Error = [" << stopCondition << "]\n";
        if (presentIter % saveIter == 0) {
            wf.writeResult(
                    "/home/jjhong0608/docker/Navier-Stokes_Result/3D/1.Lid-driven_cavity_flow/Re_1000/AGM_Result_"
                    + std::to_string(pointHeat::getTime()));
        }
    };
    findSaveIter();
    copyPointInformation();
    assignInitial();
    makeMatrixVelocity();
    auto matrixVelocity{AGM::matrixMulti<pointHeat>(&uvel, &vvel, &wvel)};
//    auto matrixPressure{AGM::matrix<point>(pts)};
    auto matrixPressure{AGM::matrixNormal<point>(pts, fixedPointIndex)};
    matrixVelocity.makeMatrix();
    matrixVelocity.factorizeMatrix();
    matrixVelocity.calculateMatrix();
    calculateDifferentiationVelocity();
//    subtractPreviousVelocity();
    makeMatrixPressure();
    matrixPressure.makeMatrix();
    matrixPressure.factorizeMatrix();
    matrixPressure.calculateMatrix();
    calculateDifferentiationPressure();
    updateValues();
    updateRHSVelocityOld1();
    matrixVelocity.calculateMatrix();
    calculateDifferentiationVelocity1();
//    subtractPreviousVelocity();
    updateValues1();
    while (pointHeat::getTime() + pointHeat::getDelta() <
           AGM::NavierStokesFunction::terminalTime() - HALFVALUE * pointHeat::getDelta()) {
        updateTime();
        assignBoundaryValue();
        updateRHSVelocity();
        matrixVelocity.calculateMatrix();
        calculateDifferentiationVelocity();
        subtractPreviousVelocity();
        updateRHSPressrue();
        matrixPressure.calculateMatrix();
        calculateDifferentiationPressure();
        updateValues();
        updateRHSVelocity1();
        matrixVelocity.calculateMatrix();
        calculateDifferentiationVelocity1();
        subtractPreviousVelocity();
        if (checkStopCondition()) {
            break;
        }
        updateValues1();
    }
    updateTime();
    matrixVelocity.releaseMatrix();
    matrixPressure.releaseMatrix();
//    wf.writeResult("/home/jjhong0608/docker/AGM3D/Navier-Stokes/Lid-driven_cavity_flow/Re1000-2/ver0");
}

AGM::solver::~solver() = default;

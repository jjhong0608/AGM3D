#include "Header/solver.h"

int main() {
    struct rlimit rlim{};
    getrlimit(RLIMIT_STACK, &rlim);
    rlim.rlim_cur = -1;
    rlim.rlim_max = -1;
    setrlimit(RLIMIT_STACK, &rlim);

    mkl_set_dynamic(0);
    mkl_set_num_threads(mkl_get_max_threads());

    auto pts = std::vector<AGM::point>{};
    auto xline = std::vector<AGM::axialLine>{};
    auto yline = std::vector<AGM::axialLine>{};
    auto zline = std::vector<AGM::axialLine>{};
    auto xyplane = std::array<std::vector<AGM::plane>, 2>{};
    auto yzplane = std::array<std::vector<AGM::plane>, 2>{};
    auto xzplane = std::array<std::vector<AGM::plane>, 2>{};
    AGM::readFile::loadData("ALG_output", &pts, &xline, &yline, &zline, &xyplane, &yzplane, &xzplane);
    AGM::point::setAxialLines(&xline, 'x');
    AGM::point::setAxialLines(&yline, 'y');
    AGM::point::setAxialLines(&zline, 'z');
    AGM::point::setPlane(&xyplane, "xy");
    AGM::point::setPlane(&yzplane, "yz");
    AGM::point::setPlane(&xzplane, "xz");

    for (auto &item: pts) {
        item.findStencil();
    }

    std::cout << "-----< information >-----" << "\n";
    std::cout << "# of the points = " << pts.size() << "\n";
    std::cout << "# of the x-axial lines = " << xline.size() << "\n";
    std::cout << "# of the y-axial lines = " << yline.size() << "\n";
    std::cout << "# of the z-axial lines = " << zline.size() << "\n";
    std::cout << "epsilon (Reynolds number) = " << pts.at(0).getMp() << "\n";
    std::cout << "-------------------------" << "\n";

    auto solver = AGM::solver(&pts);
    solver.NavierStokesSolver();

    return 0;
}

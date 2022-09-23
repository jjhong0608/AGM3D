//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#include "GreenfunctionReactionDiffusion.h"

AGM::GreenfunctionReactionDiffusion::GreenfunctionReactionDiffusion(double tm, double tau, double tp, double mpl,
                                                                    double mpr, double c) : Greenfunction(tm,
                                                                                                          tau,
                                                                                                          tp,
                                                                                                          mpl,
                                                                                                          mpr),
                                                                                            c(c),
                                                                                            alpha(std::sqrt(c / mpl)) {

}

auto AGM::GreenfunctionReactionDiffusion::E(double a, double b) const -> double {
    return UNITVALUE - std::exp(-two * alpha * (a - b));
}

auto AGM::GreenfunctionReactionDiffusion::F(double a, double b) const -> double {
    return UNITVALUE + std::exp(-two * alpha * (a - b));
}

auto AGM::GreenfunctionReactionDiffusion::l2p(double s) const -> double {
    return std::pow(s, 2) / alpha + two * s / std::pow(alpha, 2) + two / std::pow(alpha, 3);
}

auto AGM::GreenfunctionReactionDiffusion::l2m(double s) const -> double {
    return std::pow(s, 2) / alpha - two * s / std::pow(alpha, 2) + two / std::pow(alpha, 3);
}

auto AGM::GreenfunctionReactionDiffusion::l1p(double s) const -> double {
    return s / alpha + 1 / std::pow(alpha, 2);
}

auto AGM::GreenfunctionReactionDiffusion::l1m(double s) const -> double {
    return s / alpha - 1 / std::pow(alpha, 2);
}

auto AGM::GreenfunctionReactionDiffusion::l0(double a, double b) const -> double {
    return (UNITVALUE - std::exp(-alpha * (a - b))) / alpha;
}

auto AGM::GreenfunctionReactionDiffusion::green_function(double t) const -> double {
    if (t < tau || isclose(tm, t)) {
        return alpha / two / c * E(tp, tau) * E(t, tm) / E(tp, tm) * std::exp(-alpha * (tau - t));
    }
    return alpha / two / c * E(tp, t) * E(tau, tm) / E(tp, tm) * std::exp(-alpha * (t - tau));

}

auto AGM::GreenfunctionReactionDiffusion::
green_function_t(double t) const -> double {
    if (t < tau || isclose(tm, t)) {
        return std::pow(alpha, 2) / two / c * E(tp, tau) / E(tp, tm) *
               (std::exp(-alpha * (tau - t)) + std::exp(-alpha * (t + tau - two * tm)));
    }
    return -std::pow(alpha, 2) / two / c * E(tau, tm) / E(tp, tm) *
           (std::exp(-alpha * (t - tau)) + std::exp(-alpha * (two * tp - t - tau)));

}

auto AGM::GreenfunctionReactionDiffusion::green_function_tau(double t) const -> double {
    if (t < tau || isclose(tm, t)) {
        return -std::pow(alpha, 2) / two / c * E(t, tm) / E(tp, tm) *
               (std::exp(-alpha * (tau - t)) + std::exp(-alpha * (two * tp - t - tau)));
    }
    return std::pow(alpha, 2) / two / c * E(tp, t) / E(tp, tm) *
           (std::exp(-alpha * (t - tau)) + std::exp(-alpha * (t + tau - two * tm)));

}

auto AGM::GreenfunctionReactionDiffusion::green_function_ttau(double t) const -> double {
    if (t < tau || isclose(tm, t)) {
        return -std::pow(alpha, 3) / two / c / E(tp, tm) * (E(t, tm) * E(tp, tau) * std::exp(-alpha * (tau - t)) +
                                                            two * std::exp(-alpha * (two * tp - t - tau)) +
                                                            two * std::exp(-alpha * (t + tau - two * tm)));
    }
    return -std::pow(alpha, 3) / two / c / E(tp, tm) * (E(tau, tm) * E(tp, t) * std::exp(-alpha * (t - tau)) +
                                                        two * std::exp(-alpha * (two * tp - t - tau)) +
                                                        two * std::exp(-alpha * (t + tau - two * tm)));

}

auto AGM::GreenfunctionReactionDiffusion::green_function_ND(double t) const -> double {
    if (t < tau || isclose(tm, t)) {
        return alpha / two / c * E(tp, tau) * F(t, tm) / F(tp, tm) * std::exp(-alpha * (tau - t));
    }
    return alpha / two / c * E(tp, t) * F(tau, tm) / F(tp, tm) * std::exp(-alpha * (t - tau));

}

auto AGM::GreenfunctionReactionDiffusion::green_function_t_ND(double t) const -> double {
    if (t < tau || isclose(tm, t)) {
        return std::pow(alpha, 2) / two / c * E(tp, tau) / F(tp, tm) *
               (std::exp(-alpha * (tau - t)) - std::exp(-alpha * (t + tau - two * tm)));
    }
    return -std::pow(alpha, 2) / two / c * F(tau, tm) / F(tp, tm) *
           (std::exp(-alpha * (t - tau)) + std::exp(-alpha * (two * tp - t - tau)));

}

auto AGM::GreenfunctionReactionDiffusion::green_function_tau_ND(double t) const -> double {
    if (t < tau || isclose(tm, t)) {
        return -std::pow(alpha, 2) / two / c * F(t, tm) / F(tp, tm) *
               (std::exp(-alpha * (tau - t)) + std::exp(-alpha * (two * tp - t - tau)));
    }
    return std::pow(alpha, 2) / two / c * E(tp, t) / F(tp, tm) *
           (std::exp(-alpha * (t - tau)) - std::exp(-alpha * (t + tau - two * tm)));


}

auto AGM::GreenfunctionReactionDiffusion::green_function_ttau_ND(double t) const -> double {
    if (t < tau || isclose(tm, t)) {
        return -std::pow(alpha, 3) / two / c / F(tp, tm) * (F(t, tm) * E(tp, tau) * std::exp(-alpha * (tau - t)) +
                                                            two * std::exp(-alpha * (two * tp - t - tau)) -
                                                            two * std::exp(-alpha * (t + tau - two * tm)));
    }
    return -std::pow(alpha, 3) / two / c / E(tp, tm) * (F(tau, tm) * E(tp, t) * std::exp(-alpha * (t - tau)) +
                                                        two * std::exp(-alpha * (two * tp - t - tau)) -
                                                        two * std::exp(-alpha * (t + tau - two * tm)));

}

auto AGM::GreenfunctionReactionDiffusion::green_function_DN(double t) const -> double {
    if (t < tau || isclose(tm, t)) {
        return alpha / two / c * F(tp, tau) * E(t, tm) / F(tp, tm) * std::exp(-alpha * (tau - t));
    }
    return alpha / two / c * F(tp, t) * E(tau, tm) / F(tp, tm) * std::exp(-alpha * (t - tau));

}

auto AGM::GreenfunctionReactionDiffusion::green_function_t_DN(double t) const -> double {
    if (t < tau || isclose(tm, t)) {
        return std::pow(alpha, 2) / two / c * F(tp, tau) / F(tp, tm) *
               (std::exp(-alpha * (tau - t)) + std::exp(-alpha * (t + tau - two * tm)));
    }
    return -std::pow(alpha, 2) / two / c * E(tau, tm) / F(tp, tm) *
           (std::exp(-alpha * (t - tau)) - std::exp(-alpha * (two * tp - t - tau)));

}

auto AGM::GreenfunctionReactionDiffusion::green_function_tau_DN(double t) const -> double {
    if (t < tau || isclose(tm, t)) {
        return -std::pow(alpha, 2) / two / c * E(t, tm) / F(tp, tm) *
               (std::exp(-alpha * (tau - t)) - std::exp(-alpha * (two * tp - t - tau)));
    }
    return std::pow(alpha, 2) / two / c * F(tp, t) / F(tp, tm) *
           (std::exp(-alpha * (t - tau)) + std::exp(-alpha * (t + tau - two * tm)));


}

auto AGM::GreenfunctionReactionDiffusion::green_function_ttau_DN(double t) const -> double {
    if (t < tau || isclose(tm, t)) {
        return -std::pow(alpha, 3) / two / c / F(tp, tm) * (E(t, tm) * F(tp, tau) * std::exp(-alpha * (tau - t)) -
                                                            two * std::exp(-alpha * (two * tp - t - tau)) +
                                                            two * std::exp(-alpha * (t + tau - two * tm)));
    }
    return -std::pow(alpha, 3) / two / c / F(tp, tm) * (E(tau, tm) * F(tp, t) * std::exp(-alpha * (t - tau)) -
                                                        two * std::exp(-alpha * (two * tp - t - tau)) +
                                                        two * std::exp(-alpha * (t + tau - two * tm)));

}

auto AGM::GreenfunctionReactionDiffusion::integrate_square(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l2m(tau) - (l2m(tm) + l2p(tm)) * std::exp(-alpha * (tau - tm)) +
                    l2p(tau) * std::exp(-two * alpha * (tau - tm));
            return alpha / two / c * E(tp, tau) / E(tp, tm) * value;
        case 'r':
            value = l2p(tau) - (l2p(tp) + l2m(tp)) * std::exp(-alpha * (tp - tau)) +
                    l2m(tau) * std::exp(-two * alpha * (tp - tau));
            return alpha / two / c * E(tau, tm) / E(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_linear(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l1m(tau) - (l1m(tm) + l1p(tm)) * std::exp(-alpha * (tau - tm)) +
                    l1p(tau) * std::exp(-two * alpha * (tau - tm));
            return alpha / two / c * E(tp, tau) / E(tp, tm) * value;
        case 'r':
            value = l1p(tau) - (l1p(tp) + l1m(tp)) * std::exp(-alpha * (tp - tau)) +
                    l1m(tau) * std::exp(-two * alpha * (tp - tau));
            return alpha / two / c * E(tau, tm) / E(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_const(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l0(tau, tm) - l0(tau, tm) * std::exp(-alpha * (tau - tm));
            return alpha / two / c * E(tp, tau) / E(tp, tm) * value;
        case 'r':
            value = l0(tp, tau) - l0(tp, tau) * std::exp(-alpha * (tp - tau));
            return alpha / two / c * E(tau, tm) / E(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_square_t(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l2m(tau) - (l2m(tm) - l2p(tm)) * std::exp(-alpha * (tau - tm)) -
                    l2p(tau) * std::exp(-two * alpha * (tau - tm));
            value *= alpha;
            return alpha / two / c * E(tp, tau) / E(tp, tm) * value;
        case 'r':
            value = l2p(tau) - (l2p(tp) - l2m(tp)) * std::exp(-alpha * (tp - tau)) -
                    l2m(tau) * std::exp(-two * alpha * (tp - tau));
            value *= -alpha;
            return alpha / two / c * E(tau, tm) / E(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_linear_t(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l1m(tau) - (l1m(tm) - l1p(tm)) * std::exp(-alpha * (tau - tm)) -
                    l1p(tau) * std::exp(-two * alpha * (tau - tm));
            value *= alpha;
            return alpha / two / c * E(tp, tau) / E(tp, tm) * value;
        case 'r':
            value = l1p(tau) - (l1p(tp) - l1m(tp)) * std::exp(-alpha * (tp - tau)) -
                    l1m(tau) * std::exp(-two * alpha * (tp - tau));
            value *= -alpha;
            return alpha / two / c * E(tau, tm) / E(tp, tm) * value;
        default:
            return ZEROVALUE;
    }

}

auto AGM::GreenfunctionReactionDiffusion::integrate_const_t(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l0(tau, tm) + l0(tau, tm) * std::exp(-alpha * (tau - tm));
            value *= alpha;
            return alpha / two / c * E(tp, tau) / E(tp, tm) * value;
        case 'r':
            value = l0(tp, tau) + l0(tp, tau) * std::exp(-alpha * (tp - tau));
            value *= -alpha;
            return alpha / two / c * E(tau, tm) / E(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_square_tau(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l2m(tau) - (l2m(tm) + l2p(tm)) * std::exp(-alpha * (tau - tm)) +
                    l2p(tau) * std::exp(-two * alpha * (tau - tm));
            value *= -alpha * F(tp, tau);
            return alpha / two / c / E(tp, tm) * value;
        case 'r':
            value = l2p(tau) - (l2p(tp) + l2m(tp)) * std::exp(-alpha * (tp - tau)) +
                    l2m(tau) * std::exp(-two * alpha * (tp - tau));
            value *= alpha * F(tau, tm);
            return alpha / two / c / E(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_linear_tau(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l1m(tau) - (l1m(tm) + l1p(tm)) * std::exp(-alpha * (tau - tm)) +
                    l1p(tau) * std::exp(-two * alpha * (tau - tm));
            value *= -alpha * F(tp, tau);
            return alpha / two / c / E(tp, tm) * value;
        case 'r':
            value = l1p(tau) - (l1p(tp) + l1m(tp)) * std::exp(-alpha * (tp - tau)) +
                    l1m(tau) * std::exp(-two * alpha * (tp - tau));
            value *= alpha * F(tau, tm);
            return alpha / two / c / E(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_const_tau(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l0(tau, tm) - l0(tau, tm) * std::exp(-alpha * (tau - tm));
            value *= -alpha * F(tp, tau);
            return alpha / two / c / E(tp, tm) * value;
        case 'r':
            value = l0(tp, tau) - l0(tp, tau) * std::exp(-alpha * (tp - tau));
            value *= alpha * F(tau, tm);
            return alpha / two / c / E(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_square_ttau(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l1m(tau) - (l1m(tm) + l1p(tm)) * std::exp(-alpha * (tau - tm)) +
                    l1p(tau) * std::exp(-two * alpha * (tau - tm));
            value = 2 * value - std::pow(tau, 2) * E(tau, tm);
            value *= alpha * F(tp, tau);
            return alpha / two / c / E(tp, tm) * value;
        case 'r':
            value = l1p(tau) - (l1p(tp) + l1m(tp)) * std::exp(-alpha * (tp - tau)) +
                    l1m(tau) * std::exp(-two * alpha * (tp - tau));
            value = 2 * value + std::pow(tau, 2) * E(tp, tau);
            value *= -alpha * F(tau, tm);
            return alpha / two / c / E(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_linear_ttau(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l0(tau, tm) - l0(tau, tm) * std::exp(-alpha * (tau - tm)) - tau * E(tau, tm);
            value *= alpha * F(tp, tau);
            return alpha / two / c / E(tp, tm) * value;
        case 'r':
            value = l0(tp, tau) - l0(tp, tau) * std::exp(-alpha * (tp - tau)) + tau * E(tp, tau);
            value *= -alpha * F(tau, tm);
            return alpha / two / c / E(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_const_ttau(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = -alpha * F(tp, tau) * E(tau, tm);
            return alpha / two / c / E(tp, tm) * value;
        case 'r':
            value = -alpha * F(tau, tm) * E(tp, tau);
            return alpha / two / c / E(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_square_ND(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l2m(tau) - (l2m(tm) - l2p(tm)) * std::exp(-alpha * (tau - tm)) -
                    l2p(tau) * std::exp(-two * alpha * (tau - tm));
            return alpha / two / c * E(tp, tau) / F(tp, tm) * value;
        case 'r':
            value = l2p(tau) - (l2p(tp) + l2m(tp)) * std::exp(-alpha * (tp - tau)) +
                    l2m(tau) * std::exp(-two * alpha * (tp - tau));
            return alpha / two / c * F(tau, tm) / F(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_linear_ND(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l1m(tau) - (l1m(tm) - l1p(tm)) * std::exp(-alpha * (tau - tm)) -
                    l1p(tau) * std::exp(-two * alpha * (tau - tm));
            return alpha / two / c * E(tp, tau) / F(tp, tm) * value;
        case 'r':
            value = l1p(tau) - (l1p(tp) + l1m(tp)) * std::exp(-alpha * (tp - tau)) +
                    l1m(tau) * std::exp(-two * alpha * (tp - tau));
            return alpha / two / c * F(tau, tm) / F(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_const_ND(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l0(tau, tm) + l0(tau, tm) * std::exp(-alpha * (tau - tm));
            return alpha / two / c * E(tp, tau) / F(tp, tm) * value;
        case 'r':
            value = l0(tp, tau) - l0(tp, tau) * std::exp(-alpha * (tp - tau));
            return alpha / two / c * F(tau, tm) / F(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_square_t_ND(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l2m(tau) - (l2m(tm) + l2p(tm)) * std::exp(-alpha * (tau - tm)) +
                    l2p(tau) * std::exp(-two * alpha * (tau - tm));
            value *= alpha;
            return alpha / two / c * E(tp, tau) / F(tp, tm) * value;
        case 'r':
            value = l2p(tau) - (l2p(tp) - l2m(tp)) * std::exp(-alpha * (tp - tau)) -
                    l2m(tau) * std::exp(-two * alpha * (tp - tau));
            value *= -alpha;
            return alpha / two / c * F(tau, tm) / F(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_linear_t_ND(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l1m(tau) - (l1m(tm) + l1p(tm)) * std::exp(-alpha * (tau - tm)) +
                    l1p(tau) * std::exp(-two * alpha * (tau - tm));
            value *= alpha;
            return alpha / two / c * E(tp, tau) / F(tp, tm) * value;
        case 'r':
            value = l1p(tau) - (l1p(tp) - l1m(tp)) * std::exp(-alpha * (tp - tau)) -
                    l1m(tau) * std::exp(-two * alpha * (tp - tau));
            value *= -alpha;
            return alpha / two / c * F(tau, tm) / F(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_const_t_ND(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l0(tau, tm) - l0(tau, tm) * std::exp(-alpha * (tau - tm));
            value *= alpha;
            return alpha / two / c * E(tp, tau) / F(tp, tm) * value;
        case 'r':
            value = l0(tp, tau) + l0(tp, tau) * std::exp(-alpha * (tp - tau));
            value *= -alpha;
            return alpha / two / c * F(tau, tm) / F(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_square_tau_ND(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l2m(tau) - (l2m(tm) - l2p(tm)) * std::exp(-alpha * (tau - tm)) -
                    l2p(tau) * std::exp(-two * alpha * (tau - tm));
            value *= -alpha * F(tp, tau);
            return alpha / two / c / F(tp, tm) * value;
        case 'r':
            value = l2p(tau) - (l2p(tp) + l2m(tp)) * std::exp(-alpha * (tp - tau)) +
                    l2m(tau) * std::exp(-two * alpha * (tp - tau));
            value *= alpha * E(tau, tm);
            return alpha / two / c / F(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_linear_tau_ND(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l1m(tau) - (l1m(tm) - l1p(tm)) * std::exp(-alpha * (tau - tm)) -
                    l1p(tau) * std::exp(-two * alpha * (tau - tm));
            value *= -alpha * F(tp, tau);
            return alpha / two / c / F(tp, tm) * value;
        case 'r':
            value = l1p(tau) - (l1p(tp) + l1m(tp)) * std::exp(-alpha * (tp - tau)) +
                    l1m(tau) * std::exp(-two * alpha * (tp - tau));
            value *= alpha * E(tau, tm);
            return alpha / two / c / F(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_const_tau_ND(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l0(tau, tm) + l0(tau, tm) * std::exp(-alpha * (tau - tm));
            value *= -alpha * F(tp, tau);
            return alpha / two / c / F(tp, tm) * value;
        case 'r':
            value = l0(tp, tau) - l0(tp, tau) * std::exp(-alpha * (tp - tau));
            value *= alpha * E(tau, tm);
            return alpha / two / c / F(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_square_DN(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l2m(tau) - (l2m(tm) + l2p(tm)) * std::exp(-alpha * (tau - tm)) +
                    l2p(tau) * std::exp(-two * alpha * (tau - tm));
            return alpha / two / c * F(tp, tau) / F(tp, tm) * value;
        case 'r':
            value = l2p(tau) - (l2p(tp) - l2m(tp)) * std::exp(-alpha * (tp - tau)) -
                    l2m(tau) * std::exp(-two * alpha * (tp - tau));
            return alpha / two / c * E(tau, tm) / F(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_linear_DN(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l1m(tau) - (l1m(tm) + l1p(tm)) * std::exp(-alpha * (tau - tm)) +
                    l1p(tau) * std::exp(-two * alpha * (tau - tm));
            return alpha / two / c * F(tp, tau) / F(tp, tm) * value;
        case 'r':
            value = l1p(tau) - (l1p(tp) - l1m(tp)) * std::exp(-alpha * (tp - tau)) -
                    l1m(tau) * std::exp(-two * alpha * (tp - tau));
            return alpha / two / c * E(tau, tm) / F(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_const_DN(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l0(tau, tm) - l0(tau, tm) * std::exp(-alpha * (tau - tm));
            return alpha / two / c * F(tp, tau) / F(tp, tm) * value;
        case 'r':
            value = l0(tp, tau) + l0(tp, tau) * std::exp(-alpha * (tp - tau));
            return alpha / two / c * E(tau, tm) / F(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_square_t_DN(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l2m(tau) - (l2m(tm) - l2p(tm)) * std::exp(-alpha * (tau - tm)) -
                    l2p(tau) * std::exp(-two * alpha * (tau - tm));
            value *= alpha;
            return alpha / two / c * F(tp, tau) / F(tp, tm) * value;
        case 'r':
            value = l2p(tau) - (l2p(tp) + l2m(tp)) * std::exp(-alpha * (tp - tau)) +
                    l2m(tau) * std::exp(-two * alpha * (tp - tau));
            value *= -alpha;
            return alpha / two / c * E(tau, tm) / F(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_linear_t_DN(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l1m(tau) - (l1m(tm) - l1p(tm)) * std::exp(-alpha * (tau - tm)) -
                    l1p(tau) * std::exp(-two * alpha * (tau - tm));
            value *= alpha;
            return alpha / two / c * F(tp, tau) / F(tp, tm) * value;
        case 'r':
            value = l1p(tau) - (l1p(tp) + l1m(tp)) * std::exp(-alpha * (tp - tau)) +
                    l1m(tau) * std::exp(-two * alpha * (tp - tau));
            value *= -alpha;
            return alpha / two / c * E(tau, tm) / F(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_const_t_DN(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l0(tau, tm) + l0(tau, tm) * std::exp(-alpha * (tau - tm));
            value *= alpha;
            return alpha / two / c * F(tp, tau) / F(tp, tm) * value;
        case 'r':
            value = l0(tp, tau) - l0(tp, tau) * std::exp(-alpha * (tp - tau));
            value *= -alpha;
            return alpha / two / c * E(tau, tm) / F(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_square_tau_DN(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l2m(tau) - (l2m(tm) + l2p(tm)) * std::exp(-alpha * (tau - tm)) +
                    l2p(tau) * std::exp(-two * alpha * (tau - tm));
            value *= -alpha * E(tp, tau);
            return alpha / two / c / F(tp, tm) * value;
        case 'r':
            value = l2p(tau) - (l2p(tp) - l2m(tp)) * std::exp(-alpha * (tp - tau)) -
                    l2m(tau) * std::exp(-two * alpha * (tp - tau));
            value *= alpha * F(tau, tm);
            return alpha / two / c / F(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_linear_tau_DN(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l1m(tau) - (l1m(tm) + l1p(tm)) * std::exp(-alpha * (tau - tm)) +
                    l1p(tau) * std::exp(-two * alpha * (tau - tm));
            value *= -alpha * E(tp, tau);
            return alpha / two / c / F(tp, tm) * value;
        case 'r':
            value = l1p(tau) - (l1p(tp) - l1m(tp)) * std::exp(-alpha * (tp - tau)) -
                    l1m(tau) * std::exp(-two * alpha * (tp - tau));
            value *= alpha * F(tau, tm);
            return alpha / two / c / F(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

auto AGM::GreenfunctionReactionDiffusion::integrate_const_tau_DN(char i) const -> double {
    double value{};
    switch (i) {
        case 'l':
            value = l0(tau, tm) - l0(tau, tm) * std::exp(-alpha * (tau - tm));
            value *= -alpha * E(tp, tau);
            return alpha / two / c / F(tp, tm) * value;
        case 'r':
            value = l0(tp, tau) + l0(tp, tau) * std::exp(-alpha * (tp - tau));
            value *= alpha * F(tau, tm);
            return alpha / two / c / F(tp, tm) * value;
        default:
            return ZEROVALUE;
    }
}

AGM::GreenfunctionReactionDiffusion::~GreenfunctionReactionDiffusion() = default;


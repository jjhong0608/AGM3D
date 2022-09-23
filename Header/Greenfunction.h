//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#ifndef AGM3D_GREENFUNCTION_H
#define AGM3D_GREENFUNCTION_H

#include "util.h"

static const double twelve = 12.0;
static const double six = 6.0;
static const double three = 3.0;
static const double two = 2.0;
static const double tol_lin_suq = 0.3;

namespace AGM {
    class Greenfunction {
    protected:
        double tm{}, tau{}, tp{}, mpl{}, mpr{};

    public:
        Greenfunction(double tm, double tau, double tp, double mpl, double mpr);

        virtual ~Greenfunction();

        [[nodiscard]] virtual auto integrate_square(char i) const -> double;

        [[nodiscard]] virtual auto integrate_linear(char i) const -> double;

        [[nodiscard]] virtual auto integrate_const(char i) const -> double;

        [[nodiscard]] virtual auto integrate_square_t(char i) const -> double;

        [[nodiscard]] virtual auto integrate_linear_t(char i) const -> double;

        [[nodiscard]] virtual auto integrate_const_t(char i) const -> double;

        [[nodiscard]] virtual auto integrate_square_tau(char i) const -> double;

        [[nodiscard]] virtual auto integrate_linear_tau(char i) const -> double;

        [[nodiscard]] virtual auto integrate_const_tau(char i) const -> double;

        [[nodiscard]] virtual auto integrate_square_ttau(char i) const -> double;

        [[nodiscard]] virtual auto integrate_linear_ttau(char i) const -> double;

        [[nodiscard]] virtual auto integrate_const_ttau(char i) const -> double;

        [[nodiscard]] virtual auto green_function(double t) const -> double;

        [[nodiscard]] virtual auto green_function_t(double t) const -> double;

        [[nodiscard]] virtual auto green_function_tau(double t) const -> double;

        [[nodiscard]] virtual auto green_function_ttau(double t) const -> double;

        [[nodiscard]] virtual auto green_integral(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_t(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_tau(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_ttau(char pos) const -> double;

        [[nodiscard]] auto green_integral_square(char pos) const -> double;

        [[nodiscard]] auto green_integral_t_square(char pos) const -> double;

        [[nodiscard]] auto green_integral_tau_square(char pos) const -> double;

        [[nodiscard]] auto green_integral_ttau_square(char pos) const -> double;

        [[nodiscard]] auto green_integral_linear(char pos) const -> double;

        [[nodiscard]] auto green_integral_t_linear(char pos) const -> double;

        [[nodiscard]] auto green_integral_tau_linear(char pos) const -> double;

        [[nodiscard]] auto green_integral_ttau_linear(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_square_ND(char i) const -> double;

        [[nodiscard]] virtual auto integrate_linear_ND(char i) const -> double;

        [[nodiscard]] virtual auto integrate_const_ND(char i) const -> double;

        [[nodiscard]] virtual auto integrate_square_t_ND(char i) const -> double;

        [[nodiscard]] virtual auto integrate_linear_t_ND(char i) const -> double;

        [[nodiscard]] virtual auto integrate_const_t_ND(char i) const -> double;

        [[nodiscard]] virtual auto integrate_square_tau_ND(char i) const -> double;

        [[nodiscard]] virtual auto integrate_linear_tau_ND(char i) const -> double;

        [[nodiscard]] virtual auto integrate_const_tau_ND(char i) const -> double;

        [[nodiscard]] static auto integrate_square_ttau_ND(char i) -> double;

        [[nodiscard]] static auto integrate_linear_ttau_ND(char i) -> double;

        [[nodiscard]] static auto integrate_const_ttau_ND(char i) -> double;

        [[nodiscard]] virtual auto green_function_ND(double t) const -> double;

        [[nodiscard]] virtual auto green_function_t_ND(double t) const -> double;

        [[nodiscard]] virtual auto green_function_tau_ND(double t) const -> double;

        [[nodiscard]] virtual auto green_function_ttau_ND(double t) const -> double;

        [[nodiscard]] virtual auto green_integral_ND(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_t_ND(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_tau_ND(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_ttau_ND(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_square_ND(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_t_square_ND(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_tau_square_ND(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_ttau_square_ND(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_linear_ND(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_t_linear_ND(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_tau_linear_ND(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_ttau_linear_ND(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_square_DN(char i) const -> double;

        [[nodiscard]] virtual auto integrate_linear_DN(char i) const -> double;

        [[nodiscard]] virtual auto integrate_const_DN(char i) const -> double;

        [[nodiscard]] virtual auto integrate_square_t_DN(char i) const -> double;

        [[nodiscard]] virtual auto integrate_linear_t_DN(char i) const -> double;

        [[nodiscard]] virtual auto integrate_const_t_DN(char i) const -> double;

        [[nodiscard]] virtual auto integrate_square_tau_DN(char i) const -> double;

        [[nodiscard]] virtual auto integrate_linear_tau_DN(char i) const -> double;

        [[nodiscard]] virtual auto integrate_const_tau_DN(char i) const -> double;

        [[nodiscard]] static auto integrate_square_ttau_DN(char i) -> double;

        [[nodiscard]] static auto integrate_linear_ttau_DN(char i) -> double;

        [[nodiscard]] static auto integrate_const_ttau_DN(char i) -> double;

        [[nodiscard]] virtual auto green_function_DN(double t) const -> double;

        [[nodiscard]] virtual auto green_function_t_DN(double t) const -> double;

        [[nodiscard]] virtual auto green_function_tau_DN(double t) const -> double;

        [[nodiscard]] virtual auto green_function_ttau_DN(double t) const -> double;

        [[nodiscard]] virtual auto green_integral_DN(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_t_DN(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_tau_DN(char pos) const -> double;

        [[nodiscard]] auto green_integral_ttau_DN(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_square_DN(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_t_square_DN(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_tau_square_DN(char pos) const -> double;

        [[nodiscard]] auto green_integral_ttau_square_DN(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_linear_DN(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_t_linear_DN(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_tau_linear_DN(char pos) const -> double;

        [[nodiscard]] auto green_integral_ttau_linear_DN(char pos) const -> double;
    };
}

#endif //AGM3D_GREENFUNCTION_H

//
// Created by NIMS-JUNHONG on 2022/01/03.
//

#ifndef AGM3D_GREENFUNCTIONREACTIONDIFFUSION_H
#define AGM3D_GREENFUNCTIONREACTIONDIFFUSION_H

#include "Greenfunction.h"

namespace AGM {
    class GreenfunctionReactionDiffusion : public Greenfunction {
    protected:
        double c{};
        double alpha{};
    public:
        GreenfunctionReactionDiffusion(double tm, double tau, double tp, double mpl, double mpr, double c);

        ~GreenfunctionReactionDiffusion() override;

        [[nodiscard]] double E(double a, double b) const;

        [[nodiscard]] double F(double a, double b) const;

        [[nodiscard]] double l2p(double s) const;

        [[nodiscard]] double l2m(double s) const;

        [[nodiscard]] double l1p(double s) const;

        [[nodiscard]] double l1m(double s) const;

        [[nodiscard]] double l0(double a, double b) const;

        [[nodiscard]] double green_function(double t) const override;

        [[nodiscard]] double green_function_t(double t) const override;

        [[nodiscard]] double green_function_tau(double t) const override;

        [[nodiscard]] double green_function_ttau(double t) const override;

        [[nodiscard]] double green_function_ND(double t) const override;

        [[nodiscard]] double green_function_t_ND(double t) const override;

        [[nodiscard]] double green_function_tau_ND(double t) const override;

        [[nodiscard]] double green_function_ttau_ND(double t) const override;

        [[nodiscard]] double green_function_DN(double t) const override;

        [[nodiscard]] double green_function_t_DN(double t) const override;

        [[nodiscard]] double green_function_tau_DN(double t) const override;

        [[nodiscard]] double green_function_ttau_DN(double t) const override;

        [[nodiscard]] double integrate_square(char i) const override;

        [[nodiscard]] double integrate_linear(char i) const override;

        [[nodiscard]] double integrate_const(char i) const override;

        [[nodiscard]] double integrate_square_t(char i) const override;

        [[nodiscard]] double integrate_linear_t(char i) const override;

        [[nodiscard]] double integrate_const_t(char i) const override;

        [[nodiscard]] double integrate_square_tau(char i) const override;

        [[nodiscard]] double integrate_linear_tau(char i) const override;

        [[nodiscard]] double integrate_const_tau(char i) const override;

        [[nodiscard]] double integrate_square_ttau(char i) const override;

        [[nodiscard]] double integrate_linear_ttau(char i) const override;

        [[nodiscard]] double integrate_const_ttau(char i) const override;

        [[nodiscard]] double integrate_square_ND(char i) const override;

        [[nodiscard]] double integrate_linear_ND(char i) const override;

        [[nodiscard]] double integrate_const_ND(char i) const override;

        [[nodiscard]] double integrate_square_t_ND(char i) const override;

        [[nodiscard]] double integrate_linear_t_ND(char i) const override;

        [[nodiscard]] double integrate_const_t_ND(char i) const override;

        [[nodiscard]] double integrate_square_tau_ND(char i) const override;

        [[nodiscard]] double integrate_linear_tau_ND(char i) const override;

        [[nodiscard]] double integrate_const_tau_ND(char i) const override;

        [[nodiscard]] double integrate_square_DN(char i) const override;

        [[nodiscard]] double integrate_linear_DN(char i) const override;

        [[nodiscard]] double integrate_const_DN(char i) const override;

        [[nodiscard]] double integrate_square_t_DN(char i) const override;

        [[nodiscard]] double integrate_linear_t_DN(char i) const override;

        [[nodiscard]] double integrate_const_t_DN(char i) const override;

        [[nodiscard]] double integrate_square_tau_DN(char i) const override;

        [[nodiscard]] double integrate_linear_tau_DN(char i) const override;

        [[nodiscard]] double integrate_const_tau_DN(char i) const override;
    };
}


#endif //AGM3D_GREENFUNCTIONREACTIONDIFFUSION_H

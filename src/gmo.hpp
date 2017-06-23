#ifndef SRC_GMO_HPP
#define SRC_GMO_HPP

#include <stan/math/mix/mat.hpp>
#include <util.hpp>
#include <iostream>

namespace mygmo {

  template <typename F>
  void GMO(const F& log_prob, int num_locals, int num_globals) {
    using namespace stan::math;
    using Eigen::Matrix;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    // construct a function to hand to the optimizer that:
    // 1. Given data (y, sigma) and a guess of eta:

    double fx; // result not used
    MatrixXd H; // Also set by outer_opt
    VectorXd locals = VectorXd::Zero(num_locals); // set by outer_opt

    auto outer_opt = [&](VectorXd& globals) -> VectorXd {
      VectorXd grad;
      VectorXd grad_tr_MH;

      //      1. optimizes and finds the argmax of eta
      auto local_prob = [&](auto& locals) {
        return log_prob(locals, globals);
      };
      gradient_descent(local_prob, locals);

      std::cout << "mu: " << globals[0] << "\ttau: " << globals[1] << std::endl;

      //      2. returns the approximated gradient (the g1-g5 stuff)
      auto mixed_prob = [&](auto& locals_and_globals) {
        return log_prob(locals_and_globals.head(locals.size()).eval(),
                        locals_and_globals.tail(globals.size()).eval());
      };
      VectorXd locals_and_globals(locals.size() + globals.size());
      locals_and_globals << locals, globals;

      hessian(mixed_prob, locals_and_globals, fx, grad, H);
      VectorXd g1 = grad.head(num_locals);
      VectorXd g3 = grad.tail(num_globals);

      H = -H;
      MatrixXd locals_hessian_inv = H.topLeftCorner(num_locals, num_locals).inverse();
      grad_tr_mat_times_hessian(mixed_prob, locals_and_globals, locals_hessian_inv, grad_tr_MH);
      grad_tr_MH = -grad_tr_MH;
      VectorXd g2 = -0.5 * grad_tr_MH.head(num_locals);
      VectorXd g4 = -0.5 * grad_tr_MH.tail(num_globals);
      MatrixXd g5 = -(locals_hessian_inv * H.topRightCorner(locals.size(), globals.size()));
      return ((g5.transpose() * (g1 + g2)) + g3 + g4);
    };


    // 2. Use gradient descent to iterate the above function and find the
    //    max for the global vars
    VectorXd globals = VectorXd::Zero(num_globals);
    gradient_descent_manual_gradient(outer_opt, globals);

    // 3. Compute the curvature 2nd derivative at the max for globals
    std::cout << "Curvature matrix: " << std::endl;
    std::cout << H.bottomRightCorner(globals.size(), globals.size()) << std::endl;
  }
}
#endif

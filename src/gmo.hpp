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

      // std::cout << "g1: "; print_mat(g1);
      // std::cout << "g2: "; print_mat(g2);
      // std::cout << "g3: " << g3[0] << ", " << g3[1] << std::endl;
      // std::cout << "g4: " << g4[0] << ", " << g4[1] << std::endl;
      // std::cout << "g5: " << std::endl;
      // std::cout << g5.transpose() << std::endl;

      VectorXd grad_approx =  ((g5.transpose() * (g1 + g2)) + g3 + g4);


      // ====================================== TEST G5 ----
      // auto g5_fd = [&](const auto& globals) {
      //   VectorXd locals_ = locals;
      //   gradient_descent([&](const auto& locals) {
      //       return log_prob(locals, globals);
      //     }, locals_);
      //   return locals_;
      // };
      // std::cout << "finite diff of g5" << std::endl;
      // std::cout << finite_diff(g5_fd, globals) << std::endl;


      // ============================================================

      // // g1 and g3
      // std::cout << "gradient of log prob" << std::endl;
      // std::cout << grad << std::endl;
      // std::cout << "finite diff of log prob" << std::endl;
      // std::cout << finite_diff(mixed_prob, locals_and_globals) << std::endl;

      // total function
      // std::cout << "grad approx dmu: " << grad_approx[0] << "\tdtau: " << grad_approx[1] << std::endl;
      // auto reopt_for_globals_prob = [&](const auto& globals) {
      //   VectorXd locals_ = locals;
      //   gradient_descent(local_prob, locals_);

      //   MatrixXd H;
      //   VectorXd locals_and_globals(locals.size() + globals.size());
      //   locals_and_globals << locals_, globals;
      //   hessian(mixed_prob, locals_and_globals, fx, grad, H);
      //   H = -H;
      //   MatrixXd locals_hessian = H.topLeftCorner(num_locals, num_locals);

      //   return log_prob(locals_, globals) - 0.5*log_det(locals_hessian);
      // };
      // auto fdlp = finite_diff(reopt_for_globals_prob, globals);
      // std::cout << "finite diff dmu: " << fdlp[0] << "\tdtau: " << fdlp[1] << std::endl;

      // =========================================================TEST
      // auto log_det_test_hessian = [&](const auto& globals) {
      //   MatrixXd H;
      //   VectorXd locals_and_globals(locals.size() + globals.size());
      //   locals_and_globals << locals, globals;
      //   hessian(mixed_prob, locals_and_globals, fx, grad, H);
      //   H = -H;
      //   MatrixXd locals_hessian = H.topLeftCorner(num_locals, num_locals);
      //   return log_det(locals_hessian);
      // };
      // auto fdlogdet = finite_diff(log_det_test_hessian, globals) ;
      // std::cout << "finite diff of log det: " << fdlogdet[0] << ", " << fdlogdet[1] << std::endl;
      // std::cout << "    grad tr mh globals: " << grad_tr_MH[8] << ", " << grad_tr_MH[9]<< std::endl;


      // auto locals_log_det_test_hessian = [&](const auto& locals) {
      //   MatrixXd H;
      //   VectorXd locals_and_globals(locals.size() + globals.size());
      //   locals_and_globals << locals, globals;
      //   hessian(mixed_prob, locals_and_globals, fx, grad, H);
      //   H = -H;
      //   MatrixXd locals_hessian = H.topLeftCorner(num_locals, num_locals);
      //   return log_det(locals_hessian);
      // };
      // std::cout << "LBLOCK finite diff of log det" << std::endl;
      // print_mat(finite_diff(locals_log_det_test_hessian, locals));
      // std::cout << "LBLOCK grad tr mh" << std::endl;
      // print_mat(grad_tr_MH);
      // std::cout << "LBLOCK done" << std::endl;
      // //=======================================================END TEST



      return grad_approx;
    };


    // 2. Use gradient descent to iterate the above function and find the
    //    max for mu and tau
    VectorXd globals = VectorXd::Zero(num_globals);
    gradient_descent_manual_gradient(outer_opt, globals);

    // 3. Compute the curvature 2nd derivative at the max for globals mu, tau
    std::cout << "Curvature matrix: " << std::endl;
    std::cout << H.bottomRightCorner(globals.size(), globals.size()) << std::endl;
  }
}
#endif

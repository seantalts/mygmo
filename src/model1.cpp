#include <stan/math/mix/mat.hpp>
#include <iostream>

namespace {
  using namespace stan::math;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  const double alpha = 0.1;

  template <typename F, typename T>
  Matrix<T, -1, 1> gradient_descent(const F& f, Matrix<T, -1, 1> theta) {
    T fx;
    Matrix<T, -1, 1> grad;
    for (int i = 0; i < 1000; ++i) {
      gradient(f, theta, fx, grad);
      theta += -alpha * grad;
    }
    return theta;
  }

  template <typename F, typename T>
  Matrix<T, -1, 1> gradient_descent_manual_gradient(const F& f, Matrix<T, -1, 1> theta) {
    for (int i = 0; i < 1000; ++i) {
      Matrix<T, -1, 1> grad = f(theta);
      theta += -alpha * grad;
      if (0 == (i % 100)) {
        std::cout << i << std::endl;
      }
    }
    std::cout << "done with gdmg" << std::endl;
    return theta;
  }

  template <typename T_local, typename T_global,
            typename T_return = typename stan::return_type<T_local, T_global>::type>
  T_return log_prob(const Matrix<T_local, -1, 1>& eta,
                    const T_global& mu,
                    const T_global& tau,
                    const VectorXd& y,
                    const VectorXd& sigma) {
    T_return lp(0.0);

    Matrix<T_return, -1, 1> theta = tau * eta;
    for (int i = 0; i < theta.size(); ++i)
      theta(i) += mu;

    lp += normal_lpdf<true>(eta, 0, 1);
    lp += normal_lpdf<true>(y, theta, sigma);

    return lp;
  }

  VectorXd inner_optimize(double mu, double tau,
                          const VectorXd& eta,
                          const VectorXd& y, const VectorXd& sigma) {
    auto f = [&](const auto& eta) {
      return log_prob(eta, mu, tau, y, sigma);
    };

    return gradient_descent(f, eta);
  }

  void GMO(const VectorXd& y, const VectorXd& sigma) {
    // construct a function to hand to the optimizer that:
    // 1. Given data (y, sigma) and a guess of eta:
    //      1. optimizes and finds the argmax of eta
    //      2. returns the approximated gradient (the g1-g5 stuff)
    auto inner_opt = [&](VectorXd globals) -> VectorXd {
      VectorXd eta = VectorXd::Zero(y.size());
      eta = inner_optimize(globals[0], globals[1], eta, y, sigma);

      // g1 & g3
      double fx; // not used
      VectorXd g1;
      Matrix<double, -1, -1> H;
      auto local_prob = [&](const auto& eta) {
        return log_prob(eta, globals[0], globals[1], y, sigma);
      };
      hessian(local_prob, eta, fx, g1, H);
      g1 = -g1;
      VectorXd grad_tr_MH;
      grad_tr_mat_times_hessian(local_prob, eta, H, grad_tr_MH);
      auto g3 = -0.5 * grad_tr_MH;

      // g2
      VectorXd g2;
      auto global_prob = [&](const auto& globals) {
        return log_prob(eta, globals[0], globals[1], y, sigma);
      };
      gradient(global_prob, globals, fx, g2);
      g2 = -g2;

      MatrixXd g5 = MatrixXd::Ones(globals.size(), eta.size());

      //std::cout << "g1: " << g1.rows() << "x" << g1.cols() << std::endl;
      //std::cout << "g2: " << g2.rows() << "x" << g2.cols() << std::endl;
      //std::cout << "g3: " << g3.rows() << "x" << g3.cols() << std::endl;
      //std::cout << "g5: " << g5.rows() << "x" << g5.cols() << std::endl;

      auto term1 = g1 + g3;
      auto term2 = g5 * term1;
      auto result = term2 + g2;
      return result;
    };


    // 2. Use gradient descent to iterate the above function and find the
    //    max for mu and tau
    VectorXd globals = VectorXd::Zero(2);
    globals = gradient_descent_manual_gradient(inner_opt, globals);
    std::cout << "mu: " << globals[0] << " tau: " << globals[1] << std::endl;
    // 3. Compute the curvature 2nd derivative at the max for globals mu, tau
    // 4. print out mu, tau, eta, curvature
  }
}

int main(int argc, const char* argv[]) {
  // Set up data
  int J = 8;
  VectorXd y(J);
  y << 28, 8, -3, 7, -1, 1, 18, 12;
  VectorXd sigma(J);
  sigma << 15, 10, 16, 11,  9, 11, 10, 18;

  std::cout << "Calling GMO " << std::endl;
  GMO(y, sigma);

  return 0;
}

/*
  Questions:
  1. Why does reader.hpp have scalar_constrain(lp) that always does nothing? when does it ever do something? why do we often pass it something only if jacobian is true?
  2. Why does the model initialize and then fill a variable?
 */

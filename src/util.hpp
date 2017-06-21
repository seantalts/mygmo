#ifndef SRC_UTIL_HPP
#define SRC_UTIL_HPP

#include <stan/math/mix/mat.hpp>

namespace mygmo {
  using namespace stan::math;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  const double alpha = 1.0;
  const double epsilon = 1e-6;
  const int iterations = 20000;

  template <typename F>
  void gradient_descent_manual_gradient(const F& f, VectorXd& theta) {
    VectorXd grad = VectorXd::Ones(theta.size());
    int i = 0;
    for (; i < iterations && grad.norm() > epsilon; ++i) {
      grad = f(theta);
      theta += alpha * grad; // ASCENT
    }
    if (i == iterations)
      std::cout << "Did not converge." << std::endl;
  }
  template <typename F>
  void gradient_descent(const F& f, VectorXd& theta) {
    gradient_descent_manual_gradient([&f](const auto& theta){
        double fx;
        VectorXd grad;
        gradient(f, theta, fx, grad);
        return grad;
      }, theta);
  }

  MatrixXd finite_diff(const std::function<VectorXd(VectorXd&)>& f, const VectorXd& theta) {
    VectorXd perturbed = theta;
    MatrixXd grad;
    for (int k = 0; k < theta.size(); ++k) {
      perturbed(k) = theta(k) + epsilon;
      VectorXd fx_plus = f(perturbed);
      if (k == 0) {
        grad.resize(theta.size(), fx_plus.size());
      }

      perturbed(k) = theta(k) - epsilon;
      VectorXd fx_minus = f(perturbed);

      grad.row(k) = (fx_plus - fx_minus) / (2 * epsilon);
      perturbed(k) = theta(k);
    }
    return grad;
  }

  VectorXd finite_diff(const std::function<double(VectorXd&)> & f, const VectorXd& theta) {
    VectorXd perturbed = theta;
    VectorXd grad(theta.size());
    for (int k = 0; k < theta.size(); ++k) {
      perturbed(k) = theta(k) + epsilon;
      double fx_plus = f(perturbed);

      perturbed(k) = theta(k) - epsilon;
      double fx_minus = f(perturbed);

      grad(k) = (fx_plus - fx_minus) / (2 * epsilon);
      perturbed(k) = theta(k);
    }
    return grad;
  }

  double log_det(MatrixXd A) {
    Eigen::LLT<MatrixXd> LLT_of_A(A);
    MatrixXd L = LLT_of_A.matrixL();
    double sum = 0;
    for (int i = 0; i < L.rows(); ++i) {
      sum += log(L(i, i));
    }
    return 2 * sum;
  }

  template <typename T, int R, int C>
  void print_mat(Matrix<T, R, C>& mat) {
    for (int i = 0; i < mat.size(); ++i)
      std::cout << mat(i) << ", ";
    std::cout << std::endl;
  }
}
#endif

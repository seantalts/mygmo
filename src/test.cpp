#include <util.hpp>
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;

namespace mygmo {
  template <typename T_local, typename T_global>
  T_local test_f_1(const Matrix<T_local, -1, 1>& locals,
                   const Matrix<T_global, -1, 1>& globals) {
    T_local sum(0.0);
    for (int i = 0; i < locals.size(); ++i)
      sum += locals(i) * locals(i) * locals(i);
    for (int i = 0; i < globals.size(); ++i)
      sum += globals(i) * globals(i) * globals(i);
    return sum;
  }
}

int main() {
  using namespace mygmo;
  using namespace stan::math;

  VectorXd globals = VectorXd::Zero(2);
  auto f = [&](const auto& locals) {
    return test_f_1(locals, globals);
  };
  VectorXd locals(3);
  locals << -0.5, 0.2, 0.5;
  auto grad = finite_diff(f, locals);
  std::cout << "grad 1 "<< std::endl;
  std::cout << grad << std::endl;

  auto ldh = [&](const auto& locals) {
    double fx;
    VectorXd grad;
    MatrixXd H;
    hessian(f, locals, fx, grad, H);
    std::cout << H << std::endl;
    return log_det(H);
  };
  grad = finite_diff(ldh, locals);
  std::cout << "grad 2 "<< std::endl;
  std::cout << grad << std::endl;
  return 0;
}

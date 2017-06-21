#include <stan/math/mix/mat.hpp>
#include <util.hpp>
#include <gmo.hpp>
#include <iostream>

namespace {
  using namespace stan::math;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  template <typename T_local, typename T_global,
            typename T_return = typename stan::return_type<T_local, T_global>::type>
  T_return log_prob(const Matrix<T_local, -1, 1>& eta,
                    const T_global& mu,
                    const T_global& log_tau,
                    const VectorXd& y,
                    const VectorXd& sigma) {
    T_return lp(0.0);
    Matrix<T_return, -1, 1> theta(eta.size());

    for (int i = 0; i < theta.size(); ++i)
      theta(i) = eta(i) * exp(log_tau) + mu;

    lp += normal_lpdf(eta, 0, 10);
    lp += normal_lpdf(mu, 0, 10);
    lp += normal_lpdf(exp(log_tau), 0, 10);
    lp += normal_lpdf(y, theta, sigma);

    // adjust for change of vars:
    // log | d/d(log_tau) exp(log_tau) | = log_tau
    lp += log_tau;

    return lp;
  }
}

int main(int argc, const char* argv[]) {
  using Eigen::VectorXd;
  // Set up data
  int J = 8;
  VectorXd y(J); y << 28, 8, -3, 7, -1, 1, 18, 12;
  VectorXd sigma(J); sigma << 15, 10, 16, 11,  9, 11, 10, 18;

  std::cout << "Calling GMO " << std::endl;
  auto f = [&](const auto& locals, const auto& globals) {
    return log_prob(locals, globals[0], globals[1], y, sigma);
  };
  mygmo::GMO(f, J, 2);

  return 0;
}

/*
  Questions:
  1. Why does reader.hpp have scalar_constrain(lp) that always does nothing? when does it ever do something? why do we often pass it something only if jacobian is true?
  2. Why does the model initialize and then fill a variable?
 */

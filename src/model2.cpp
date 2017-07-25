#include <stan/math/mix/mat.hpp>
#include <util.hpp>
#include <gmo.hpp>
#include <iostream>

namespace {
  using namespace stan::math;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using std::vector;

  struct hls {
    vector<int> group_for_indiv;
    MatrixXd indiv_predictors;
    vector<VectorXd> group_predictors;
    VectorXd outcomes;
    int num_groups;
    size_t num_indivs;

    hls(vector<int> group_for_indiv,
        MatrixXd indiv_predictors,
        vector<VectorXd> group_predictors,
        VectorXd outcomes):
      group_for_indiv(group_for_indiv),
      indiv_predictors(indiv_predictors),
      group_predictors(group_predictors),
      outcomes(outcomes),
      num_groups(group_predictors.size()),
      num_indivs(group_for_indiv.size()) {}
    template <typename T_local, typename T_global,
              typename T_return = typename stan::return_type<T_local, T_global>::type>
    T_return operator()(const Matrix<T_local, -1, -1>& group_coeffs,
                        const Matrix<T_local, -1, -1>& indiv_coeffs_by_group,
                        const Matrix<T_global, -1, -1>& prior_correlation,
                        const Matrix<T_global, -1, 1>& prior_scale,
                        double sigma) {
      T_return lp(0.0);
      lp += cauchy_lpdf(prior_scale, 0, 2.5);
      lp += lkj_corr_lpdf(prior_correlation, 2);
      lp += normal_lpdf(to_vector(group_coeffs, 0, 5));
      {
        VectorXd u_gamma(num_groups);
        for (int j = 0; j < num_groups; ++j)
          u_gamma[j] = group_predictors[j] * group_coeffs;
        lp += multi_normal_lpdf(indiv_coeffs_by_group,
                                u_gamma,
                                quad_form_diag(prior_correlation, prior_scale));
      }
      for (size_t n = 0; n < num_indivs; ++n)
        lp += normal_lpdf(indiv_predictors[n]
                          * indiv_coeffs_by_group[group_for_indiv[n]],
                          sigma);
    }
  };
}

int main(int argc, const char* argv[]) {
  using Eigen::VectorXd;
  // Set up data
  int N = 8;
  int K = 2;
  int J = 2;
  int L = 3;
  MatrixXd x(N, K);
  x << -6.3228490, 13.5336893, -1.5675845, 0.3730190, -4.4344975, 5.0613612,
    -2.4631817, -5.0399788, -8.3380768, 0.1271512, 19.9426983, 6.1103013, -3.0379787,
    -11.0154858, -25.3439396, -3.0595938;
  VectorXd u1(L);
  VectorXd u2(L);
  u1 << -3.8335517, -11.9532094, 4.4341937;
  u2 << -2.2207795, 10.4194008, -5.0440609;
  vector<VectorXd> u {u1, u2};
  VectorXd y(N);
  y << 1527.6682088, -1267.6828239, 1314.5835858, 488.7249873, 758.1054782,
    777.1681065, -2091.1838994, -704.2232294;
  hls model({1, 1, 2, 2, 1, 1, 2, 2},
            x, u, y);

  std::cout << "Calling GMO for hierarchical logistic regr" << std::endl;
  auto f = [&](const auto& locals, const auto& globals) {
    // Let's consider the coeffs to be locals we want to marginalize over.
    // locals
    MatrixXd group_coeffs(L, K);
    MatrixXd indiv_coeffs_by_group(L, K);
    // globals
    VectorXd prior_corr(K);
    VectorXd prior_scale(K);
  };
  //mygmo::GMO(model, );
  // mygmo::GMO(f, 50, 2);

  return 0;
}

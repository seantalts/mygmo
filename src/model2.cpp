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

  void gen_data() {
    using std::cout;
    using std::endl;

    int num_indivs = 100;
    int num_groups = 10;
    int num_group_features = 5;
    int num_indiv_features = 5;
    vector<int> group_for_indiv(num_indivs);
    for (int i = 0; i < num_groups; ++i)
      for (int j = 0; j < num_indivs / num_groups; ++j)
        group_for_indiv[i*10+j] = i;
    cout << group_for_indiv[23] << endl;
  }
}

int main(int argc, const char* argv[]) {
  using Eigen::VectorXd;
  // Set up data

  gen_data();

  // hls model(vector<int>());
  // std::cout << "Calling GMO for hierarchical logistic regr" << std::endl;
  // auto f = [&](const auto& locals, const auto& globals) {
  //   return model(locals, globals[0], globals[1], y, sigma);
  // };
  // mygmo::GMO(f, 50, 2);

  return 0;
}

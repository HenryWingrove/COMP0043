#include<iostream>
#include<vector>
#include<cmath>

// https://www.quantstart.com/articles/Asian-option-pricing-with-C-via-Monte-Carlo-Methods/
double gaussian_box_muller() {
  double x = 0.0;
  double y = 0.0;
  double euclid_sq = 0.0;

  do {
    x = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
    y = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
    euclid_sq = x*x + y*y;
  } while (euclid_sq >= 1.0);

  return x*sqrt(-2*log(euclid_sq)/euclid_sq);
}


class GeometricBrownianMotion {
public:
    GeometricBrownianMotion(double spot, double vol, double r, double q, double dt)
                                                        : spot_(spot), vol_(vol), r_(r), q_(q), dt_(dt)

    std::vector<double> get_sde(int);

private:
    double spot_;
    double vol_;
    double r_;
    double q_;
    double dt_;
};

GeometricBrownianMotion::get_sde(int nsteps) {
    std::vector<double> dX(nsteps, 1.0);
    dX[0] = spot_;
    for (size_t i = 1; i < dX.size(); i++){
        double W = gaussian_box_muller();
        dX[i] = dX[i-1] * exp(dt_*(r_ - 0.5*pow(v_, 2)) * sqrt(pow(v_, 2)*W));
    }
}
#include "openmc/wmp.h"

#include "openmc/constants.h"
#include "openmc/hdf5_interface.h"
#include "openmc/math_functions.h"

#include <cmath>

namespace openmc {

WindowedMultipole::WindowedMultipole(hid_t group)
{
  // Get name of nuclide from group, removing leading '/'
  name_ = object_name(group).substr(1);

  // Read scalar values.
  read_dataset(group, "spacing", spacing_);
  read_dataset(group, "sqrtAWR", sqrt_awr_);
  read_dataset(group, "E_min", E_min_);
  read_dataset(group, "E_max", E_max_);

  // Read the "data" array.  Use its shape to figure out the number of poles
  // and residue types in this data.
  read_dataset(group, "data", data_);
  int n_residues = data_.shape()[1] - 1;

  // Check to see if this data includes fission residues.
  fissionable_ = (n_residues == 3);

  // Read the "windows" array and use its shape to figure out the number of
  // windows.
  read_dataset(group, "windows", windows_);
  int n_windows = windows_.shape()[0];

  // Read the "broaden_poly" arrays.
  read_dataset(group, "broaden_poly", broaden_poly_);
  if (n_windows != broaden_poly_.shape()[0]) {
    fatal_error("broaden_poly array shape is not consistent with the windows "
      "array shape in WMP library for " + name_ + ".");
  }

  // Read the "curvefit" array.
  read_dataset(group, "curvefit", curvefit_);
  if (n_windows != broaden_poly_.shape()[0]) {
    fatal_error("curvefit array shape is not consistent with the windows "
      "array shape in WMP library for " + name_ + ".");
  }
  fit_order_ = curvefit_.shape()[1] - 1;
}

std::tuple<double, double, double>
WindowedMultipole::evaluate(double E, double sqrtkT)
{
  using namespace std::complex_literals;

  // ==========================================================================
  // Bookkeeping

  // Define some frequently used variables.
  double sqrtE = std::sqrt(E);
  double invE = 1.0 / E;

  // Locate window containing energy
  int i_window = (sqrtE - std::sqrt(E_min_)) / spacing_;
  int startw = windows_(i_window, 0) - 1;
  int endw = windows_(i_window, 1) - 1;

  // Initialize the ouptut cross sections
  double sig_s = 0.0;
  double sig_a = 0.0;
  double sig_f = 0.0;

  // ==========================================================================
  // Add the contribution from the curvefit polynomial.

  if (sqrtkT > 0.0 && broaden_poly_(i_window)) {
    // Broaden the curvefit.
    double dopp = sqrt_awr_ / sqrtkT;
    std::vector<double> broadened_polynomials(fit_order_ + 1);
    broaden_wmp_polynomials(E, dopp, fit_order_ + 1, broadened_polynomials.data());
    for (int i_poly = 0; i_poly < fit_order_ + 1; ++i_poly) {
      sig_s += curvefit_(i_window, i_poly, FIT_S) * broadened_polynomials[i_poly];
      sig_a += curvefit_(i_window, i_poly, FIT_A) * broadened_polynomials[i_poly];
      if (fissionable_) {
        sig_f += curvefit_(i_window, i_poly, FIT_F) * broadened_polynomials[i_poly];
      }
    }
  } else {
    // Evaluate as if it were a polynomial
    double temp = invE;
    for (int i_poly = 0; i_poly < fit_order_ + 1; ++i_poly) {
      sig_s += curvefit_(i_window, i_poly, FIT_S) * temp;
      sig_a += curvefit_(i_window, i_poly, FIT_A) * temp;
      if (fissionable_) {
        sig_f += curvefit_(i_window, i_poly, FIT_F) * temp;
      }
      temp *= sqrtE;
    }
  }

  // ==========================================================================
  // Add the contribution from the poles in this window.

  if (sqrtkT == 0.0) {
    // If at 0K, use asymptotic form.
    for (int i_pole = startw; i_pole <= endw; ++i_pole) {
      std::complex<double> psi_chi = -1.0i / (data_(i_pole, MP_EA) - sqrtE);
      std::complex<double> c_temp = psi_chi / E;
      sig_s += (data_(i_pole, MP_RS) * c_temp).real();
      sig_a += (data_(i_pole, MP_RA) * c_temp).real();
      if (fissionable_) {
        sig_f += (data_(i_pole, MP_RF) * c_temp).real();
      }
    }
  } else {
    // At temperature, use Faddeeva function-based form.
    double dopp = sqrt_awr_ / sqrtkT;
    if (endw >= startw) {
      for (int i_pole = startw; i_pole <= endw; ++i_pole) {
        std::complex<double> z = (sqrtE - data_(i_pole, MP_EA)) * dopp;
        std::complex<double> w_val = faddeeva(z) * dopp * invE * SQRT_PI;
        sig_s += (data_(i_pole, MP_RS) * w_val).real();
        sig_a += (data_(i_pole, MP_RA) * w_val).real();
        if (fissionable_) {
          sig_f += (data_(i_pole, MP_RF) * w_val).real();
        }
      }
    }
  }

  return std::make_tuple(sig_s, sig_a, sig_f);
}

std::tuple<double, double, double>
WindowedMultipole::evaluate_deriv(double E, double sqrtkT)
{
  // ==========================================================================
  // Bookkeeping

  // Define some frequently used variables.
  double sqrtE = std::sqrt(E);
  double invE = 1.0 / E;
  double T = sqrtkT*sqrtkT / K_BOLTZMANN;

  if (sqrtkT == 0.0) {
    fatal_error("Windowed multipole temperature derivatives are not implemented"
      " for 0 Kelvin cross sections.");
  }

  // Locate window containing energy
  int i_window = (sqrtE - std::sqrt(E_min_)) / spacing_;
  int startw = windows_(i_window, 0) - 1;
  int endw = windows_(i_window, 1) - 1;

  // Initialize the ouptut cross sections.
  double sig_s = 0.0;
  double sig_a = 0.0;
  double sig_f = 0.0;

  // ==========================================================================
  // Add the contribution from the curvefit polynomial.  The math is nasty so
  // the polynomial derivatives are approximated using a finite difference
  // between T + 1K and T - 1K.

  if (broaden_poly_(i_window)) {
    double dopp1 = sqrt_awr_ / std::sqrt(K_BOLTZMANN * (T - 1.0));
    double dopp2 = sqrt_awr_ / std::sqrt(K_BOLTZMANN * (T + 1.0));
    double polys1[fit_order_ + 1];
    double polys2[fit_order_ + 1];
    broaden_wmp_polynomials(E, dopp1, fit_order_+1, polys1);
    broaden_wmp_polynomials(E, dopp2, fit_order_+1, polys2);

    for (int i_poly = 0; i_poly < fit_order_ + 1; ++i_poly) {
      double d_poly = (polys2[i_poly] - polys1[i_poly]) / 2.;
      sig_s += curvefit_(i_window, i_poly, FIT_S) * d_poly;
      sig_a += curvefit_(i_window, i_poly, FIT_A) * d_poly;
      if (fissionable_) {
        sig_f += curvefit_(i_window, i_poly, FIT_F) * d_poly;
      }
    }
  }

  // ==========================================================================
  // Add the contribution from the poles in this window.

  double dopp = sqrt_awr_ / sqrtkT;
  if (endw >= startw) {
    for (int i_pole = startw; i_pole <= endw; ++i_pole) {
      std::complex<double> z = (sqrtE - data_(i_pole, MP_EA)) * dopp;
      std::complex<double> w_val = -invE * SQRT_PI * 0.5 * w_derivative(z, 2);
      sig_s += (data_(i_pole, MP_RS) * w_val).real();
      sig_a += (data_(i_pole, MP_RA) * w_val).real();
      if (fissionable_) {
        sig_f += (data_(i_pole, MP_RF) * w_val).real();
      }
    }
    sig_s *= -0.5*sqrt_awr_ / std::sqrt(K_BOLTZMANN) * std::pow(T, -1.5);
    sig_a *= -0.5*sqrt_awr_ / std::sqrt(K_BOLTZMANN) * std::pow(T, -1.5);
    sig_f *= -0.5*sqrt_awr_ / std::sqrt(K_BOLTZMANN) * std::pow(T, -1.5);
  }

  return std::make_tuple(sig_s, sig_a, sig_f);
}

} // namespace openmc

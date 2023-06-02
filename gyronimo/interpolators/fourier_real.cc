// ::gyronimo:: - gyromotion for the people, by the people -
// An object-oriented library for gyromotion applications in plasma physics.
// Copyright (C) 2021 Paulo Rodrigues.

// ::gyronimo:: is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// ::gyronimo:: is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with ::gyronimo::.  If not, see <https://www.gnu.org/licenses/>.

// @fourier_real.cc, this file is part of ::gyronimo::

#include <numeric>
#include <iostream>
#include <complex>
#include <gyronimo/core/error.hh>
#include <gyronimo/core/dblock.hh>
#include <gyronimo/interpolators/fourier_real.hh>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>

namespace gyronimo {

fourier_real::fourier_real(
    const narray_t& u, const narray_t& dcos, const narray_t& dsin,
    const int mi, const int mf, const interpolator1d_factory* ifactory)
    : m_(mf - mi + 1),
      Acos_(mf - mi + 1, nullptr),
      Asin_(mf - mi + 1, nullptr) {
  std::iota(std::begin(m_), std::end(m_), mi);
  this->build_interpolators(u, dcos, dsin, ifactory);
}
fourier_real::fourier_real(
    const narray_t& u, const narray_t& dcos, const narray_t& dsin,
    const narray_t& _m, const interpolator1d_factory* ifactory)
    : m_(std::begin(_m), std::end(_m)),
      Acos_(_m.size(), nullptr),
      Asin_(_m.size(), nullptr) {
  this->build_interpolators(u, dcos, dsin, ifactory);
}

fourier_real::fourier_real(
      const narray_t& u, const narray_t& v,
      const narray_t& d, const interpolator1d_factory* ifactory)
      : m_(std::floor((v.size()-1)/2)+1), 
        Acos_(std::floor((v.size()-1)/2)+1, nullptr),
        Asin_(std::floor((v.size()-1)/2)+1, nullptr){
  std::iota(m_.begin(), m_.end(), 0.);
  
  narray_t rescos(m_.size()*u.size());
  narray_t ressin(m_.size()*u.size());

  size_t Nsamples = v.size()-1;
  narray_t data(Nsamples);
  double *arr = new double[Nsamples];
  
  gsl_fft_real_wavetable *real;
  gsl_fft_real_workspace *work;

  work = gsl_fft_real_workspace_alloc(Nsamples);
  real = gsl_fft_real_wavetable_alloc(Nsamples);

  for(size_t k = 0; k<u.size(); k++){
    data = d[std::slice(k*v.size(), Nsamples, 1)];

    std::copy(std::begin(data), std::end(data), arr);

    gsl_fft_real_transform(arr, 1, Nsamples, real, work);
    
    std::copy(arr, arr + Nsamples, std::begin(data));

    data *= 2./Nsamples;

    rescos[k] = data[0]*.5;
    ressin[k] = 0.;
    rescos[std::slice(k+u.size(), m_.size()-1, u.size())] = data[std::slice(1,m_.size()-1,2)];
    ressin[std::slice(k+u.size(), m_.size()-2 + Nsamples%2, u.size())] = data[std::slice(2,m_.size()-2 + Nsamples%2,2)];
  }

  ressin *= -1.;

  gsl_fft_real_wavetable_free(real);
  gsl_fft_real_workspace_free(work);

  delete[] arr;

  this->build_interpolators(u, rescos, ressin, ifactory);
}

fourier_real::~fourier_real() {
  for (size_t p = 0; p < m_.size(); p++) {
    if (Acos_[p]) delete Acos_[p];
    if (Asin_[p]) delete Asin_[p];
  }
}
double fourier_real::operator()(double u, double v) const {
  double sum = 0.0;
  for (size_t p = 0; p < m_.size(); p++){
    sum += ((*Acos_[p])(u))*std::cos(m_[p]*v);
    sum += ((*Asin_[p])(u))*std::sin(m_[p]*v);
  }
  return sum;
}
double fourier_real::partial_u(double u, double v) const {
  double sum = 0.0;
  for (size_t p = 0; p < m_.size(); p++){
    sum += ((*Acos_[p]).derivative(u))*std::cos(m_[p]*v);
    sum += ((*Asin_[p]).derivative(u))*std::sin(m_[p]*v);
  }
  return sum;
}
double fourier_real::partial_v(double u, double v) const {
  double sum = 0.0;
  for (size_t p = 0; p < m_.size(); p++){
    sum += -m_[p]*((*Acos_[p])(u))*std::sin(m_[p]*v);
    sum += m_[p]*((*Asin_[p])(u))*std::cos(m_[p]*v);
  }
  return sum;
}
double fourier_real::partial2_uu(double u, double v) const {
  double sum = 0.0;
  for (size_t p = 0; p < m_.size(); p++){
    sum += ((*Acos_[p]).derivative2(u))*std::cos(m_[p]*v);
    sum += ((*Asin_[p]).derivative2(u))*std::sin(m_[p]*v);
  }
  return sum;
}
double fourier_real::partial2_uv(double u, double v) const {
  double sum = 0.0;
  for (size_t p = 0; p < m_.size(); p++){
    sum += -m_[p]*((*Acos_[p]).derivative(u))*std::sin(m_[p]*v);
    sum += m_[p]*((*Asin_[p]).derivative(u))*std::cos(m_[p]*v);
  }
  return sum;
}
double fourier_real::partial2_vv(double u, double v) const {
  double sum = 0.0;
  for (size_t p = 0; p < m_.size(); p++){
    sum += -m_[p]*m_[p]*((*Acos_[p])(u))*std::cos(m_[p]*v);
    sum += -m_[p]*m_[p]*((*Asin_[p])(u))*std::sin(m_[p]*v);
  }
  return sum;
}
void fourier_real::build_interpolators(
    const narray_t& u, const narray_t& dcos, const narray_t& dsin,
    const interpolator1d_factory* ifactory) {
  if (dsin.size() != dcos.size())
    error(__func__, __FILE__, __LINE__, "mismatched dcos and dsin.", 1);
  if (dsin.size() != u.size()*m_.size())
    error(__func__, __FILE__, __LINE__, "mismatched dcos, dsin, or u.", 1);
  for (size_t p = 0; p < m_.size(); p++) {
    std::slice index_range(p*u.size(), u.size(), 1);
    Acos_[p] = ifactory->interpolate_data(
        dblock_adapter(u), dblock_adapter<narray_t>(dcos[index_range]));
    Asin_[p] = ifactory->interpolate_data(
        dblock_adapter(u), dblock_adapter<narray_t>(dsin[index_range]));
  }
}

} // end namespace gyronimo.

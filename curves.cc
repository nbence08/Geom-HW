#include "curves.hh"

#include <algorithm>
#include <cassert>
#include <limits>

double BezierCurve::bernstein(size_t i, size_t n, double u)
{
  DoubleVector tmp(n + 1, 0.0);
  tmp[n-i] = 1.0;
  double u1 = 1.0 - u;
  for (size_t k = 1; k <= n; ++k)
    for (size_t j = n; j >= k; --j)
      tmp[j] = tmp[j] * u1 + tmp[j-1] * u;
  return tmp[n];
}

Point BezierCurve::evaluateOneByOne(double u) const
{
  Point p(0.0, 0.0, 0.0);
  for (size_t k = 0; k <= n; ++k)
    p += cp[k] * bernstein(k, n, u);
  return p;
}

void BezierCurve::bernsteinAll(size_t n, double u, DoubleVector &coeff)
{
  coeff.clear(); coeff.reserve(n + 1);
  coeff.push_back(1.0);
  double u1 = 1.0 - u;
  for (size_t j = 1; j <= n; ++j) {
    double saved = 0.0;
    for (size_t k = 0; k < j; ++k) {
      double tmp = coeff[k];
      coeff[k] = saved + tmp * u1;
      saved = tmp * u;
    }
    coeff.push_back(saved);
  }
}

Point BezierCurve::evaluate(double u) const
{
  DoubleVector coeff; bernsteinAll(n, u, coeff);
  Point p(0.0, 0.0, 0.0);
  for (size_t k = 0; k <= n; ++k)
    p += cp[k] * coeff[k];
  return p;
}

Point BezierCurve::evaluateWithCachedCofficients(const DoubleVector &coeff) const
{
  Point p(0.0, 0.0, 0.0);
  for (size_t k = 0; k <= n; ++k)
    p += cp[k] * coeff[k];
  return p;
}

Point BezierCurve::evaluateByDeCasteljau(double u) const
{
  PointVector tmp = cp;
  double u1 = 1.0 - u;
  for (size_t k = 1; k <= n; ++k)
    for (size_t i = 0; i <= n - k; ++i)
      tmp[i] = tmp[i] * u1 + tmp[i+1] * u;
  return tmp[0];
}

void BezierCurve::derivativeControlPoints(size_t d, PointMatrix &dcp) const
{
  dcp.clear(); dcp.resize(d + 1);
  dcp[0] = cp;
  for (size_t k = 1; k <= d; ++k) {
    size_t tmp = n - k + 1;
    dcp[k].reserve(tmp);
    for (size_t i = 0; i <= n - k; ++i)
      dcp[k].push_back((dcp[k-1][i+1] - dcp[k-1][i]) * tmp);
  }
}

void BezierCurve::bernsteinAll(size_t n, double u, DoubleMatrix &coeff)
{
  coeff.clear(); coeff.resize(n + 1);
  coeff[0].push_back(1.0);
  double u1 = 1.0 - u;
  for (size_t j = 1; j <= n; ++j) {
    coeff[j].reserve(j + 1);
    double saved = 0.0;
    for (size_t k = 0; k < j; ++k) {
      double tmp = coeff[j-1][k];
      coeff[j].push_back(saved + tmp * u1);
      saved = tmp * u;
    }
    coeff[j].push_back(saved);
  }
}

Point BezierCurve::derivativesByControlPoints(double u, size_t d, VectorVector &der) const
{
  size_t du = std::min(d, n);
  der.clear(); der.reserve(d + 1);
  DoubleMatrix coeff; bernsteinAll(n, u, coeff);
  PointMatrix dcp; derivativeControlPoints(du, dcp);
  for (size_t k = 0; k <= du; ++k) {
    der.emplace_back(0.0, 0.0, 0.0);
    for (size_t j = 0; j <= n - k; ++j)
      der[k] += dcp[k][j] * coeff[n-k][j];
  }
  for (size_t k = n + 1; k <= d; ++k)
    der.emplace_back(0.0, 0.0, 0.0);
  return der[0];
}



#include "vector.hh"

#include <cstddef>

using Point        = Vector;
using PointVector  = std::vector<Point>;
using VectorVector = std::vector<Vector>;
using DoubleVector = std::vector<double>;
using DoubleMatrix = std::vector<DoubleVector>;
using PointMatrix  = std::vector<PointVector>;

struct BezierCurve
{
  size_t n;                     // degree => n + 1 control points
  PointVector cp;               // control points

  static double bernstein(size_t i, size_t n, double u);
  Point evaluateOneByOne(double u) const;
  static void bernsteinAll(size_t n, double u, DoubleVector &coeff);
  Point evaluate(double u) const;
  Point evaluateWithCachedCofficients(const DoubleVector &coeff) const;
  Point evaluateByDeCasteljau(double u) const;
  void derivativeControlPoints(size_t d, PointMatrix &dcp) const;
  static void bernsteinAll(size_t n, double u, DoubleMatrix &coeff);
  Point derivativesByControlPoints(double u, size_t d, VectorVector &der) const;
};

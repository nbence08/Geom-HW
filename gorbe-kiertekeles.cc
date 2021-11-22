/******************************************************************************







                    Görbekiértékelési algoritmusok

                        a "NURBS Book" alapján


                    Dr. Várady Tamás, Salvi Péter
               BME Villamosmérnöki és Informatikai Kar
               Irányítástechnika és Informatika Tanszék






******************************************************************************/


// [1] Bézier görbék reprezentációja

// Vector típus a szokásos műveletekkel
#include "vector.hh"

using Point        = Vector;
using PointVector  = std::vector<Point>;
using VectorVector = std::vector<Vector>;
using DoubleVector = std::vector<double>;
using DoubleMatrix = std::vector<DoubleVector>;
using PointMatrix  = std::vector<PointVector>;

struct BezierCurve
{
    size_t n;                     // fokszám => n + 1 kontrollpont
    PointVector cp;               // kontrollpontok

    double bernstein(size_t i, size_t n, double u) const; // member függvények
    Point BezierCurve::evaluateOneByOne(double u) const;
    void BezierCurve::bernsteinAll(size_t n, double u, DoubleVector& coeff) const;
    Point BezierCurve::evaluate(double u) const;
    size_t binomial(size_t n, size_t k);
    Point evaluateWithCachedCofficients(const DoubleVector& coeff) const;
    Point evaluateByDeCasteljau(double u) const;
    void derivativeControlPoints(size_t d, PointMatrix& dcp) const;
    void bernsteinAll(size_t n, double u, DoubleMatrix& coeff);
    Point derivativesByControlPoints(double u, size_t d, VectorVector& der);


};


// [2] Bézier kiértékelés a rekurzív Bernstein képlet alapján

double BezierCurve::bernstein(size_t i, size_t n, double u) const
{
  DoubleVector tmp(n + 1, 0.0);
  tmp[n-i] = 1.0;
  double u1 = 1.0 - u;
  for (size_t k = 1; k <= n; ++k)
    for (size_t j = n; j >= k; --j)
      tmp[j] = tmp[j-1] * u + tmp[j] * u1;
  return tmp[n];
}

Point BezierCurve::evaluateOneByOne(double u) const
{
  Point p(0.0, 0.0, 0.0);
  for (size_t k = 0; k <= n; ++k)
    p += cp[k] * bernstein(k, n, u);
  return p;
}


// [3] Az összes Bernstein polinom kiértékelése egyszerre

void BezierCurve::bernsteinAll(size_t n, double u, DoubleVector &coeff) const
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

// [4] Görbekiértékelés előre kiszámolt Bernstein polinomokkal

Point BezierCurve::evaluate(double u) const
{
  DoubleVector coeff; bernsteinAll(n, u, coeff);
  Point p(0.0, 0.0, 0.0);
  for (size_t k = 0; k <= n; ++k)
    p += cp[k] * coeff[k];
  return p;
}

Point BezierCurve::
evaluateWithCachedCofficients(const DoubleVector &coeff) const
{
  Point p(0.0, 0.0, 0.0);
  for (size_t k = 0; k <= n; ++k)
    p += cp[k] * coeff[k];
  return p;
}


// [5] deCasteljau algoritmus

Point BezierCurve::evaluateByDeCasteljau(double u) const
{
  PointVector tmp = cp;
  double u1 = 1.0 - u;
  for (size_t k = 1; k <= n; ++k)
    for (size_t i = 0; i <= n - k; ++i)
      tmp[i] = tmp[i] * u1 + tmp[i+1] * u;
  return tmp[0];
}


// [6] Bézier deriváltak kontrollpontjainak kiszámítása

// Feltételezi, hogy d <= n
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


// [7] Az összes Bernstein polinom (minden kisebb fokszámra is)

void BezierCurve::bernsteinAll(size_t n, double u, DoubleMatrix &coeff)
{
  coeff.clear(); coeff.resize(n + 1);
  coeff[0].push_back(1.0);
  double u1 = 1.0 - u;
  for (size_t j = 1; j <= n; ++j) {
    coeff[j].reserve(j + 1);
    double saved = 0.0;
    for (size_t k = 0; k < j; ++k) {
      double tmp = coeff[j-1][k];     // ... = coeff[k]  helyett
      coeff[j].push_back(saved + tmp * u1); // coeff[k] = ...  helyett
      saved = tmp * u;
    }
    coeff[j].push_back(saved);
  }
}


// [8] Deriváltszámítás

Point BezierCurve::
derivativesByControlPoints(double u, size_t d, VectorVector &der)
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

// Vége

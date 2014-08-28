/******************************************************************************
 * libKrylov
 * =========
 * Math.hxx contians the core mathematical data structures used in libKrylov
 *
 * 19 August 2014
 * ~ ry
 * ***************************************************************************/
#ifndef KRY_MATH_HXX
#define KRY_MATH_HXX

#include "kry/Utility.hxx"
#include "kry/Runtime.hxx"
#include <initializer_list>
#include <memory>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <random>

namespace kry
{

class Vector;
class Matrix;
class Column;
class ColumnRange;
class SparseMatrix;

class Rotator;

template<class A> A T(A a) 
{ 
  a.transpose();
  return a;
}

class Vector
{
  public:
    explicit Vector(size_t n);
    Vector(std::initializer_list<double>);
    Vector(Column c);
    static Vector Zero(size_t n);
    static Vector Random(size_t n, double mean, double variance);
    ~Vector();

    double & operator()(size_t i) const;
    Vector operator!() const;

    Vector operator=(const Column x);

    Vector operator+=(Vector x) const;
    Vector operator-=(Vector x) const;
    Vector operator*=(double x) const;
    Vector operator/=(double x) const;

    size_t n() const;

  private:
    size_t _n;
    std::shared_ptr<double> _data;
};

std::ostream & operator<<(std::ostream&, const Vector);

Vector operator+(Vector a, Vector b);
Vector operator-(Vector a, Vector b);
double operator*(Vector a, Vector b);
Vector operator*(Vector a, double b);
Vector operator*(double a, Vector b);
Vector operator/(Vector a, double b);
bool operator==(const Vector a, const Vector b);

double norm(const Vector a);

class Matrix
{
  public:
    Matrix(size_t m, size_t n);
    Matrix(size_t m, size_t n, std::initializer_list<double>);
    static Matrix Zero(size_t m, size_t n);
    static Matrix Identity(size_t m, size_t n);
    static Matrix Random(size_t m, size_t n, double mean, double variance);

    double & operator()(size_t i, size_t j) const;
    Matrix operator!() const;

    Matrix operator*=(const Matrix A);

    Column C(size_t idx);
    ColumnRange C(size_t begin, size_t end);

    size_t m() const, n() const;

    void transpose();
    bool transposed{false};

  private:
    size_t _m, _n;
    std::shared_ptr<double> _data;
};

bool operator==(Matrix A, Matrix B);
Vector operator*(Matrix A, Vector x);
Matrix operator*(Matrix A, Matrix B);

Vector back_substitute(Matrix, Vector);

std::ostream & operator<<(std::ostream &o, const Matrix x);

class Column
{
  public:
    Column(Matrix M, size_t idx);

    size_t idx() const;
    Matrix M() const;

    Column operator=(const Vector);
    Column operator-=(const Vector);
    Column operator+=(const Vector);
    Column operator()(size_t rbegin, size_t rend);
    double & operator()(size_t idx);

    size_t rbegin, rend;

  private:
    size_t _idx;
    Matrix _M;
};

bool operator==(const Column, const Vector);

class ColumnRange
{
  public:
    ColumnRange(Matrix M, size_t begin, size_t end);

    size_t begin() const, end() const;
    Matrix M() const;

    void transpose();

  private:
    size_t _begin, _end;
    Matrix _M;
};

Vector operator*(const ColumnRange, const Vector);
Vector operator*(const ColumnRange, const Column);


class SparseMatrix
{

  public:
    SparseMatrix(size_t m, size_t n, size_t z);
    SparseMatrix(size_t m, size_t n, size_t z, 
        std::vector<size_t> rs, std::vector<size_t> cs, 
        std::vector<double> vs);
    static SparseMatrix Identity(size_t m, size_t n, size_t z);

    double & operator()(size_t i, size_t j);
    double get(size_t i, size_t j);

    size_t m() const, n() const, z() const;
    size_t r(size_t i) const;
    size_t c(size_t i, size_t j) const, c(size_t k) const;
    double v(size_t i, size_t j) const, v(size_t k) const;

  private:
    size_t _m, _n, _z;
    std::shared_ptr<size_t> _r, _c;
    std::shared_ptr<double> _v;
};

Vector operator*(const SparseMatrix, const Vector);
Vector operator*(const SparseMatrix, const Column);
std::ostream & operator<<(std::ostream &, SparseMatrix);


class Rotator
{
  public:
    Rotator(Vector, size_t, size_t);
    Rotator(Matrix, size_t, size_t);

    Vector apply(Vector);
    Matrix apply_left(Matrix);
    Matrix apply_right(Matrix);

    size_t i() const, j() const;
    double c() const, s() const;

  private:
    size_t _i, _j;
    double xi, xj;
    double _c, _s;

    void compute_c_s();
};

Vector operator* (const Rotator, const Vector);
Matrix operator* (const Rotator, const Matrix);
Matrix operator* (const Matrix, const Rotator);

}

#endif
 

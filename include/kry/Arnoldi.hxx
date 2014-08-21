/******************************************************************************
 * libKrylov
 * =========
 * Arnoldi.hxx contians the data structures associated with the Arnoldi 
 * algorithm
 *
 * 19 August 2014
 * ~ ry
 * ***************************************************************************/
#ifndef KRY_ARNOLDI_HXX
#define KRY_ARNOLDI_HXX

#include "kry/Math.hxx"
#include "kry/Utility.hxx"

namespace kry 
{

class Arnoldi
{
  public:
    Arnoldi() = delete;
    Arnoldi(size_t n, SparseMatrix A, Vector x0, Vector b);

    Matrix Q, H;
    SparseMatrix A;
    Vector x0, b, r0, d, t, xn;
    double r0_norm;

    size_t N() const, n() const;
    void operator()();

  private:
    size_t _n;
    void rotate_h2t();

};

}

#endif

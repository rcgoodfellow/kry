/******************************************************************************
 * libKrylov
 * =========
 * Math.cxx contians the implementation of the  core mathematical data 
 * structures used in libKrylov
 *
 * 19 August 2014
 * ~ ry
 * ***************************************************************************/

#include "kry/Arnoldi.hxx"

using namespace kry;

Arnoldi::Arnoldi(size_t n, SparseMatrix A, Vector x0, Vector b)
  : Q{Matrix::Zero(A.m(),n)}, 
    H{Matrix::Zero(n,n)}, 
    A{A}, x0{x0}, b{b}, r0(A.m()),
    _n{n}
{
  if(A.n() != x0.n() || A.n() != b.n()) { throw NON_CONFORMAL; }

}

size_t Arnoldi::N() const { return A.m(); }
size_t Arnoldi::n() const { return _n; }

void Arnoldi::operator()()
{
  r0 = b - A*x0;  
  r0 /= norm(r0);
  Q.C(0) = r0;

  for(size_t i=0; i<n()-1; ++i)
  {
    //Next subspace element
    Q.C(i+1) = A*Q.C(i);

    //ortho
    H.C(i) = T(Q.C(0,i+1)) * Q.C(i+1);
    Q.C(i+1) -= Q.C(0,i+1) * H.C(i)(0,i+1);

    //reortho
    Vector s = T(Q.C(0,i+1)) * Q.C(i+1);
    Q.C(i+1) -= Q.C(0,i+1) * s;
    H.C(i) += s;

    //normo
    Vector x = Q.C(i+1);
    double xn = norm(x);
    H.C(i)(i+1) = xn;
    Q.C(i+1) = x / xn;
  }
}

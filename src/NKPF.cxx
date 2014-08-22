/******************************************************************************
 * libKrylov
 * =========
 * NKPF.cxx
 *
 * This file contains the implementation of the Newton-Krylov power flow
 *
 * 21 August 2014
 * ~ ry
 * ***************************************************************************/
#include "kry/NKPF.hxx"

using namespace kry;

size_t JacobiMap::size() { return j0_sz + j1_sz; }

NKPF::NKPF(SparseMatrix Y, SparseMatrix YA, JacobiMap jmap, size_t n)
  : 
    Q(jmap.size(), n),
    H(n,n),
    ve(Y.n()*2),
    dve(jmap.size()),
    ps(Y.n()*2),
    pc(Y.n()*2),
    dp(jmap.size()),
    dv0(jmap.size()),
    dr0(jmap.size()),
    qdp(n),
    qdv(n),
    Y(Y), YA(YA), jmap(jmap)
{ }

double NKPF::g(size_t i) { return Y(i,i)*cos(YA(i,i)); }
double NKPF::b(size_t i) { return Y(i,i)*sin(YA(i,i)); }

double NKPF::p(size_t i)
{
  double pi = pow(v(i),2)*g(i);
  for(size_t k=0; k<Y.r(i); ++k)
  {
    size_t j = Y.c(i,k);
    pi += Y(i,j) * v(i) * v(j) * cos( YA(i,j) + va(j) - va(i) );
  }
  return pi;
}

double NKPF::q(size_t i)
{
  double qi = -pow(v(i),2)*b(i);
  for(size_t k=0; k<Y.r(i); ++k)
  {
    size_t j = Y.c(i,k);
    qi -= Y(i,j) * v(i) * v(j) * sin( YA(i,j) + va(j) - va(i) );
  }
  return qi;
}

//real-power gradient ---------------------------------------------------------
double NKPF::jdp(size_t i) { return jdp_va(i) + jdp_v(i); }

double NKPF::jdp_va(size_t i)
{
  double ods{0}, s{0};
  for(size_t k=0; k<Y.r(i); ++k)
  {
    size_t j = Y.c(i, k);
    ods += dp_dva(i, j); //implicitly compute diagonal coefficient
    s += ods * dva(j);
  }
  s += -ods * dva(i);
  return s;
}

double NKPF::jdp_v(size_t i)
{
  double ods{0}, s{0};
  for(size_t k=0; k<Y.r(i); ++k)
  {
    size_t j = Y.c(i,k);
    ods += dp_dv(i, j);
    s += ods * dv(j);
  }
  s += (ods + 2*pow(v(i),2)*g(i)) * dv(i);
  return s;
}

double NKPF::dp_dva(size_t i, size_t j)
{
  return -v(i) * v(j) * Y(i,j) * sin(YA(i,j) + va(j) - va(i));
}

double NKPF::dp_dv(size_t i, size_t j)
{
  return v(i) * v(j) * Y(i,j) * cos(YA(i,j) + va(j) - va(i));
}

//reactive-power gradient -----------------------------------------------------
double NKPF::jdq(size_t i) { return jdq_va(i) + jdq_v(i); }

double NKPF::jdq_va(size_t i)
{
  double ods{0}, s{0};
  for(size_t k=0; k<Y.r(i); ++k)
  {
    size_t j = Y.c(i,k);
    ods += dq_dva(i,j);
    s += ods * dva(j);
  }
  s += -ods * dva(i);
  return s;
}

double NKPF::jdq_v(size_t i)
{
  double ods{0}, s{0};
  for(size_t k=0; k<Y.r(i); ++k)
  {
    size_t j = Y.c(i,k);
    ods += dq_dv(i,j);
    s += ods * dv(j);
  }
  s += (-ods - 2*pow(v(i),2) * b(i)) * dv(i); 
  return s;
}

double NKPF::dq_dva(size_t i, size_t j)
{
  return -v(i) * v(j) * Y(i,j) * cos(YA(i,j) + va(j) - va(i));
}

double NKPF::dq_dv(size_t i, size_t j)
{
  return -v(j) * v(i) * Y(i,j) * sin(YA(i,j) + va(j) - va(i));
}

double NKPF::v(size_t i) { return ve( jmap.map[i].j0 ); }
double NKPF::va(size_t i) { return ve( jmap.map[i].j1 ); }
double NKPF::dv(size_t i) { return dve( jmap.map[i].j0 ); }
double NKPF::dva(size_t i) { return dve( jmap.map[i].j1 ); }

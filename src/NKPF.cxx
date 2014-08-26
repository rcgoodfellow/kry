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

jidx::jidx(int j0, int j1) : j0{j0}, j1{j1} {}

size_t JacobiMap::size() { return j0_sz + j1_sz; }

NKPF::NKPF(SparseMatrix Y, SparseMatrix YA, JacobiMap jmap, size_t n,
    Vector initial, Vector ps)
  : 
    N(Y.m()),
    n(n),
    Q(jmap.size(), n),
    H(n,n),
    ve(initial),
    dve(jmap.size()),
    ps(ps),
    pc(N*2),
    dp(jmap.size()),
    dv0(jmap.size()),
    dr0(jmap.size()),
    qdp(n),
    qdv(n),
    Jdv(jmap.size()),
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
double NKPF::jdp(size_t i) 
{ 
  double val = jdp_va(i) + jdp_v(i); 

  if(val != val) throw std::runtime_error("unreal jacobi response");

  return val;
}

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
  double val = -v(i) * v(j) * Y(i,j) * sin(YA(i,j) + va(j) - va(i));
  return val;
}

double NKPF::dp_dv(size_t i, size_t j)
{
  double val = v(i) * v(j) * Y(i,j) * cos(YA(i,j) + va(j) - va(i));
  return val;
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

double NKPF::v(size_t i) { return ve(i); }
double NKPF::va(size_t i) { return ve(i); }
double NKPF::dv(size_t i) 
{ 
  int idx = jmap.map[i].j0;
  if(idx == -1) return 0;
  return dve( idx ); 
}

double NKPF::dva(size_t i) 
{ 
  int idx = jmap.map[i].j1;
  if(idx == -1) return 0;
  return dve( idx ); 
}

void NKPF::compute_pc()
{
  for(size_t i=0, ii=N; i<pc.n()/2; ++i, ++ii) 
  { 
    pc(i) = p(i); 
    pc(ii) = q(i);
  }
}

void NKPF::compute_dp()
{
  size_t jn = jmap.j0_sz;
  for(size_t i=0; i<N; ++i)
  {
    jidx j = jmap.map[i];
    if(j.j0 != -1) { dp(j.j0) = ps(i) - pc(i); }
    if(j.j1 != -1) { dp(jn + j.j1) = ps(N + i) - pc(N + i); }
  }
}

void NKPF::compute_dve()
{
  dve = Vector::Zero(dve.n());
  for(size_t i=0; i<dve.n(); ++i) { dve(i) = 1; }
  compute_Jdv();  
  std::cout << "Jdv:" << Jdv << std::endl;
}

void NKPF::compute_Jdv()
{
  Jdv = Vector::Zero(Jdv.n());
  for(size_t i=0; i<N; ++i) 
  { 
    int idx = jmap.map[i].j0;
    if(idx != -1) Jdv(idx) = jdp(i); 

    idx = jmap.j0_sz + jmap.map[i].j1;
    if(idx != -1) Jdv(idx) = jdq(i);
  }
}

//compute ---------------------------------------------------------------------
NKPF & NKPF::operator()()
{
  std::cout << "ps:" << ps << std::endl;
  compute_pc();
  std::cout << "pc:" << pc << std::endl;
  compute_dp();
  std::cout << "dp:" << dp << std::endl;
  compute_dve();

  return *this;
}

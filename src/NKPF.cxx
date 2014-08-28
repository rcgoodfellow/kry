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
    Q(Matrix::Zero(jmap.size(), n+1)),
    H(Matrix::Zero(n,n)),
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
    Y(Y), YA(YA), 
    J(jmap.size(), jmap.size(), jmap.size()),
    jmap(jmap)
{ }

double NKPF::g(size_t i) { return Y(i,i)*cos(YA(i,i)); }
double NKPF::b(size_t i) { return Y(i,i)*sin(YA(i,i)); }

double NKPF::p(size_t i)
{
  double pi = pow(v(i),2)*g(i);
  for(size_t k=0; k<Y.r(i); ++k)
  {
    size_t j = Y.c(i,k);
    if(i == j) continue;
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
    if(i == j) continue;
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
  double c{0}, ods{0}, s{0};
  for(size_t k=0; k<Y.r(i); ++k)
  {
    size_t j = Y.c(i, k);
    if(j == i) continue;
    c = dp_dva(i, j); 
    ods += c; //implicitly compute diagonal coefficient
    s += c * dva(j);
  }
  s += -ods * dva(i);
  std::cout << "J11(" << i << ") = " << ods << std::endl;
  return s;
}

double NKPF::jdp_v(size_t i)
{
  double c{0}, ods{0}, s{0};
  for(size_t k=0; k<Y.r(i); ++k)
  {
    size_t j = Y.c(i,k);
    if(j == i) continue;
    c = dp_dv(i, j);
    ods += c;
    s += c * dv(j);
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
  double c{0}, ods{0}, s{0};
  for(size_t k=0; k<Y.r(i); ++k)
  {
    size_t j = Y.c(i,k);
    if(j == i) continue;
    c = dq_dva(i,j);
    ods += c;
    s += c * dva(j);
  }
  s += -ods * dva(i);
  return s;
}

double NKPF::jdq_v(size_t i)
{
  double c{0}, ods{0}, s{0};
  for(size_t k=0; k<Y.r(i); ++k)
  {
    size_t j = Y.c(i,k);
    if(j == i) continue;
    c = dq_dv(i,j);
    ods += c;
    s += c * dv(j);
  }
  double coeff = (-ods - 2*pow(v(i),2) * b(i));
  s +=  coeff * dv(i); 
  std::cout << "J22(" << i << ") = " << coeff << std::endl;
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

double NKPF::v(size_t i) { return ve(i+N); }
double NKPF::va(size_t i) { return ve(i); }
double NKPF::dv(size_t i) 
{ 
  int idx = jmap.map[i].j1;
  if(idx == -1) return 0;
  return dve( jmap.j0_sz + idx ); 
}

double NKPF::dva(size_t i) 
{ 
  int idx = jmap.map[i].j0;
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

void NKPF::j11()
{
  for(size_t i=0; i<N; ++i)
  {
    int _i = jmap.map[i].j0;
    if(_i == -1) continue;
   
    double ii{0};
    for(size_t k=0; k<Y.r(i); ++k)
    {
      size_t j = Y.c(i,k);
      if(i == j) continue;
      double t = dp_dva(i,j);
      ii += t;
      int _j = jmap.map[j].j0;
      if(_j == -1) continue;

      J(_i,_j) = t;
    }
    J(_i,_i) = -ii;
  }
}

void NKPF::j22()
{
  for(size_t i=0; i<N; ++i)
  {
    int _i = jmap.map[i].j1;
    if(_i == -1) continue;

    double ii{0};
    for(size_t k=0; k<Y.r(i); ++k)
    {
      size_t j = Y.c(i,k);
      if(i == j) continue;
      double t = dq_dv(i,j);
      ii += t;
      int _j = jmap.map[j].j1;
      if(_j == -1) continue;

      J(jmap.j0_sz + _i, jmap.j0_sz + _j) = t;
    }
    size_t iidx = jmap.j0_sz + _i;
    J(iidx, iidx) = ii - 2*pow(v(i),2)*b(i);
  }
}

void NKPF::j21()
{
  for(size_t i=0; i<N; ++i)
  {
    int _i = jmap.map[i].j1;
    if(_i == -1) continue;

    double ii{0};
    for(size_t k=0; k<Y.r(i); ++k)
    {
      size_t j = Y.c(i,k);
      if(i == j) continue;
      double t = dq_dva(i,j);
      ii += t;
      int _j = jmap.map[j].j0;
      if(_j == -1) continue;

      J(jmap.j0_sz + _i, _j) = t;
    }
    
    int _ii = jmap.map[i].j0;
    if(_ii != -1) J(jmap.j0_sz + _i, _ii) = -ii;
  }
}

void NKPF::j12()
{
  for(size_t i=0; i<N; ++i)
  {
    int _i = jmap.map[i].j0;
    if(_i == -1) continue;

    double ii{0};
    for(size_t k=0; k<Y.r(i); ++k)
    {
      size_t j = Y.c(i,k);
      if(i == j) continue;
      double t = dp_dv(i,j);
      ii += t;
      int _j = jmap.map[j].j1;
      if(_j == -1) continue;

      J(_i, jmap.j0_sz + _j) = t;
    }

    int _ii = jmap.map[i].j1;
    if(_ii != -1) J(_i, jmap.j0_sz + _ii) = ii + 2*pow(v(i),2)*g(i);
  }
}

void NKPF::build_Jacobi()
{
  j11();
  j22();
  j21();
  j12();
}

using std::cout;
using std::endl;

void NKPF::compute_dve()
{
  std::cout << "pc: " << pc << std::endl;
  dve = Vector::Zero(dve.n());
  build_Jacobi();

  dr0 = dp - J*dve;
  dr0_norm = norm(dr0);
  dr0 /= dr0_norm;
  Q.C(0) = dr0;

  for(size_t i=0; i<jmap.size(); ++i)
  {
    LOG_FUNC();
    Q.C(i+1) = J*Q.C(i);

    H.C(i) = T(Q.C(0,i+1)) * Q.C(i+1);
    Q.C(i+1) -= Q.C(0,i+1) * H.C(i)(0,i+1);
    
    //Vector s = T(Q.C(0,i+1)) * Q.C(i+1);
    //Q.C(i+1) -= Q.C(0,i+1) * s;
    //H.C(i) += s;

    Vector x = Q.C(i+1);
    double xn = norm(x);
    H.C(i)(i+1) = xn;
    LOG_KVP(i);
    LOG_KVP(xn);
    //dve = x / xn;
    Q.C(i+1) = x / xn;

  }

  std::fstream fs ("Q.matrix", std::fstream::out);
  fs << Q;
  fs.close();

  fs.open("H.matrix", std::fstream::out);
  fs << H;
  fs.close();

  qdp = Vector::Zero(qdp.n());
  qdp(0) = dr0_norm;

  for(size_t i=0; i<n-1; ++i)
  {
    Rotator r(H, i, i+1);
    r.apply_left(H);
    r.apply(qdp);
  }

  fs.open("Ht.matrix", std::fstream::out);
  fs << H;
  fs.close();

  qdv = back_substitute(H, qdp);
  std::cout << "qdv:" << qdv << endl; 

  dve = Q.C(0,n) * qdv;
  std::cout << "dve:" << dve << endl;

  std::cout << "dp:" << dp << endl;


  //std::cout << "Q:" << endl << Q << std::endl;
  //std::cout << "H:" << endl << H << std::endl;

  //std::cout << J << std::endl;

  /*
  compute_Jdv();  
  std::cout << "Jdv:" << Jdv << std::endl;

  dr0 = dp - Jdv;
  std::cout << "dr0:" << dr0 << std::endl;

  dr0_norm = norm(dr0);
  std::cout << "dr0_norm:" << dr0_norm << std::endl;

  dr0 /= dr0_norm;
  std::cout << "dr0:" << dr0 << std::endl;

  dve = dr0;

  Q.C(0) = dr0;
  std::cout << "Q:" << Q << std::endl;

  for(size_t i=0; i<n; ++i)
  {
    LOG_FUNC();
    compute_Jdv();
    Q.C(i+1) = Jdv;
    std::cout << "Jdv:" << Jdv << std::endl;

    H.C(i) = T(Q.C(0,i+1)) * Q.C(i+1);
    Q.C(i+1) -= Q.C(0,i+1) * H.C(i)(0,i+1);
    
    //Vector s = T(Q.C(0,i+1)) * Q.C(i+1);
    //Q.C(i+1) -= Q.C(0,i+1) * s;
    //H.C(i) += s;

    Vector x = Q.C(i+1);
    double xn = norm(x);
    H.C(i)(i+1) = xn;
    LOG_KVP(i);
    LOG_KVP(xn);
    dve = x / xn;
    Q.C(i+1) = dve;
  }
  std::cout << "Q:" << Q << std::endl;
  std::cout << "H:" << H << std::endl;
  */
}

void NKPF::compute_Jdv()
{
  Jdv = Vector::Zero(Jdv.n());
  for(size_t i=0; i<N; ++i) 
  { 
    int idx = jmap.map[i].j0;
    if(idx != -1) Jdv(idx) = jdp(i); 

    idx = jmap.map[i].j1;
    if(idx != -1) Jdv(jmap.j0_sz + idx) = jdq(i);
  }
}

//compute ---------------------------------------------------------------------
NKPF & NKPF::operator()()
{
  //std::cout << "ps:" << ps << std::endl;
  compute_pc();
  //std::cout << "pc:" << pc << std::endl;
  compute_dp();
  //std::cout << "dp:" << dp << std::endl;
  compute_dve();

  return *this;
}

/******************************************************************************
 * libKrylov
 * =========
 * Math.cxx contians the implementation of the  core mathematical data 
 * structures used in libKrylov
 *
 * 19 August 2014
 * ~ ry
 * ***************************************************************************/

#include "kry/Math.hxx"

using namespace kry;

using std::runtime_error;
using std::min;
using std::max;
using std::ceil;
using std::ostream;
using std::setprecision;
using std::fixed;

Vector::Vector(size_t n)
  : _n{n}, _data(alloc<double>(n), _mm_free)
{ }

Vector::Vector(std::initializer_list<double> lst)
  : Vector(lst.size())
{
  size_t i{0};
  double *_d = _data.get();
  for(double d : lst) { _d[i++] = d; };
}

Vector Vector::Zero(size_t n)
{
  Vector x(n);
  for(size_t i=0; i<n; ++i){ x(i) = 0; }
  return x;
}

Vector Vector::Random(size_t n, double mean, double variance)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> d(mean, variance);

  Vector x(n);
  for(size_t i=0; i<n; ++i){ x(i) = d(gen); }
  return x;
}

Vector::~Vector()
{ }

double & Vector::operator()(size_t i) const
{
  return _data.get()[i];
}

Vector Vector::operator!() const
{
  Vector x(n());
  for(size_t i=0; i<n(); ++i){ x(i) = this->operator()(i); }
  return x;
}
    
Vector Vector::operator+=(Vector x) const
{
  for(size_t i=0; i<n(); ++i){ this->operator()(i) += x(i); }
  return *this;
}

Vector Vector::operator-=(Vector x) const
{
  for(size_t i=0; i<n(); ++i){ this->operator()(i) -= x(i); }
  return *this;
}
    
Vector Vector::operator*=(double x) const
{

  for(size_t i=0; i<n(); ++i){ this->operator()(i) *= x; }
  return *this;
}

Vector Vector::operator/=(double x) const
{

  for(size_t i=0; i<n(); ++i){ this->operator()(i) /= x; }
  return *this;
}

size_t Vector::n() const
{
  return _n;
}

Vector kry::operator+(Vector a, Vector b)
{
  if(a.n() != b.n()) { throw NON_CONFORMAL; }

  Vector c(a.n());
  for(size_t i=0; i<a.n(); ++i) { c(i) = a(i) + b(i); }

  return c;
}

Vector kry::operator-(Vector a, Vector b)
{
  if(a.n() != b.n()) { throw NON_CONFORMAL; }

  Vector c(a.n());
  for(size_t i=0; i<a.n(); ++i) { c(i) = a(i) - b(i); }

  return c;
}

double kry::operator*(Vector a, Vector b)
{
  if(a.n() != b.n()) { throw NON_CONFORMAL; }
  
  double d{0};
  for(size_t i=0; i<a.n(); ++i) { d += a(i) * b(i); }

  return d;
}

bool kry::operator==(const Vector a, const Vector b)
{
  if(a.n() != b.n()) { throw NON_CONFORMAL; }

  for(size_t i=0; i<a.n(); ++i) { if(a(i) != b(i)) return false; }
  return true;
}

Vector kry::operator*(Vector a, double b)
{
  Vector c(a.n());
  for(size_t i=0; i<a.n(); ++i){ c(i) = a(i) * b; }
  return c;
}

Vector kry::operator*(double a, Vector b)
{
  return b * a;
}

Vector kry::operator/(Vector a, double b)
{
  Vector c(a.n());
  for(size_t i=0; i<a.n(); ++i){ c(i) = a(i) / b; }
  return c;
}

double kry::norm(const Vector a)
{
  return sqrt(a * a);
}

//Matrix ----------------------------------------------------------------------

Matrix::Matrix(size_t m, size_t n)
  : _m{m}, _n{n}, _data(alloc<double>(m*n), _mm_free)
{}

Matrix::Matrix(size_t m, size_t n, std::initializer_list<double> lst)
  : Matrix(m, n)
{
  if(lst.size() != m*n) { throw NON_CONFORMAL; }

  double *_d = _data.get();
  size_t i{0};
  for(double d : lst) { _d[i++] = d; }
}

Matrix Matrix::Zero(size_t m, size_t n)
{
  Matrix A(m,n);
  for(size_t i=0; i<m; ++i){ for(size_t j=0; j<n; ++j) 
  {
    A(i,j) = 0;
  }}
  return A;
}

Matrix Matrix::Identity(size_t m, size_t n)
{
  Matrix A = Matrix::Zero(m,n);
  size_t last = min(m,n);
  for(size_t i=0; i<last; ++i){ A(i,i) = 1; }
  return A;
}
    
Matrix Matrix::Random(size_t m, size_t n, double mean, double variance)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> d(mean, variance);

  Matrix A(m, n);
  for(size_t i=0; i<A.m(); ++i) {
    for(size_t j=0; j<A.m(); ++j) { A(i,j) = d(gen); }
  }
  return A;
}

size_t Matrix::m() const { return _m; }
size_t Matrix::n() const { return _n; }

void Matrix::transpose() { transposed = true; }
    
double & Matrix::operator()(size_t i, size_t j) const
{
  return _data.get()[i*n() + j];
}

Matrix Matrix::operator!() const
{
  Matrix A(m(),n());
  for(size_t i=0; i<m(); ++i)
  {
    for(size_t j=0; j<n(); ++j)
    {
      A(i,j) = (*this)(i,j);
    }
  }
  return A;
}

Matrix Matrix::operator*=(const Matrix A)
{
  Matrix T = !(*this);
  *this = T * A;
  return *this;
}

Column Matrix::C(size_t idx)
{
  if(idx >= n()) { throw BAD_COLREF }
  return Column(*this, idx);
}

ColumnRange Matrix::C(size_t begin, size_t end)
{
  return ColumnRange(*this, begin, end);
}

bool kry::operator==(Matrix A, Matrix B)
{
  if(A.n() != B.n() || A.m() != B.m()) return false;

  for(size_t i=0; i<A.m(); ++i)
  {
    for(size_t j=0; j<A.n(); ++j)
    {
      if(A(i,j) != B(i,j)) return false;
    }
  }
  return true;
}

Vector mul_mv_t(Matrix A, Vector x, size_t cbegin, size_t cend)
{
  size_t n = cend - cbegin ;
  if(cbegin > cend || cend > A.n()) { throw BAD_COLRANGE; }
  if(x.n() != A.m()) { throw NON_CONFORMAL; }

  LOG_FUNC();
  Vector Ax = Vector::Zero(n);

  size_t 
    dpt = max((size_t)(A.m()*n/(double)RT::thread_count()), 
              RT::MEPT),                         //data per thd
    rpt = ceil(dpt/(double)A.m()),               //row per thd
    nt = ceil(n/(double)rpt);                    //# of thds

  LOG_KVP(dpt);
  LOG_KVP(rpt);
  LOG_KVP(nt);

  CountdownLatch cl(nt);
  for(size_t tid=0; tid<nt; ++tid)
  {
    Task t = [tid,A,x,Ax,rpt,&cl,cbegin,cend]()
    {
      size_t jbegin = tid*rpt + cbegin,
             jend = jbegin + rpt;
      
      jend = min(jend, cend); //don't overshoot

      for(size_t i=0; i<A.m(); ++i)
      {
        //the hope is that the compiler vectorizes this loop
        for(size_t j=jbegin; j<jend; ++j)
        {
          Ax(j) += A(i, j) * x(i);
        }
      }
      --cl;
    };
    RT::get().workers[tid]->mtx->lock();
    RT::get().workers[tid]->task = t;
    RT::get().workers[tid]->mtx->unlock();
    RT::get().workers[tid]->cnd->notify_all();
  }
  cl.wait();

  return Ax;
}

Vector mul_mv(Matrix A, Vector x, size_t cbegin, size_t cend)
{
  size_t n = cend - cbegin ;
  if(cbegin > cend || cend > A.n()) { throw BAD_COLRANGE; }
  if(n != x.n()) { throw NON_CONFORMAL; }

  LOG_FUNC();
  Vector Ax = Vector::Zero(A.m());

  size_t 
    dpt = max((size_t)ceil(A.m()*n/(double)RT::thread_count()), 
              RT::MEPT),                             //data per thd
    rpt = ceil(dpt/(double)n),                       //row per thd
    nt = ceil(A.m()/(double)rpt);                    //# of thds

  LOG_KVP(dpt);
  LOG_KVP(rpt);
  LOG_KVP(nt);
  LOG_KVP(cbegin);
  LOG_KVP(cend);

  CountdownLatch cl(nt);
  for(size_t tid=0; tid<nt; ++tid)
  {
    Task t = [tid,A,x,Ax,rpt,&cl,cbegin,cend]()
    {
      size_t jbegin = tid*rpt,
             jend = jbegin + rpt;
      
      jend = min(jend, A.m()); //don't overshoot

      for(size_t i=cbegin; i<cend; ++i)
      {
        //the hope is that the compiler vectorizes this loop
        for(size_t j=jbegin; j<jend; ++j)
        {
          Ax(j) += A(j, i) * x(i);
        }
      }
      --cl;
    };
    RT::get().workers[tid]->mtx->lock();
    RT::get().workers[tid]->task = t;
    RT::get().workers[tid]->mtx->unlock();
    RT::get().workers[tid]->cnd->notify_all();
  }
  cl.wait();

  return Ax;
}

Vector kry::operator*(Matrix A, Vector x)
{
  if(A.transposed){ return mul_mv_t(A, x, 0, A.n()); }
  else { return mul_mv(A, x, 0, A.n()); }
}

Matrix kry::operator*(Matrix A, Matrix B)
{
  if(A.n() != B.m()) { throw NON_CONFORMAL; }
  LOG_FUNC();
  Matrix AB = Matrix::Zero(A.m(), B.n());

  size_t
    rpt = ceil(AB.m()/(double)RT::thread_count()),
    nt = min(RT::thread_count(), A.m());

  LOG_KVP(rpt);
  LOG_KVP(nt);

  CountdownLatch cl(nt);
  for(size_t tid=0; tid<nt; ++tid)
  {
    Task t = [tid,A,B,AB,rpt,&cl]()
    {
      size_t ibegin = tid*rpt,
             iend = ibegin + rpt;

      iend = min(iend, AB.m());

      for(size_t i=ibegin; i<iend; ++i)
      {
        for(size_t j=0; j<AB.n(); ++j)
        {
          for(size_t k=0; k<A.n(); ++k)
          {
            AB(i,j) += A(i,k) * B(k,j);
          }
        }
      }
      --cl;
    };
    RT::get().workers[tid]->mtx->lock();
    RT::get().workers[tid]->task = t;
    RT::get().workers[tid]->mtx->unlock();
    RT::get().workers[tid]->cnd->notify_all();
  }
  cl.wait();

  return AB;
}

Vector kry::back_substitute(Matrix A, Vector b)
{
  Vector x(b.n());
  for(size_t i=0; i<x.n(); ++i){ x(i) = 1; }

  double s;
  for(int i=x.n()-1; i>=0; --i)
  {
    s=0;
    for(int j=x.n()-1; j>i; --j)
    {
      s += A(i,j) * x(j);
    }
    x(i) = (b(i) - s)/A(i,i);
  }

  return x;
}

ostream & kry::operator<<(ostream &o, const Vector x)
{
  o << setprecision(6) << fixed;
  o << "[";
  for(size_t i=0; i<x.n()-1; ++i)
  {
    o << x(i) << " ";
  }
  o << x(x.n()-1);
  o << "]";
  return o;
}

ostream & kry::operator<<(ostream &o, const Matrix A)
{
  o << setprecision(6) << fixed;
  for(size_t i=0; i<A.m(); ++i)
  {
    for(size_t j=0; j<A.n(); ++j)
    {
      o << A(i,j) << ",";
    }
    o << std::endl;
  }
  return o;
}

Column::Column(Matrix M, size_t idx) 
  :  rbegin{0}, rend{M.m()}, _idx{idx}, _M{M} {}

size_t Column::idx() const { return _idx; }

Matrix Column::M() const { return _M; }

bool kry::operator==(const Column c, const Vector x)
{
  if(c.M().m() != x.n()) return false;

  for(size_t i=0; i<x.n(); ++i)
  {
    if(c.M()(i,c.idx()) != x(i)) return false;
  }

  return true;
}

Column Column::operator=(const Vector x)
{
  if(M().m() < x.n()) { throw NON_CONFORMAL; }
  
  for(size_t i=0; i<x.n(); ++i)
  {
    M()(i,idx()) = x(i);
  }

  return *this;
}

Column Column::operator-=(const Vector x)
{
  if(M().m() < x.n()) { throw NON_CONFORMAL; }
  
  for(size_t i=0; i<x.n(); ++i)
  {
    M()(i,idx()) -= x(i);
  }

  return *this;
}

Column Column::operator+=(const Vector x)
{
  if(M().m() < x.n()) { throw NON_CONFORMAL; }
  
  for(size_t i=0; i<x.n(); ++i)
  {
    M()(i,idx()) += x(i);
  }

  return *this;
}

Column Column::operator()(size_t rbegin, size_t rend)
{
  this->rbegin = rbegin;
  this->rend = rend;
  return *this;
}

double & Column::operator()(size_t idx)
{
  return M()(idx, _idx);
}

Vector::Vector(const Column c)
  : Vector(c.rend - c.rbegin)
{
  for(size_t i=c.rbegin; i<c.rend; ++i) { (*this)(i) = c.M()(i,c.idx()); }
}

Vector Vector::operator=(const Column c)
{
  if(c.M().m() > n()) { throw NON_CONFORMAL; }
  
  for(size_t i=0; i<n(); ++i)
  {
      (*this)(i) = c.M()(i,c.idx());
  }

  return *this; 
}

SparseMatrix::SparseMatrix(size_t m, size_t n, size_t z)
  : _m{m}, _n{n}, _z{z},
    _r(alloc<size_t>(m), _mm_free),
    _c(alloc<size_t>(m*z), _mm_free),
    _v(alloc<double>(m*z), _mm_free)
{
  size_t *_rp = _r.get();
  for(size_t i=0; i<m; ++i){ _rp[i] = 0; }
  size_t *_cp = _c.get();
  for(size_t i=0; i<m*z; ++i){ _cp[i] = 111; }
  double *_vp = _v.get();
  for(size_t i=0; i<m*z; ++i){ _vp[i] = 0; }
}

SparseMatrix::SparseMatrix(size_t m, size_t n, size_t z,
        std::vector<size_t> rs, std::vector<size_t> cs, 
        std::vector<double> vs)
  : SparseMatrix(m,n,z)
{
  if(rs.size() > m   ||
     cs.size() > m*n ||
     vs.size() > m*n ||
     vs.size() != cs.size()) 
  { throw UNREAL; }
 
  size_t *_rp = _r.get(), *_cp = _c.get();
  double *_vp = _v.get();
  size_t j=0;
  for(size_t i=0; i<rs.size(); ++i)
  {
    if(rs[i] > z) { throw UNREAL; }
    _rp[i] = rs[i];
    for(size_t k=0; k<rs[i]; ++k, ++j)
    {
      if(cs[j] > n) { throw UNREAL; }
      _cp[i*z + k] = cs[j];
      _vp[i*z + k] = vs[j];
    }
  }

  if(j != vs.size()) { throw UNREAL; }

}

SparseMatrix SparseMatrix::Identity(size_t m, size_t n, size_t z)
{
  SparseMatrix A(m, n, z);

  size_t *_rp = A._r.get(), *_cp = A._c.get();
  double *_vp = A._v.get();
  for(size_t i=0; i<max(m,n); ++i)
  {
    _rp[i] = 1;
    _cp[i*z] = i;
    _vp[i*z] = 1;
  }
  return A;
}
    
double & SparseMatrix::operator()(size_t i, size_t j)
{
  if(i >= _m || j >= _n) { throw UNREAL; }
  
  size_t *_rp = _r.get(), *_cp = _c.get();
  double *_vp = _v.get();
  
  for(size_t k=0; k<_rp[i]; ++k)
  {
    if(_cp[_z*i + k] == j){return _vp[_z*i + k];}
  }

  //if we could not find an entry, add one

  //if there is no room throw
  if(r(i) == z()) throw OVERSTUFFED;

  //add a unit entry into the requested row column pair
  size_t k = r(i);
  _cp[_z*i + k] = j;
  _vp[_z*i + k] = 0;
  ++_rp[i];
 
  //return a reference to the newly added entry
  return _vp[_z*i + k];
}

double SparseMatrix::get(size_t i, size_t j)
{
  if(i >= _m || j >= _n) { throw UNREAL; }
  
  size_t *_rp = _r.get(), *_cp = _c.get();
  double *_vp = _v.get();
  
  for(size_t k=0; k<_rp[i]; ++k)
  {
    if(_cp[_z*i + k] == j){return _vp[_z*i + k];}
  }

  return 0;
}

size_t SparseMatrix::m() const { return _m; }
size_t SparseMatrix::n() const { return _n; }
size_t SparseMatrix::z() const { return _z; }

size_t SparseMatrix::r(size_t i) const { return _r.get()[i]; }
size_t SparseMatrix::c(size_t i, size_t j) const { return _c.get()[z()*i+j]; }
size_t SparseMatrix::c(size_t k) const { return _c.get()[k]; }
double SparseMatrix::v(size_t i, size_t j) const { return _v.get()[z()*i+j]; }
double SparseMatrix::v(size_t k) const { return _v.get()[k]; }

Vector kry::operator*(const SparseMatrix A, const Vector x)
{
  if(A.n() != x.n()) { throw NON_CONFORMAL; }
  LOG_FUNC();
  Vector Ax = Vector::Zero(A.m());

  size_t
    rpt = ceil(A.m()/(double)RT::thread_count()),
    nt = min(RT::thread_count(), A.m());

  LOG_KVP(rpt);
  LOG_KVP(nt);

  CountdownLatch cl(nt);
  for(size_t tid=0; tid<nt; ++tid)
  {
    Task t = [tid,A,x,Ax,rpt,&cl]()
    {
      size_t ibegin = tid*rpt,
             iend = ibegin + rpt;
      iend = min(iend, A.m());

      for(size_t i=ibegin; i<iend; ++i)
      {
        for(size_t j=0; j<A.r(i); ++j)
        {
          Ax(i) += A.v(i,j) * x(A.c(i,j));
        }
      }

      --cl;
    };
    RT::get().workers[tid]->mtx->lock();
    RT::get().workers[tid]->task = t;
    RT::get().workers[tid]->mtx->unlock();
    RT::get().workers[tid]->cnd->notify_all();
  }
  cl.wait();
  
  return Ax;
}

Vector kry::operator*(const SparseMatrix A, const Column c)
{
  Vector x = c;
  return A * x;
}

ostream & kry::operator<<(ostream &o, SparseMatrix A)
{
  /*
  for(size_t i=0; i<A.m(); ++i)
  {
    for(size_t j=0; j<A.n(); ++j)
    {
      o << A.get(i,j) << " ";
    }
    o << "\n";
  }
  return o;
  */
  for(size_t i=0; i<A.m(); ++i) { o << A.r(i) << " "; }
  o << std::endl;

  for(size_t i=0; i<A.m()*A.z(); ++i) { o << A.c(i) << " "; }
  o << std::endl;

  for(size_t i=0; i<A.m()*A.z(); ++i) { o << A.v(i) << " "; }
  o << std::endl;

  return o;
}

ColumnRange::ColumnRange(Matrix M, size_t begin, size_t end)
  : _begin{begin}, _end{end}, _M{M}
{}

size_t ColumnRange::begin() const { return _begin; }
size_t ColumnRange::end() const { return _end; }
Matrix ColumnRange::M() const { return _M; }
void ColumnRange::transpose() { _M.transpose(); }

Vector kry::operator*(const ColumnRange C, const Vector x)
{
  if(C.M().transposed) { return mul_mv_t(C.M(), x, C.begin(), C.end()); }
  else { return mul_mv(C.M(), x, C.begin(), C.end()); }
}

Vector kry::operator*(const ColumnRange C, const Column c)
{
  Vector x = c;
  return C * x;
}

//Rotator ---------------------------------------------------------------------

Rotator::Rotator(Vector x, size_t i, size_t j)
  : _i{i}, _j{j}
{
  xi = x(i);
  xj = x(j);

  compute_c_s();
}

Rotator::Rotator(Matrix A, size_t i, size_t j)
  : _i{i}, _j{j}
{
  xi = A(i, i); //annihilator
  xj = A(j, i); //annihilatee

  compute_c_s();
}

void Rotator::compute_c_s()
{
  double B = fmax(fabs(xi), fabs(xj));
  if(B == 0) { _c = 1, _s = 0; }
  else
  {
    double _xi = xi/B, _xj = xj/B;
    double v = sqrt(_xi*_xi + _xj*_xj);
    _c = _xi/v, _s = _xj/v;
  }
}

size_t Rotator::i() const { return _i; }
size_t Rotator::j() const { return _j; }
double Rotator::c() const { return _c; }
double Rotator::s() const { return _s; }

Vector kry::operator* (const Rotator r, const Vector v)
{
  Vector x = !v;
  double xi = x(r.i()), xj = x(r.j());

  x(r.i()) = r.c()*xi + r.s()*xj;
  x(r.j()) = (-r.s())*xi + r.c()*xj;

  return x;
}

Vector Rotator::apply(Vector x)
{
  double xi = x(i()), xj = x(j());

  x(i()) = c()*xi + s()*xj;
  x(j()) = (-s())*xi + c()*xj;

  return x;
}
    
Matrix Rotator::apply_left(Matrix M)
{
  double xi, xj;
  for(size_t k=0; k<M.n(); ++k)
  {
    xi = M(i(), k),
    xj = M(j(), k);

    M(i(), k) = c()*xi + s()*xj;
    M(j(), k) = (-s())*xi + c()*xj;
  }

  return M;
}

Matrix Rotator::apply_right(Matrix)
{
  //TODO
}

Matrix kry::operator* (const Rotator r, const Matrix M)
{
  Matrix A = !M;
  double xi, xj;

  for(size_t k=0; k<A.n(); ++k)
  {
    xi = M(r.i(), k),
    xj = M(r.j(), k);

    A(r.i(), k) = r.c()*xi + r.s()*xj;
    A(r.j(), k) = (-r.s())*xi + r.c()*xj;
  }

  return A;
}

Matrix kry::operator* (const Matrix, const Rotator)
{
  //TODO
}

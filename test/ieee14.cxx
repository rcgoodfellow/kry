#include "gtest/gtest.h"
#include "kry/PSSE.hxx"

using std::cout;
using std::endl;
using std::array;
using std::vector;

using namespace kry;

TEST(ieee18, go)
{
  input::psse::Source muffin14("ieee14.psse");
  cout << muffin14.basicReport() << endl;

  array<SparseMatrix, 2> YY = muffin14.ymatrix();
  cout << "Y:" << endl;
  cout << YY[0] << endl;

  cout << "YA:" << endl;
  cout << YY[1] << endl;

  JacobiMap jmap = muffin14.jmap();

  //flat start
  Vector initial = Vector::Zero(14*2);
  for(size_t i=14; i<14*2; ++i) { initial(i) = 1; }
  initial(14+0) = 1.06;
  initial(14+1) = 1.045;
  initial(14+2) = 1.01;
  initial(14+5) = 1.07;
  initial(14+7) = 1.09;

  Vector ps = muffin14.psch();

  NKPF nkpf(YY[0], YY[1], jmap, 14, initial, ps);
  nkpf();
}


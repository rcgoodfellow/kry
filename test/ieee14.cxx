#include "gtest/gtest.h"
#include "kry/PSSE.hxx"

using std::cout;
using std::endl;
using std::array;
using std::vector;
using std::fstream;

using namespace kry;

TEST(ieee18, go)
{
  input::psse::Source ieee14_src("ieee14.psse");
  
  fstream fs("ieee14_basic_report.txt", fstream::out);
  fs << ieee14_src.basicReport();
  fs.close();

  array<SparseMatrix, 2> YY = ieee14_src.ymatrix();

  JacobiMap jmap = ieee14_src.jmap();

  fs.open("Y.smatrix", fstream::out);
  fs << YY[0];
  fs.close();
  fs.open("YA.smatrix", fstream::out);
  fs << YY[1];
  fs.close();

  //flat start
  Vector initial = Vector::Zero(14*2);
  for(size_t i=14; i<14*2; ++i) { initial(i) = 1; }
  initial(14+0) = 1.06;
  initial(14+1) = 1.045;
  initial(14+2) = 1.01;
  initial(14+5) = 1.07;
  initial(14+7) = 1.09;

  Vector ps = ieee14_src.psch();

  std::cout << "ps: " << ps << endl;

  size_t n = jmap.size();

  NKPF nkpf(YY[0], YY[1], jmap, n, initial, ps);
  nkpf();
}


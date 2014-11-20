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
  input::psse::Source ieee14_src("systems/ieee14.psse");
  
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

  size_t n = jmap.size();

  NKPF nkpf(YY[0], YY[1], jmap, n, initial, ps);
  nkpf();
}

TEST(NewEngland, go)
{
  input::psse::Source ne_src("systems/case39.psse");
  cout << ne_src.basicReport();

  fstream fs("ne_basic_report.txt", fstream::out);
  fs << ne_src.basicReport();
  fs.close();

  array<SparseMatrix, 2> YY = ne_src.ymatrix();
  JacobiMap jmap = ne_src.jmap();

  //flat start
  Vector initial = Vector::Zero(39*2);
  for(size_t i=39; i<39*2; ++i) { initial(i) = 1; }
  initial(39+28) = 1.05048;
  initial(39+29) = 1.04750;
  initial(39+30) = 0.98200;
  initial(39+31) = 0.98310;
  initial(39+32) = 0.99720;
  initial(39+33) = 1.01230;
  initial(39+34) = 1.04930;
  initial(39+35) = 1.06350;
  initial(39+36) = 1.02780;
  initial(39+37) = 1.02650;
  initial(39+38) = 1.03000;

  Vector ps = ne_src.psch();
  size_t n = jmap.size();
  NKPF nkpf(YY[0], YY[1], jmap, n, initial, ps);
  nkpf();
}


#include "gtest/gtest.h"
#include "kry/PSSE.hxx"

using std::cout;
using std::endl;
using std::array;

using namespace kry;

TEST(ieee18, go)
{
  input::psse::Source muffin14("ieee14.psse");
  cout << muffin14.basicReport() << endl;

  array<SparseMatrix, 2> YY = input::psse::ymatrix(muffin14);

  cout << "Y:" << endl;
  cout << YY[0] << endl;

  cout << "YA:" << endl;
  cout << YY[1] << endl;
}

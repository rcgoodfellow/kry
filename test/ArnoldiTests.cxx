#include "kry/Math.hxx"
#include "kry/Arnoldi.hxx"
#include "gtest/gtest.h"

using namespace kry;

TEST(Arnoldi, Small)
{
  SparseMatrix A(5,5,4,
      {2,        3,            3,           4,                 2      },
      {0,  1,    0,   3,  2,   3,  2,  1,   2,  1,   4,  3,    3,  4  },
      {0.4,0.24, 0.74,0.4,0.3, 0.5,0.9,0.7, 0.5,0.83,0.7,0.65, 0.7,0.7});

  Vector b{0.47, 0.32, 0.34, 0.41, 0.28};
  Vector x0{0,0,0,0,0};

  Arnoldi arnoldi(5, A, x0, b);
  arnoldi();
  
  std::cout << "Q" << std::endl;
  std::cout << arnoldi.Q << std::endl;
   
  std::cout << "H" << std::endl;
  std::cout << arnoldi.H << std::endl;

  std::cout << "d" << std::endl;
  std::cout << arnoldi.d << std::endl;
  
  std::cout << "t" << std::endl;
  std::cout << arnoldi.t << std::endl;

  std::cout << "xn" << std::endl;
  std::cout << arnoldi.xn << std::endl;

}

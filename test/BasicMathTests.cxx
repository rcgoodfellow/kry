#include "kry/Math.hxx"
#include "gtest/gtest.h"

using namespace kry;
using std::cout;
using std::endl;

TEST(Vector, Basics)
{
  Vector x = Vector::Zero(10);
  EXPECT_EQ(10UL, x.n());
  for(size_t i=0; i<10; ++i)
  {
    EXPECT_DOUBLE_EQ(0.0, x(i));
  }
  
  Vector y{1,2,3,4,5}, z{2,4,6};

  EXPECT_EQ(5UL, y.n());
  EXPECT_EQ(3UL, z.n());
  EXPECT_DOUBLE_EQ(3, y(2));
  EXPECT_DOUBLE_EQ(4, z(1));
}

TEST(Vector, Add)
{
  Vector x{1,2,3,4,5}, y{2,4,6,8,10};

  Vector z = x + y;

  EXPECT_DOUBLE_EQ(12, z(3));

  Vector xx = !x;
  xx += y;

  EXPECT_TRUE(xx == z);

  Vector xxx = x;
  xxx += y;

  EXPECT_TRUE(xxx == xx);
  EXPECT_TRUE(xx == x);
  EXPECT_TRUE(x == z);
}

TEST(Vector, Subtract)
{
  Vector x{2,4,6,8,10}, y{1,2,3,4,5};

  Vector z = x - y;

  EXPECT_TRUE(z == y);

  Vector xx = !x;
  xx -= y;

  EXPECT_TRUE(xx == z);

}

TEST(Vector, Dot)
{
  Vector x{1,3,5,7,9}, y{0,2,4,6,8};

  double z = x * y;

  EXPECT_EQ(1*0+2*3+4*5+6*7+8*9, z); 
}

TEST(Vector, DivScale)
{
  Vector x{2,4,6,8,10};
  double s{2};

  Vector y = x / s;

  EXPECT_DOUBLE_EQ(1, y(0));
  EXPECT_DOUBLE_EQ(2, y(1));
  EXPECT_DOUBLE_EQ(3, y(2));
  EXPECT_DOUBLE_EQ(4, y(3));
  EXPECT_DOUBLE_EQ(5, y(4));
  
  Vector z = !x;
  z /= s;

  EXPECT_TRUE(z == y);
}

TEST(Vector, MulScale)
{
  Vector x{2,4,6,8,10};
  double s{2};

  Vector y = x * s;

  Vector xx{4,8,12,16,20};

  EXPECT_TRUE(y == xx);

  Vector z = x;
  z *= s;
  EXPECT_TRUE(z == y);
}

TEST(Vector, Norm)
{
  Vector x{1,3,5,7,9};

  double nx = norm(x);
  double _nx_{12.845232578665129};

  EXPECT_TRUE(_nx_ == nx);
}

TEST(Matrix, MulVecBasic)
{
  Matrix A = Matrix::Identity(5,5);
  Vector x{1,2,3,4,5};

  Vector Ax = A * x;

  EXPECT_TRUE(Ax == x);

}

TEST(Matrix, MulVec)
{
  Matrix A = Matrix(5,5,
      {98,3,40,4,5,
       2,52,15,6,45,
       87,63,5,89,83,
       32,79,82,98,74,
       9,28,73,4,7});

  Vector x{3, 6, 90, 23, 64};

  Vector Ax = A * x;

  Vector _Ax_{4324,  4686,  8448, 14940,  7305};

  EXPECT_TRUE(_Ax_ == Ax);
}

TEST(Matrix, MulId)
{
  Matrix A = Matrix::Identity(5,5);
  
  Matrix B = Matrix(5,5,
      {98,3,40,4,5,
       2,52,15,6,45,
       87,63,5,89,83,
       32,79,82,98,74,
       9,28,73,4,7});

  Matrix AB = A * B;

  EXPECT_TRUE(AB == B);
}

TEST(Matrix, Mul)
{
  Matrix A = Matrix(5,5,
      {98,3,40,4,5,
       2,52,15,6,45,
       87,63,5,89,83,
       32,79,82,98,74,
       9,28,73,4,7});
  
  Matrix B = Matrix(5,5,
      {98,3,40,4,5,
       2,52,15,6,45,
       87,63,5,89,83,
       32,79,82,98,74,
       9,28,73,4,7});

  Matrix AB = A * B;

  Matrix _AB_(5,5,
    {13263,  3426,  4858,  4382,  4276,
     2202,  5389,  4712,  2423,  4354,
     12682, 13207, 17807, 10225, 10852,
     14230, 19184, 16313, 17800, 18291,
     7480,  6594,  1984,  7121,  7709});

  EXPECT_TRUE(AB == _AB_);

  EXPECT_DOUBLE_EQ(16313, _AB_(3, 2));

  Matrix Z = !A;
  Z *= B;

  EXPECT_TRUE(Z == AB);
}

TEST(Matrix, ColRef)
{
  Matrix A = Matrix::Identity(5,5);    
  Vector x1{0,0,1,0,0},
         x2{2,2,2,2,2};
  
  Column c2 = A.C(2);
  EXPECT_TRUE(c2 == x1);

  A.C(2) = Vector{2,2,2,2,2};

  EXPECT_TRUE(c2 == x2);

  Vector x3(5);
  x3 = A.C(2);

  EXPECT_TRUE(x3 == x2);
}

TEST(SparseMatrix, Basics)
{
  SparseMatrix A = SparseMatrix::Identity(5,5,3);

  EXPECT_EQ(A.m(), 5UL);
  EXPECT_EQ(A.n(), 5UL);
  EXPECT_EQ(A.z(), 3UL);

  EXPECT_DOUBLE_EQ(1, A(4, 4));
}

TEST(SparseMatrix, SparseMatrixVecMul2)
{
  SparseMatrix A(5, 5, 3,
      {1,   2,       2,       3,           2      },
      {0,   1,  4,   2,  4,   1,  3,  4,   0,  4  },
      {1.0, 1.0,7.0, 1.0,3.3, 2.2,1.0,1.1, 0.4,1.0});
  
  Vector x{1,2,3,4,5};

  Vector Ax = A * x;

  Vector _Ax_{1, 37, 19.5, 13.9, 5.4};

  EXPECT_TRUE(Ax == _Ax_);
}

TEST(Rotator, Vector)
{
  //Random vector of size 10 with mean 50 and variance 10
  Vector x = Vector::Random(10, 50, 10);

  Rotator r(x, 4, 7);

  Vector rx = r * x;

  EXPECT_NEAR(0, rx(7), 1e-14);  //TODO: lower this error bound
  
  r.apply(x);

  EXPECT_NEAR(0, x(7), 1e-14); 
}

TEST(Rotator, Matrix)
{
  Matrix A = Matrix::Random(10, 10, 50, 10);
  Rotator r(A, 4 ,7); //in column 4 eliminate row 7

  Matrix rA = r * A;
  
  EXPECT_NEAR(0, rA(7,4), 1e-14);

  cout << A << endl;
  r.apply_left(A);
  cout << A << endl;
  
  EXPECT_NEAR(0, A(7,4), 1e-14);
}



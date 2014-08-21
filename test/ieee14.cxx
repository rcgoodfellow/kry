#include "gtest/gtest.h"
#include "kry/PSSE.hxx"

TEST(ieee18, go)
{
  kry::input::psse::Source muffin14("ieee14.psse");
  std::cout << muffin14.basicReport() << std::endl;
}

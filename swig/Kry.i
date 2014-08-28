/******************************************************************************
 * libKrylov
 * =========
 * Kry.i is the swig interface file for libKrylov
 *
 * 28 August 2014
 * ~ ry
 *****************************************************************************/
 %module Kry
 %include "std_string.i"
 %{
 #include "kry/PSSE.hxx"
 #include "kry/Math.hxx"
 %}

namespace kry { namespace input { namespace psse {
struct Source
{
  Source(std::string fn); 
  std::string basicReport();
};
}}}

namespace kry {


}



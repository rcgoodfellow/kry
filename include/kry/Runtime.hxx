/******************************************************************************
 * libKrylov
 * =========
 * Utility.hxx contians the classes and functions that make up the libKrylov
 * runtime system
 *
 * 19 August 2014
 * ~ ry
 * ***************************************************************************/
#ifndef KRY_RUNTIME_HXX
#define KRY_RUNTIME_HXX

#include <vector>
#include <functional>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <memory>
#include <fstream>
#include <iostream>

#include "kry/Utility.hxx"

namespace kry
{

class RT;
class Worker;

using Task = std::function<void()>;
  
class RT
{
  public:
    RT();
    static RT& get();
    static size_t thread_count();

    std::vector<Worker*> workers;
    static constexpr size_t MEPT{5}; //minimum elements per thread
    std::fstream logger{"kry.log", std::fstream::out};
};

#define BLU "\e[0;34m"
#define RED "\e[0;31m"
#define YEL "\e[0;33m"
#define PRP "\e[0;35m"
#define GRN "\e[0;32m"
#define CYN "\e[0;36m"
#define GRY "\e[0;30m"
#define NRM "\e[0m"

template<class T>
void log_kvp(const char* key, T value)
{
  RT::get().logger 
    << key << "=" << std::to_string(value) 
    << std::endl;
}

#ifdef DEBUG
#define LOG_FUNC() RT::get().logger << __PRETTY_FUNCTION__ << std::endl;
#define LOG_KVP(__X__) log_kvp(#__X__, __X__)
#else
#define LOG_FUNC() ;;
#define LOG_KVP(__X__) ;;
#endif

template <class MatrixType>
void log_matrix(const char* name, const MatrixType &M)
{
  std::string fn{name};
  fn += ".matrix";
  std::fstream fs(fn, std::fstream::out);
  fs << M;
  fs.close();
}

template <class VectorType>
void log_vector(const char* name, const VectorType &v)
{
  std::string fn{name};
  fn += ".vector";
  std::fstream fs(fn, std::fstream::out);
  fs << v;
  fs.close();
}

template <class VectorType, class SuffixType>
void log_vector(const char* name, SuffixType sx, const VectorType &v)
{
  std::string fn{name};
  fn += std::to_string(sx) + ".vector";
  std::fstream fs(fn, std::fstream::out);
  fs << v;
  fs.close();
}

#define LOG_MATRIX(__M__) log_matrix(#__M__, __M__)
#define LOG_VECTOR(__V__) log_vector(#__V__, __V__)

#define LOG_VECTOR_SX(__V__, __SX__) log_vector(#__V__, __SX__, __V__)

class Worker
{
  public:
    Worker();
    Worker(Worker&&) = default;
    void start();
    void work();

    Task task;
    std::unique_ptr<std::thread> thd;
    std::unique_ptr<std::mutex> mtx{new std::mutex()};
    std::unique_ptr<std::condition_variable> cnd{new std::condition_variable()};
};

}

#endif

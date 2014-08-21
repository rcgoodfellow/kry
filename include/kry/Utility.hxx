/******************************************************************************
 * libKrylov
 * =========
 * Utility.hxx contians usefull stuff used throught libKrylov
 *
 * 19 August 2014
 * ~ ry
 * ***************************************************************************/
#ifndef KRY_UTILITY_HXX
#define KRY_UTILITY_HXX

#include <string>
#include <mm_malloc.h>
#include <mutex>
#include <condition_variable>
#include <iostream>

#define KRY_DEFAULT_ALIGN 64

#define CRIME_SCENE \
  std::string(__FILE__) + std::string(":") + \
  std::to_string(__LINE__) + \
  std::string(":") + std::string(__func__)

#define NON_CONFORMAL \
  std::runtime_error("non-conformal operation:" + CRIME_SCENE);

#define UNREAL \
  std::runtime_error("unreal sparse matrix coordinates:" + CRIME_SCENE);

#define BAD_COLRANGE \
  std::runtime_error("bad column range:" + CRIME_SCENE);

#define BAD_COLREF \
  std::runtime_error("bad column reference:" + CRIME_SCENE);

namespace kry
{

//Allocate a typed, aligned chunk of memory without havint to do the typecast
//and sizeof math bit each time
template<class T>
T* alloc(size_t sz)
{
  return (T*)_mm_malloc(sizeof(T)*sz, KRY_DEFAULT_ALIGN);
}

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
  std::cout 
    << GRN << key << GRY "="
    << NRM << std::to_string(value)
    << std::endl;
}

#ifdef DEBUG
#define LOG_FUNC() printf("%s%s%s\n", BLU, __PRETTY_FUNCTION__, NRM);
#define LOG(__X__) printf __X__
#define LOG_KVP(__X__) log_kvp(#__X__, __X__)
#else
#define LOG_FUNC ;;
#define LOG(__X__) ;;
#define LOG_KVP(__X__) ;;
#endif

class CountdownLatch
{
  public:
    CountdownLatch(int);
    void wait();
    void set(int);
    void operator--();
    void operator--(int);
    void operator++();
    void operator++(int);
    int operator()() const;

  private:
    std::atomic<int> _cnt;
    std::shared_ptr<std::mutex> _mtx{new std::mutex};
    std::shared_ptr<std::condition_variable> _cnd{new std::condition_variable};
};

}

#endif

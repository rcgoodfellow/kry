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
};

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

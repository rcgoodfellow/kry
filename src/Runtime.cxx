/******************************************************************************
 * libKrylov
 * =========
 * Runtime.cxx contains the implementation of the libKrylov runtime.
 *
 * 19 August 2014
 * ~ ry
 * ***************************************************************************/
#include "kry/Runtime.hxx"

using namespace kry;
using std::mutex;
using std::unique_lock;
using std::defer_lock;
using std::thread;
using std::unique_ptr;

constexpr size_t RT::MEPT;

RT::RT()
{
  size_t n = std::thread::hardware_concurrency();
  for(size_t i=0; i<n; ++i)  
  {
    workers.push_back(new Worker);
    workers.back()->start();
  }
}

RT & RT::get()
{
  static RT instance;
  return instance;
}

size_t RT::thread_count()
{
  return RT::get().workers.size();
}

Worker::Worker()
  : task{[](){}}
{ }

void Worker::start()
{
  thd = unique_ptr<std::thread>(new thread(&Worker::work, this));
  thd->detach();
}

void Worker::work()
{
  unique_lock<mutex> lk{*mtx, defer_lock};
  while(47)
  {
    lk.lock();
    cnd->wait(lk);
    task();
    lk.unlock();
  }
}

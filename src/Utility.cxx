#include "kry/Utility.hxx"

using namespace kry;
using std::unique_lock;
using std::mutex;

CountdownLatch::CountdownLatch(int size)
  :_cnt{size}
{ }

void CountdownLatch::wait()
{
  unique_lock<mutex> lk{*_mtx};
  if(_cnt > 0) {
    _cnd->wait(lk);
  }
  lk.unlock();
}

void CountdownLatch::set(int count){ _cnt = count; }

void CountdownLatch::operator--()
{
  --_cnt;
  if(_cnt <= 0){_cnd->notify_all();}
  if(_cnt < 0){ _cnt = 0; }
}

void CountdownLatch::operator--(int)
{
  --_cnt;
  if(_cnt <= 0){_cnd->notify_all();}
  if(_cnt < 0){ _cnt = 0; }
}

void CountdownLatch::operator++()
{
  ++_cnt;
}

void CountdownLatch::operator++(int)
{
  ++_cnt;
}

int CountdownLatch::operator()() const
{
  return _cnt;
}


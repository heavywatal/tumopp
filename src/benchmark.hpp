/*! @file benchmark.hpp
    @brief Utilities for benchmarking
*/
#pragma once
#ifndef TUMOPP_BENCHMARK_HPP_
#define TUMOPP_BENCHMARK_HPP_

#ifndef _WIN32
  #include <wtl/resource.hpp>
#endif // _WIN32

#include <sstream>

namespace tumopp {

class Benchmark {
  public:
    void append(std::size_t size) {
        sst_ << size;
        append_ru();
        sst_ << "\n";
    }
    std::stringstream stringstream() const {
        std::stringstream sst;
        sst << "size\tmemory\tutime\tstime\n" << sst_.rdbuf();
        return sst;
    }
  private:
#ifdef _WIN32
    void append_ru() {}
#else
    void append_ru() {
        auto ru = wtl::getrusage<std::milli, std::mega>(epoch_);
        sst_ << "\t" << ru.maxrss
             << "\t" << ru.utime
             << "\t" << ru.stime;
    }
    const rusage epoch_ = wtl::getrusage();
#endif // _WIN32
    std::stringstream sst_;
};

} // namespace tumopp

#endif // TUMOPP_BENCHMARK_HPP_

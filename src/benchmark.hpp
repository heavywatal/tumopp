/*! @file benchmark.hpp
    @brief Utilities for benchmarking
*/
#pragma once
#ifndef TUMOPP_BENCHMARK_HPP_
#define TUMOPP_BENCHMARK_HPP_

#ifdef BENCHMARK
  #include <wtl/resource.hpp>
  #define MAKE_EPOCH_IF_BENCHMARK() wtl::getrusage()
  #define APPEND_IF_BENCHMARK(ost, size, epoch) WRITE_BENCHMARK(ost, size, epoch)
#else
  #define MAKE_EPOCH_IF_BENCHMARK() 0
  #define APPEND_IF_BENCHMARK(ost, size, epoch) epoch += 0
#endif // BENCHMARK

namespace tumopp {

#ifdef BENCHMARK
  //! Save current resource usage to #benchmark_
  template <class RU> inline
  void WRITE_BENCHMARK(std::ostream& ost, std::size_t size, const RU& epoch) {
      auto ru = wtl::getrusage<std::milli, std::mega>(epoch);
      ost << size << "\t"
          << ru.maxrss << "\t"
          << ru.utime << "\t"
          << ru.stime << "\n";
  }
#endif // BENCHMARK

constexpr const char* BENCHMARK_HEADER = "size\tmemory\tutime\tstime\n";

} // namespace tumopp

#endif // TUMOPP_BENCHMARK_HPP_

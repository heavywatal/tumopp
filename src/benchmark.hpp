/*! @file benchmark.hpp
    @brief Defines Benchmark class
*/
#pragma once
#ifndef TUMOPP_BENCHMARK_HPP_
#define TUMOPP_BENCHMARK_HPP_

#ifndef _WIN32
  #include <wtl/resource.hpp>
#endif // _WIN32

#include <sstream>

namespace tumopp {

//! Utility class for benchmarking
class Benchmark {
  public:
    Benchmark() {
        sst_ << "size"
#ifndef _WIN32
             << "\t" << wtl::rusage_header()
#endif // _WIN32
             << "\n";
    }
    //! Append current state to #sst_
    void append(std::size_t size) {
        sst_ << size
#ifndef _WIN32
             << "\t" << wtl::getrusage<std::milli, std::kilo>()
#endif // _WIN32
             << "\n";
    }
    //! Get TSV stream buffer
    std::streambuf* rdbuf() const {return sst_.rdbuf();}
  private:
    //! String in TSV format
    std::stringstream sst_{};
};

} // namespace tumopp

#endif // TUMOPP_BENCHMARK_HPP_

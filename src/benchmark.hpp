/*! @file benchmark.hpp
    @brief Defines Benchmark class
*/
#pragma once
#ifndef TUMOPP_BENCHMARK_HPP_
#define TUMOPP_BENCHMARK_HPP_

#include <wtl/resource.hpp>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <sstream>

namespace tumopp {

//! Utility class for benchmarking
class Benchmark {
  public:
    Benchmark() {
        sst_ << "size\t"
             << fmt::format("{}\n", fmt::join(wtl::rusage_names, "\t"));
    }
    //! Append current state to #sst_
    void append(std::size_t size) {
        sst_ << size
             << "\t"
             << fmt::format("{}\n", fmt::join(wtl::get_rusage<std::milli, std::kilo>(), "\t"));
    }
    //! Get TSV stream buffer
    std::streambuf* rdbuf() const {return sst_.rdbuf();}
  private:
    //! String in TSV format
    std::stringstream sst_{};
};

} // namespace tumopp

#endif // TUMOPP_BENCHMARK_HPP_

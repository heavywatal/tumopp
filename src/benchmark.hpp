/*! @file benchmark.hpp
    @brief Defines Benchmark class
*/
#pragma once
#ifndef TUMOPP_BENCHMARK_HPP_
#define TUMOPP_BENCHMARK_HPP_

#include <wtl/resource.hpp>
#include <fmt/base.h>
#include <fmt/ranges.h>

#include <iterator>

namespace tumopp {

//! Utility class for benchmarking
class Benchmark {
  public:
    Benchmark() {
        fmt::format_to(std::back_inserter(buffer_), "size\t{}\n",
                       fmt::join(wtl::rusage_names, "\t"));
    }
    //! Append current state to #buffer_
    void append(std::size_t size) {
        fmt::format_to(std::back_inserter(buffer_), "{}\t{}\n",
                       size, fmt::join(wtl::get_rusage<std::milli, std::kilo>(), "\t"));
    }
    //! Get string in TSV format
    const std::string& str() const {return buffer_;}
  private:
    //! String in TSV format
    std::string buffer_;
};

} // namespace tumopp

#endif // TUMOPP_BENCHMARK_HPP_

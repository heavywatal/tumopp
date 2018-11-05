/*! @file random.hpp
    @brief Utilities for random number generation
*/
#pragma once
#ifndef TUMOPP_RANDOM_HPP_
#define TUMOPP_RANDOM_HPP_

#ifdef SFMT_FOUND
  #include <sfmt.hpp>
#endif
#include <random>

namespace tumopp {

#ifdef SFMT_FOUND
    using urbg_t = wtl::sfmt19937_64;
#else
    using urbg_t = std::mt19937_64;
#endif

} // namespace tumopp

#endif // TUMOPP_RANDOM_HPP_

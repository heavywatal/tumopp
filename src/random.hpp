/*! @file random.hpp
    @brief Defines type alias `urbg_t`
*/
#pragma once
#ifndef TUMOPP_RANDOM_HPP_
#define TUMOPP_RANDOM_HPP_

#ifdef SFMT_FOUND
  #include <sfmt.hpp>
#endif
#include <random>

namespace tumopp {

//! Type alias of random number generator
#ifdef SFMT_FOUND
    using urbg_t = wtl::sfmt19937_64;
#else
    using urbg_t = std::mt19937_64;
#endif

} // namespace tumopp

#endif // TUMOPP_RANDOM_HPP_

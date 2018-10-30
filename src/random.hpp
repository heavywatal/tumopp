/*! @file random.hpp
    @brief Utilities for random number generation
*/
#pragma once
#ifndef TUMOPP_RANDOM_HPP_
#define TUMOPP_RANDOM_HPP_

#ifdef SFMT_FOUND
  #include <sfmt.hpp>
#endif
#include <wtl/random.hpp>

namespace tumopp {

//! Getter of static instance
inline auto& engine64() {
#ifdef SFMT_FOUND
    return wtl::sfmt64();
#else
    return wtl::mt64();
#endif
}

}

#endif // TUMOPP_RANDOM_HPP_

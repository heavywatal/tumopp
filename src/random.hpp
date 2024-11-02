/*! @file random.hpp
    @brief Defines type alias `urbg_t`
*/
#pragma once
#ifndef TUMOPP_RANDOM_HPP_
#define TUMOPP_RANDOM_HPP_

#include <pcglite/pcglite.hpp>
#include <random>

namespace tumopp {

//! Type alias of random number generator
using urbg_t = pcglite::pcg64;

} // namespace tumopp

#endif // TUMOPP_RANDOM_HPP_

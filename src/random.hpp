#pragma once
#ifndef TUMOPP_RANDOM_HPP
#define TUMOPP_RANDOM_HPP

#ifdef SFMT_FOUND
  #include <sfmt.hpp>
#endif
#include <wtl/random.hpp>

namespace tumopp {

inline auto& engine64() {
#ifdef SFMT_FOUND
    return wtl::sfmt64();
#else
    return wtl::mt64();
#endif
}

}

#endif//TUMOPP_RANDOM_HPP

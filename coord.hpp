// -*- mode: c++; coding: utf-8 -*-
/*! @file coord.hpp
    @brief Coordinate system
*/
#pragma once
#ifndef COORD_HPP_
#define COORD_HPP_

#include <vector>
#include <algorithm>
#include <numeric>


/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

template <class T>
std::vector<T>& operator+=(std::vector<T>& lhs, const std::vector<T>& rhs) {
    assert(lhs.size() == rhs.size());
    std::transform(lhs.begin(), lhs.end(), rhs.begin(),
                   lhs.begin(), std::plus<T>());
    return lhs;
}

template <class T>
std::vector<T>& operator-=(std::vector<T>& lhs, const std::vector<T>& rhs) {
    assert(lhs.size() == rhs.size());
    std::transform(lhs.begin(), lhs.end(), rhs.begin(),
                   lhs.begin(), std::minus<T>());
    return lhs;
}

template <class T>
std::vector<T> operator+(std::vector<T> lhs, const std::vector<T>& rhs) {
    return lhs += rhs;
}

template <class T>
std::vector<T> operator-(std::vector<T> lhs, const std::vector<T>& rhs) {
    return lhs -= rhs;
}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

class Coord {
  public:
    Coord() = default;
    virtual ~Coord() = default;
    std::vector<std::vector<int>> neighbors(const std::vector<int>& v) const {
        std::vector<std::vector<int>> output = directions(v.size());
        for (auto& d: output) {
            d += v;
        }
        return output;
    }
    std::vector<int> outward(const std::vector<int>& v) const {
        const auto candidates = neighbors(v);
        return *std::max_element(candidates.begin(), candidates.end(),
                                 [this](const std::vector<int>& lhs, const std::vector<int>& rhs) {
            return distance(lhs) < distance(rhs);
        });
    }
    virtual std::vector<std::vector<int>> directions(const size_t dimensions) const = 0;
    virtual size_t distance(const std::vector<int>& v) const = 0;
};

class Orthogonal: public Coord {
  public:
    Orthogonal() = default;
    size_t distance(const std::vector<int>& v) const {
        return std::accumulate(v.begin(), v.end(), 0, [](const int lhs, const int rhs) {
            return std::abs(lhs) + std::abs(rhs);
        });
    }
    std::vector<std::vector<int>> directions(const size_t dimensions) const {
        std::vector<std::vector<int>> output;
        output.reserve(2 * dimensions);
        std::vector<int> v(dimensions, 0);
        v.back() += 1;
        do {
            output.push_back(v);
        } while (std::next_permutation(v.begin(), v.end()));
        v.assign(dimensions, 0);
        v.front() -= 1;
        do {
            output.push_back(v);
        } while (std::next_permutation(v.begin(), v.end()));
        return output;
    }
};


class Lattice: public Coord {
  public:
    Lattice() = default;
    size_t distance(const std::vector<int>& v) const {
        return std::abs(*std::max_element(v.begin(), v.end(), [](const int lhs, const int rhs) {
            return std::abs(lhs) < std::abs(rhs);
        }));
    }
    std::vector<std::vector<int>> directions(const size_t dimensions) const {
        std::vector<std::vector<int>> output;
        output.reserve(std::pow(3, dimensions) - 1);
        for (const int x: {-1, 0, 1}) {
            for (const int y: {-1, 0, 1}) {
                if (dimensions == 2) {
                    if (x == 0 && y == 0) continue;
                    output.push_back({x, y});
                    continue;
                }
                for (const int z: {-1, 0, 1}) {
                    if (x == 0 && y == 0 && z == 0) continue;
                    output.push_back({x, y, z});
                }
            }
        }
        return output;
    }
};

class Hex: public Coord {
  public:
    Hex() = default;
    size_t distance(const std::vector<int>& v) const {
        std::vector<size_t> absv;
        absv.reserve(v.size() * 2);
        for (auto x: v) {
            absv.push_back(std::abs(x));
        }
        absv.push_back(std::abs(v[0] + v[1]));
        if (v.size() == 3) {
            absv.push_back(std::abs(v[0] + v[2]));
        }
        return *std::max_element(absv.begin(), absv.end());
    }

    std::vector<std::vector<int>> directions(const size_t dimensions) const {
        std::vector<std::vector<int>> output;
        std::vector<int> v{-1, 0, 1};
        if (dimensions == 2) {
            output.reserve(6);
            do {
                output.push_back({v[0], v[1]});
            } while (std::next_permutation(v.begin(), v.end()));
        }
        else {
            output.reserve(12);
            do {
                output.push_back({v[0], v[1], 0});
            } while (std::next_permutation(v.begin(), v.end()));
            output.push_back({0, 0, -1});
            output.push_back({1, 0, -1});
            output.push_back({1, -1, -1});
            output.push_back({0, 0, 1});
            output.push_back({-1, 0, 1});
            output.push_back({-1, 1, 1});
        }
        return output;
    }
};



#endif /* COORD_HPP_ */

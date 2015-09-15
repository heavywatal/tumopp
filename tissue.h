// -*- mode: c++; coding: utf-8 -*-
/*! @file tissue.h
    @brief Interface of Tissue class
*/
#pragma once
#ifndef TISSUE_H_
#define TISSUE_H_
#include <cmath>
#include <cassert>
#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <string>
#include <numeric>

#include <boost/functional/hash.hpp>

#include "gland.h"

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


namespace std {
  template <> struct hash<std::vector<int>> {
    size_t operator() (const std::vector<int>& v) const {
        return boost::hash_range(v.begin(), v.end());
    }
  };
}

namespace boost {
    namespace program_options {
        class options_description;
    }
}


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

class Tissue {
  public:
    static size_t DIMENSIONS() {return DIMENSIONS_;}

    //! Constructor
    Tissue(const std::vector<int>& origin=std::vector<int>(DIMENSIONS_));

    //! Set coordinate function object
    template <class FuncObj>
    void set_coord() {coord_func_.reset(new FuncObj());}

    //! Initiate tumor with neighboring glands regularly
    void init_regularly();
    //! Initiate tumor with neighboring glands randomly
    void init_randomly();
    //! Mark first cells with mutation
    void stain();

    //! Big bang model
    void grow_random(const size_t max_size);
    //! Another extreme
    void grow_even(const size_t max_size);
    //! in-between
    void grow_poisson(const size_t max_size);

    //! Return tumor state as TSV string
    std::string snapshot() const;
    //! Return mutation_coords_ as TSV string
    std::string mutation_history() const;

    std::string snapshot_header() const;
    std::string evolution_history() const;

    friend std::ostream& operator<< (std::ostream&, const Tissue&);

    //! Unit test
    static void unit_test();
    static boost::program_options::options_description& opt_description();

  private:
    //! Dimensions: {1, 2, 3}
    static size_t DIMENSIONS_;

    //! Emplace daughter gland at the specified coord
    void emplace(const std::vector<int>& current_coord, Gland&& daughter);
    //! Emplace daughter gland and push other glands outward
    void push(Gland&& daughter, std::vector<int>* current_coord,
              const std::vector<int>& direction={});
    //! Push if not exists an empty neighbor
    void push_fill(Gland&& daughter, const std::vector<int>& current_coord,
                   const std::vector<int>& direction={});
    //! Change direction before push_fill
    void walk_fill(Gland&& daughter, const std::vector<int>& current_coord);

    std::unique_ptr<Coord> coord_func_;

    //! key: coords, value: gland
    std::unordered_map<std::vector<int>, Gland> tumor_;
    //! The coordinates of the existing glands
    std::vector<std::vector<int>> coords_;
    //! The coordinates of the past mutations
    std::vector<std::vector<int>> mutation_coords_;
    //! Timing of mutations (tumor size as proxy)
    std::vector<size_t> mutation_stages_;
    std::vector<std::string> evolution_history_;

    const std::string sep_ = "\t";
};

#endif /* TISSUE_H_ */

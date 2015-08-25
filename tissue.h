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

#include <boost/functional/hash.hpp>

#include "gland.h"

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace boost {
    namespace program_options {
        class options_description;
    }
}

class Hex {
  public:
    Hex(const int x, const int y): x_(x), y_(y) {}
    Hex(const int x, const int y, int z): Hex{x, y} {assert(z += y == -x);}
    Hex(const std::vector<int>& v): x_(v[0]), y_(v[1]) {}
    int x() const {return x_;}
    int y() const {return y_;}
    int z() const {return -x_ -y_;}
    std::vector<int> vec() const {return {x_, y_};}

    int radius() const {
        int d = std::abs(x_);
        d += std::abs(y_);
        d += std::abs(z());
        return d /= 2;
    }

    std::vector<Hex> neighbors() const {
        std::vector<Hex> output;
        output.reserve(6);
        for (auto& v: directions()) {
            output.push_back(v += *this);
        }
        return output;
    }

    static int distance(const Hex& lhs, const Hex& rhs) {
        return (lhs - rhs).radius();
    }
    static std::vector<Hex> directions() {
        std::vector<Hex> output;
        output.reserve(6);
        std::vector<int> v{-1, 0, 1};
        do {
            output.push_back(Hex(v[0], v[1]));
        } while (std::next_permutation(v.begin(), v.end()));
        return output;
    }
    Hex& operator+= (const Hex& rhs) {
        x_ += rhs.x_;
        y_ += rhs.y_;
        return *this;
    }
    Hex operator+ (const Hex& rhs) const {
        return Hex(x_ + rhs.x_, y_ + rhs.y_);
    }
    Hex operator- (const Hex& rhs) const {
        return Hex(x_ - rhs.x_, y_ - rhs.y_);
    }
    bool operator= (const Hex& rhs) const {
        return x_ == rhs.x_ && y_ == rhs.y_;
    }
    friend std::ostream& operator<< (std::ostream& ost, const Hex& hex) {
        return ost << "[" << hex.x_ << ", " << hex.y_ << ", " << hex.z() << "]";
    }

  private:
    int x_;
    int y_;
};

namespace std {
  template <> struct hash<Hex> {
    size_t operator() (const Hex& hex) const {
        size_t h = 0;
        boost::hash_combine(h, hex.x());
        boost::hash_combine(h, hex.y());
        return h;
    }
  };
  template <> struct hash<std::vector<int>> {
    size_t operator() (const std::vector<int>& v) const {
        return boost::hash_range(v.begin(), v.end());
    }
  };
}

class Tissue {
    typedef std::vector<int> coord_t;
  public:
    //! Constructor
    Tissue(const std::vector<int>& origin): coords_{1, origin} {
        tumor_.emplace(origin, Gland());
        init_regularly();
    }
    //! Constructor
    Tissue(const size_t dimensions=2): Tissue{std::vector<int>(dimensions)} {};

    //! Initiate tumor with 2^D glands regularly
    void init_regularly();
    //! Initiate tumor with 2^D glands randomly
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
    //! Emplace daughter gland at the specified coord
    void emplace(const std::vector<int>& current_coord, Gland&& daughter);
    //! Emplace daughter gland and push other glands outward
    void push(Gland&& daughter, std::vector<int>* current_coord, const std::vector<int>& direction);
    //! Emplace daughter gland at a neighboring site and push if needed
    void fill_neighbor(Gland&& daughter, const std::vector<int>& current_coord);
    void push_neighbor(Gland&& daughter, const std::vector<int>& current_coord,
                       const std::vector<int>& direction={});
    void hex_neighbor(Gland&& daughter, const std::vector<int>& current_coord,
                       const std::vector<int>& direction={});

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

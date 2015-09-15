// -*- mode: c++; coding: utf-8 -*-
/*! @file tissue.cpp
    @brief Implementation of Tissue class
*/
#include "tissue.h"

#include <cmath>
#include <iostream>
#include <sstream>

#include <boost/program_options.hpp>

#include "cxxwtils/prandom.hpp"
#include "cxxwtils/iostr.hpp"
#include "cxxwtils/gz.hpp"
#include "cxxwtils/debug.hpp"
#include "cxxwtils/algorithm.hpp"


size_t Tissue::DIMENSIONS_ = 2;

//! Program options
/*! @return Program options description

    Command line option | Symbol | Variable
    --------------------| ------ | -------------------------
    `-D,--dimensions`   | -      | Tissue::DIMENSIONS_
*/
boost::program_options::options_description& Tissue::opt_description() {
    namespace po = boost::program_options;
    static po::options_description desc{"Tissue"};
    desc.add_options()
        ("dimensions,D", po::value<size_t>(&DIMENSIONS_)->default_value(DIMENSIONS_))
    ;
    return desc;
}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

class Orthogonal: public std::vector<int> {
  public:
    Orthogonal() = default;
    Orthogonal(const std::vector<int>& v): std::vector<int>(v) {}
    size_t distance() const {
        return std::accumulate(begin(), end(), 0, [](const int lhs, const int rhs) {
            return std::abs(lhs) + std::abs(rhs);
        });
    }
    std::vector<Orthogonal> neighbors() const {
        const auto dirs = directions(size());
        std::vector<Orthogonal> output;
        output.reserve(dirs.size());
        for (const auto& v: dirs) {
            output.push_back(v + *this);
        }
        return output;
    }

    static std::vector<std::vector<int>> directions(const size_t dimensions) {
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

class Lattice: public std::vector<int> {
  public:
    Lattice() = default;
    Lattice(const std::vector<int>& v): std::vector<int>(v) {}
    size_t distance() const {
        return *std::max_element(begin(), end(), [](const int lhs, const int rhs) {
            return std::abs(lhs) < std::abs(rhs);
        });
    }
    std::vector<Lattice> neighbors() const {
        const auto dirs = directions(size());
        std::vector<Lattice> output;
        output.reserve(dirs.size());
        for (const auto& v: dirs) {
            output.push_back(v + *this);
        }
        return output;
    }
    static std::vector<std::vector<int>> directions(const size_t dimensions) {
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

class Hex: public std::vector<int> {
  public:
    Hex() = default;
    Hex(const std::vector<int>& v): std::vector<int>(v) {}
    size_t distance() const {
        std::vector<int> absv;
        absv.reserve(this->size() * 2);
        for (auto x: *this) {
            absv.push_back(std::abs(x));
        }
        absv.push_back(std::abs(xy()));
        if (size() == 3) {
            absv.push_back(std::abs(xz()));
        }
        return *std::max_element(absv.begin(), absv.end());
    }

    std::vector<Hex> neighbors() const {
        const auto dirs = directions(size());
        std::vector<Hex> output;
        output.reserve(dirs.size());
        for (const auto& v: dirs) {
            output.push_back(v + *this);
        }
        return output;
    }

    static std::vector<std::vector<int>> directions(const size_t dimensions) {
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

  private:
    int xy() const {return -this->operator[](0) -this->operator[](1);}
    int xz() const {return -this->operator[](0) -this->operator[](2);}
};

//namespace std {
//  template <> struct hash<Hex> {
//    size_t operator() (const Hex& hex) const {
//        size_t h = 0;
//        boost::hash_combine(h, hex.x());
//        boost::hash_combine(h, hex.y());
//        return h;
//    }
//  };
//}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

Tissue::coord_t Tissue::outward(const coord_t& current) {
    const auto candidates = current.neighbors();
    return *std::max_element(candidates.begin(), candidates.end(),
                             [](const coord_t& lhs, const coord_t& rhs) {
        return lhs.distance() < rhs.distance();
    });
}

void Tissue::grow_random(const size_t max_size) {HERE;
    coords_.reserve(max_size);
    while (tumor_.size() < max_size) {
        auto current_coord = *wtl::prandom().choice(coords_.begin(), coords_.end());
        auto& parent = tumor_[current_coord];
        Gland daughter = parent;
        if (parent.bernoulli_apoptosis()) {
            parent = std::move(daughter);
        } else {
//            push(std::move(daughter), &current_coord);
//            push_fill(std::move(daughter), current_coord);
            walk_fill(std::move(daughter), current_coord);
        }
        if (Gland::bernoulli_mutation()) {
            tumor_[current_coord].mutate();
            mutation_coords_.push_back(std::move(current_coord));
            mutation_stages_.push_back(tumor_.size());
        }
    }
}

void Tissue::grow_even(const size_t max_size) {HERE;
    coords_.reserve(max_size);
    size_t age = 1;
    evolution_history_.reserve(max_size);
    evolution_history_.push_back(snapshot());
    for (auto it=coords_.rbegin(); tumor_.size() < max_size; ++it) {
        while (it != coords_.rend() && tumor_[*it].age() == age) {++it;}
        if (it == coords_.rend()) {
            it = coords_.rbegin();
            ++age;
        }
        auto& mother = tumor_[*it];
        mother.stamp(age);
        auto current_coord = *it;
        Gland daughter = mother;
        if (mother.bernoulli_apoptosis()) {
            mother = std::move(daughter);
        } else {
            push(std::move(daughter), &current_coord);
//            push_fill(std::move(daughter), current_coord);
//            walk_fill(std::move(daughter), current_coord);
        }
        if (Gland::bernoulli_mutation()) {
            tumor_[current_coord].mutate();
            mutation_coords_.push_back(std::move(current_coord));
            mutation_stages_.push_back(tumor_.size());
        }
        if (tumor_.size() <= 256) {
            evolution_history_.push_back(snapshot());
        }
    }
}

void Tissue::grow_poisson(const size_t max_size) {HERE;
    //! @todo
}

std::string Tissue::evolution_history() const {
    return wtl::str_join(evolution_history_, "");
}

void Tissue::emplace(const std::vector<int>& coord, Gland&& daughter) {
    coords_.push_back(coord);
    tumor_.emplace(coord, std::move(daughter));
}

void Tissue::push(Gland&& daughter, std::vector<int>* coord, const std::vector<int>& direction) {
    // TODO: rate-limiting
    auto new_dir = direction;
    if (direction.empty()) {
        const auto directions = coord_t::directions(coord->size());
        new_dir = *wtl::prandom().choice(directions.begin(), directions.end());
    }
    *coord += new_dir;
    auto it = tumor_.find(*coord);
    if (it == tumor_.end()) {
        emplace(*coord, std::move(daughter));
    } else {
        push(std::move(it->second), coord, new_dir);
        it->second = std::move(daughter);
        *coord = it->first;
    }
}

//! @todo
void Tissue::push_fill(Gland&& daughter, const std::vector<int>& current_coord,
                       const std::vector<int>& direction) {
    static const auto directions = coord_t::directions(current_coord.size());
    const auto neighbors = coord_t(current_coord).neighbors();
    std::vector<coord_t> empty_neighbors;
    empty_neighbors.reserve(neighbors.size());
    for (auto& x: neighbors) {
        if (tumor_.find(x) == tumor_.end()) {
            empty_neighbors.push_back(x);
        }
    }
    if (empty_neighbors.empty()) {
        auto new_dir = direction;
        if (direction.empty()) {
            new_dir = *wtl::prandom().choice(directions.begin(), directions.end());
        }
        auto new_coord = current_coord + new_dir;
        push_fill(std::move(tumor_[current_coord]), new_coord, new_dir);
    } else {
        auto new_coord = *std::min_element(empty_neighbors.begin(), empty_neighbors.end(),
            [](const coord_t& x, const coord_t& y) {
                return x.distance() < y.distance();
        });
        emplace(new_coord, std::move(tumor_[current_coord]));
    }
    tumor_[current_coord] = std::move(daughter);
}

//! @todo
void Tissue::walk_fill(Gland&& daughter, const std::vector<int>& current_coord) {
    static const auto directions = coord_t::directions(current_coord.size());
    const auto neighbors = coord_t(current_coord).neighbors();
    std::vector<coord_t> empty_neighbors;
    empty_neighbors.reserve(neighbors.size());
    for (auto& x: neighbors) {
        if (tumor_.find(x) == tumor_.end()) {
            empty_neighbors.push_back(x);
        }
    }
    if (empty_neighbors.empty()) {
        auto new_coord = *wtl::prandom().choice(neighbors.begin(), neighbors.end());
        walk_fill(std::move(tumor_[current_coord]), new_coord);
    } else {
        auto new_coord = *std::min_element(empty_neighbors.begin(), empty_neighbors.end(),
            [](const coord_t& x, const coord_t& y) {
                return x.distance() < y.distance();
        });
        emplace(new_coord, std::move(tumor_[current_coord]));
    }
    tumor_[current_coord] = std::move(daughter);
}

void Tissue::init_regularly() {HERE;
    for (const auto& coord: coord_t(coords_.front()).neighbors()) {
        emplace(coord, Gland());
    }
}

void Tissue::stain() {HERE;
    for (auto& item: tumor_) {
        item.second.mutate();
        mutation_coords_.push_back(item.first);
        mutation_stages_.push_back(mutation_coords_.size());
    }
}

std::string Tissue::snapshot_header() const {HERE;
    std::ostringstream ost;
    ost.precision(16);
    ost << "time" << sep_ << "size" << sep_;
    std::vector<std::string> axes{"x", "y", "z"};
    axes.resize(coords_.front().size());
    wtl::ost_join(ost, axes, sep_) << sep_ << "sites" << sep_ << "fitness\n";
    return ost.str();
}

std::string Tissue::snapshot() const {
    std::ostringstream ost;
    ost.precision(16);
    for (auto& item: tumor_) {
        ost << evolution_history_.size() << sep_ << tumor_.size() << sep_;
        wtl::ost_join(ost, item.first, sep_) << sep_;
        wtl::ost_join(ost, item.second.sites(), "|") << sep_
            << item.second.fitness() << "\n";
    }
    return ost.str();
}

std::string Tissue::mutation_history() const {HERE;
    std::ostringstream ost;
    ost.precision(16);
    std::vector<std::string> xyz{"x", "y", "z"};
    xyz.resize(mutation_coords_.front().size());
    ost << "size" << sep_ << "effect";
    for (auto x: xyz) {
        ost << sep_ << "origin_" << x;
    }
    ost << "\n";
    for (size_t i=0; i<mutation_coords_.size(); ++i) {
        ost << mutation_stages_[i] << sep_
            << Gland::MUTATION_EFFECTS()[i] << sep_;
        wtl::ost_join(ost, mutation_coords_[i], sep_) << "\n";
    }
    return ost.str();
}

//! Stream operator for debug print
std::ostream& operator<< (std::ostream& ost, const Tissue& tissue) {
    return ost << tissue.tumor_;
}

template <class T> inline
void test_coordinate(const std::vector<int>& v) {
    T coord(v);
    std::cerr << typeid(coord).name() << std::endl;
    std::cerr << coord << std::endl;
    std::cerr << coord.distance() << std::endl;
    std::cerr << coord.neighbors() << std::endl;
    std::cerr << T::directions(v.size()) << std::endl;
}

void Tissue::unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr.precision(15);
    Tissue tissue;
    tissue.stain();
    tissue.grow_even(10);
    std::cerr << tissue << std::endl;
    std::cerr << tissue.coords_ << std::endl;
    std::cerr << tissue.snapshot() << std::endl;
    std::cerr << tissue.mutation_history() << std::endl;

    test_coordinate<Orthogonal>({3, -2});
    test_coordinate<Lattice>({3, -2});
    test_coordinate<Hex>({3, -2});

    test_coordinate<Orthogonal>({3, -2, 1});
    test_coordinate<Lattice>({3, -2, 1});
    test_coordinate<Hex>({3, -2, 1});
}

// -*- mode: c++; coding: utf-8 -*-
/*! @file tissue.cpp
    @brief Implementation of Tissue class
*/
#include "tissue.h"

#include <cassert>
#include <iostream>
#include <sstream>

#include <boost/program_options.hpp>

#include <cxxwtils/prandom.hpp>
#include <cxxwtils/iostr.hpp>
#include <cxxwtils/gz.hpp>
#include <cxxwtils/debug.hpp>
#include <cxxwtils/algorithm.hpp>


size_t Tissue::DIMENSIONS_ = 2;
std::string Tissue::SCHEDULE_ = "random";
std::string Tissue::PACKING_ = "push";

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
        ("schedule,S", po::value<std::string>(&SCHEDULE_)->default_value(SCHEDULE_))
        ("packing,P", po::value<std::string>(&PACKING_)->default_value(PACKING_))
    ;
    return desc;
}


void Tissue::stain() {HERE;
    assert(tumor_.empty());
    for (const auto& coord: coord_func_->core()) {
        std::shared_ptr<Gland> founder(new Gland(coord));
        founder->mutate();
        tumor_.insert(founder);
        mutation_coords_.push_back(coord);
        mutation_stages_.push_back(mutation_coords_.size());
    }
}

void Tissue::grow(const size_t max_size) {HERE;
    evolution_history_.reserve(max_size);
    evolution_history_.push_back(snapshot());
    while (tumor_.size() < max_size) {
        std::vector<std::shared_ptr<Gland>> mothers;
        if (SCHEDULE_ == "random") {
            auto it = std::next(tumor_.begin(), wtl::prandom().randrange(tumor_.size()));
            mothers.push_back(*it);
        } else if (SCHEDULE_ == "poisson") {
            mothers.reserve(tumor_.size());
            for (auto it=tumor_.begin(); it!=tumor_.end(); ++it) {
                //! @todo fitness
                const size_t n = wtl::prandom().poisson(1.0);
                for (size_t i=0; i<n; ++i) {
                    mothers.push_back(*it);
                }
            }
        } else if (SCHEDULE_ == "even") {
            mothers.reserve(tumor_.size());
            for (auto it=tumor_.begin(); it!=tumor_.end(); ++it) {
                mothers.push_back(*it);
            }
        } else {exit(1);}
        for (auto& mother: mothers) {
            const auto neighbors = coord_func_->neighbors(mother->coord());
            if (mother->bernoulli_birth()) {
                auto daughter = std::make_shared<Gland>(*mother);  // Gland copy ctor
                bool success = true;
                if (PACKING_ == "push") {push(daughter);}
                else if (PACKING_ == "pushfill") {push_fill(daughter, coord_func_->origin());}
                else if (PACKING_ == "walkfill") {push_fill(daughter);}
                else if (PACKING_ == "empty") {success = fill_empty(daughter);}
                if (success && Gland::bernoulli_mutation()) {
                    daughter->mutate();
                    mutation_coords_.push_back(daughter->coord());
                    mutation_stages_.push_back(tumor_.size());
                }
            }
            if (mother->bernoulli_death()) {
                tumor_.erase(mother);
            } else if (Gland::bernoulli_mutation()) {
                mother->mutate();
                mutation_coords_.push_back(mother->coord());
                mutation_stages_.push_back(tumor_.size());
            }
            if (tumor_.size() <= 128) {
                evolution_history_.push_back(snapshot());
            }
        }
    }
}

std::string Tissue::evolution_history() const {
    return wtl::str_join(evolution_history_, "");
}

void Tissue::push(const std::shared_ptr<Gland>& moving, std::vector<int> direction) {
    if (direction.empty()) {
        const auto directions = coord_func_->directions();
        direction = *wtl::prandom().choice(directions.begin(), directions.end());
    }
    moving->set_coord(moving->coord() + direction);
    const auto result = tumor_.insert(moving);
    if (!result.second) {
        const std::shared_ptr<Gland> existing = *result.first;
        tumor_.erase(result.first);
        push(existing, direction);
        tumor_.insert(moving);
    }
}

void Tissue::push_fill(const std::shared_ptr<Gland>& daughter, const std::vector<int>& direction) {
    const auto present_coord = daughter->coord();
    auto neighbors = coord_func_->neighbors(present_coord);
    std::shuffle(neighbors.begin(), neighbors.end(), wtl::prandom());
//    std::sort(neighbors.begin(), neighbors.end(),
//        [this](const std::vector<int>& x, const std::vector<int>& y) {
//            return coord_func_->euclidean_distance(x) < coord_func_->euclidean_distance(y);
//        });
    for (auto& x: neighbors) {
        daughter->set_coord(x);
        if (tumor_.insert(daughter).second) {
            // end if found and filled an empty space
            return;
        }
    }
    if (direction.empty()) {
        // choose a random direction everytime
        daughter->set_coord(*wtl::prandom().choice(neighbors.begin(), neighbors.end()));
        push_fill(daughter);    
    } else if (direction == coord_func_->origin()) {
        // choose a random direction at first call
        daughter->set_coord(*wtl::prandom().choice(neighbors.begin(), neighbors.end()));
        push_fill(daughter, daughter->coord() - present_coord);
    } else {
        // keep the direction
        daughter->set_coord(present_coord + direction);
        push_fill(daughter, direction);
    }
}

//! @todo
void Tissue::push_layer(const std::shared_ptr<Gland>& daughter) {
}

bool Tissue::fill_empty(const std::shared_ptr<Gland>& daughter) {
    daughter->set_coord(coord_func_->random_neighbor(daughter->coord(), wtl::prandom()));
    const auto result = tumor_.insert(daughter);
    return result.second;
}

std::string Tissue::snapshot_header() const {HERE;
    std::ostringstream ost;
    ost.precision(16);
    ost << "time" << sep_ << "size" << sep_;
    std::vector<std::string> axes{"x", "y", "z"};
    axes.resize(DIMENSIONS_);
    wtl::ost_join(ost, axes, sep_) << sep_ << "sites" << sep_ << "fitness\n";
    return ost.str();
}

std::string Tissue::snapshot() const {
    std::ostringstream ost;
    ost.precision(16);
    for (auto& item: tumor_) {
        ost << evolution_history_.size() << sep_ << tumor_.size() << sep_;
        wtl::ost_join(ost, item->coord(), sep_) << sep_;
        wtl::ost_join(ost, item->sites(), "|") << sep_
            << item->fitness() << "\n";
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
    T coord(v.size());
    std::cerr << typeid(coord).name() << std::endl;
    std::cerr << coord.graph_distance(v) << std::endl;
    std::cerr << coord.euclidean_distance(v) << std::endl;
    std::cerr << coord.neighbors(v) << std::endl;
    std::cerr << coord.directions() << std::endl;
}

void Tissue::unit_test() {
    std::cerr << __PRETTY_FUNCTION__ << std::endl;
    std::cerr.precision(15);

    test_coordinate<Neumann>({3, -2});
    test_coordinate<Moore>({3, -2});
    test_coordinate<Hexagonal>({3, -2});

    test_coordinate<Neumann>({3, -2, 1});
    test_coordinate<Moore>({3, -2, 1});
    test_coordinate<Hexagonal>({3, -2, 1});

    std::cerr << Neumann(2).core() << std::endl;
    std::cerr << Neumann(3).core() << std::endl;
    std::cerr << Hexagonal(2).core() << std::endl;
    std::cerr << Hexagonal(3).core() << std::endl;

    Tissue tissue;
    tissue.stain();
    tissue.grow(10);
    std::cerr << tissue << std::endl;
    std::cerr << tissue.snapshot_header();
    std::cerr << tissue.snapshot() << std::endl;
    std::cerr << tissue.mutation_history() << std::endl;
}

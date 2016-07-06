// -*- mode: c++; coding: utf-8 -*-
/*! @file tissue.h
    @brief Interface of Tissue class
*/
#pragma once
#ifndef TISSUE_HPP_
#define TISSUE_HPP_

#include <iostream>
#include <vector>
#include <valarray>
#include <unordered_set>
#include <memory>
#include <string>

#include <boost/functional/hash.hpp>

#include "cell.hpp"
#include "coord.hpp"

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace boost {
    namespace program_options {
        class options_description;
    }
}

namespace std {
  template <> struct hash<std::valarray<int>> {
    size_t operator() (const std::valarray<int>& v) const {
        return boost::hash_range(std::begin(v), std::end(v));
    }
  };
  template <> struct hash<std::shared_ptr<tumopp::Cell>> {
    size_t operator() (const std::shared_ptr<tumopp::Cell>& x) const {
        return std::hash<std::valarray<int>>()(x->coord());
    }
  };
}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace tumopp {

inline bool all(const std::valarray<bool>& v) {
    return std::all_of(std::begin(v), std::end(v), [](bool x){return x;});
}

class equal_shptr_cell {
  public:
    bool operator() (const std::shared_ptr<tumopp::Cell>& lhs,
                     const std::shared_ptr<tumopp::Cell>& rhs) const {
        return all(lhs->coord() == rhs->coord());
    }
};

class Tissue {
  public:
    static size_t DIMENSIONS() {return DIMENSIONS_;}

    //! Constructor
    Tissue() = default;

    //! main function
    bool grow(const size_t max_size=MAX_SIZE_);

    std::ostream& write_segsites(std::ostream&, const std::vector<std::shared_ptr<Cell>>&) const;
    std::vector<std::shared_ptr<Cell>> sample_random(const size_t) const;
    std::vector<std::shared_ptr<Cell>> sample_section(const size_t) const;

    void collect(std::ostream&, const Cell&);
    void snap(std::ostream&);

    std::string specimens() const {return specimens_.str();}
    std::string snapshots() const {return snapshots_.str();}
    std::string drivers() const {return drivers_.str();}
    std::string pairwise_distance(const size_t) const;

    friend std::ostream& operator<< (std::ostream&, const Tissue&);

    size_t size() const {return tumor_.size();};
    double radius() const {return coord_func_->radius(tumor_.size());}
    double relative_pos(std::shared_ptr<Cell> x) const {
        return coord_func_->euclidean_distance(x->coord()) / radius();
    }

    //! Set coordinate function object
    template <class FuncObj>
    void set_coord() {coord_func_ = std::make_unique<FuncObj>(DIMENSIONS_);}
    const std::unique_ptr<Coord>& coord_func() const {return coord_func_;}

    //! Unit test
    static void unit_test();
    static boost::program_options::options_description opt_description();

  private:
    //! Dimensions: {1, 2, 3}
    static size_t DIMENSIONS_;

    //! Coordinate system
    static std::string COORDINATE_;

    //! packing method: push, push_fill, walk_fill
    static std::string PACKING_;

    //! 0: flat, +: peripheral growth
    static double SIGMA_E_;

    //! initial population size
    static size_t INITIAL_SIZE_;

    //! Max tumor size to stop simulation
    static size_t MAX_SIZE_;

    void init();

    bool insert(const std::shared_ptr<Cell>&);

    //! Swap with a random neighbor
    void migrate(const std::shared_ptr<Cell>&);

    //! Emplace daughter cell and push other cells outward
    void push(std::shared_ptr<Cell> moving, const std::valarray<int>& direction);
    void pushn_everytime(std::shared_ptr<Cell> moving);
    //! Fill empty neighbor or push to the direction
    void fill_push(std::shared_ptr<Cell> moving, const std::valarray<int>& direction);
    //! Insert x if it has an empty neighbor
    bool fill_empty(const std::shared_ptr<Cell>& x);
    //! Put new cell and return existing.
    bool swap_existing(std::shared_ptr<Cell>* x);
    size_t steps_to_empty(std::valarray<int> current, const std::valarray<int>& direction) const;
    std::valarray<int> to_nearest_empty(const std::valarray<int>& current, size_t search_max=26) const;

    size_t num_empty_neighbors(const std::valarray<int>&) const;

    void queue_push(double delta_t, const std::shared_ptr<Cell>&);

    double positional_value(const std::valarray<int>&) const;

    std::vector<size_t> generate_neutral_mutations() const;

    std::string header() const;

    /////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////
    // Data member

    //! cells
    std::unordered_set<
        std::shared_ptr<Cell>,
        std::hash<std::shared_ptr<Cell>>,
        equal_shptr_cell> tumor_;

    //! event queue
    std::multimap<double, std::shared_ptr<Cell>> queue_;

    double time_ = 0.0;
    size_t id_tail_ = 0;

    std::unique_ptr<Coord> coord_func_;

    std::ostringstream specimens_;
    std::ostringstream snapshots_;
    std::ostringstream drivers_;
    const char* sep_ = "\t";
};

} // namespace tumopp

#endif /* TISSUE_HPP_ */

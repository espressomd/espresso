#ifndef PROPAGATION_HPP
#define PROPAGATION_HPP

#include <boost/iterator/filter_iterator.hpp>

template <typename Predicate> struct PropagationPredicate {
  Predicate predicate;

  PropagationPredicate(Predicate pred) : predicate(pred) {}

  bool operator()(Particle const &p) { return predicate(p.propagation()); };
};

template <typename Predicate>
class ParticleRangeFiltered
    : public boost::iterator_range<boost::iterators::filter_iterator<
          PropagationPredicate<Predicate>, ParticleIterator<Cell **>>> {
  using base_type = boost::iterator_range<boost::iterators::filter_iterator<
      PropagationPredicate<Predicate>, ParticleIterator<Cell **>>>;

public:
  using base_type::base_type;
  auto size() const { return std::distance(this->begin(), this->end()); };
};

/*
using ParticleRangeDefault =
    ParticleRangeFiltered<PropagationMode::TRANS_SYSTEM_DEFAULT>;
using ParticleRangeLangevin =
    ParticleRangeFiltered<PropagationMode::TRANS_LANGEVIN>;
using ParticleRangeStokesian =
    ParticleRangeFiltered<PropagationMode::TRANS_STOKESIAN>;

*/
#endif
/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#define BOOST_TEST_MODULE test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "field_coupling/ForceField.hpp"
#include "field_coupling/PotentialField.hpp"

#include "field_coupling/detail/Base.hpp"
#include "field_coupling/detail/BindCoupling.hpp"
using namespace FieldCoupling;

#include <utils/Vector.hpp>

template <bool linear> struct Id {
  static constexpr const bool is_linear = linear;
  mutable int count = 0;

  template <class P, class T> T operator()(const P &, T x) const {
    count++;
    return x;
  }
};

struct Particle {};

BOOST_AUTO_TEST_CASE(BindCoupling_test) {
  using FieldCoupling::detail::BindCoupling;
  using FieldCoupling::detail::make_bind_coupling;

  /* is_linear */
  {
    static_assert(BindCoupling<Id<true>, Particle>::is_linear, "");
    static_assert(!BindCoupling<Id<false>, Particle>::is_linear, "");
  }

  /* make_bind_coupling */
  {
    /* Return type */
    using std::is_same;
    using return_type = decltype(make_bind_coupling(int{}, float{}));
    static_assert(is_same<BindCoupling<int, float>, return_type>::value, "");
  }

  /* Call */
  {
    auto const id = Id<true>{};
    auto const bc = make_bind_coupling(id, Particle{});
    const int x = 5;

    BOOST_CHECK(5 == bc(5));
    BOOST_CHECK(id.count == 1);
  }
}

struct MoveOnly {
  const int id;

  MoveOnly(int id) : id(id) {}
  MoveOnly(MoveOnly const &) = delete;
  MoveOnly(MoveOnly &&) = default;
};

BOOST_AUTO_TEST_CASE(FieldBase_test) {
  using FieldCoupling::detail::Base;

  /* ctor move */
  {
    auto base = Base<MoveOnly, MoveOnly>(MoveOnly{1}, MoveOnly{2});

    BOOST_CHECK(1 == base.coupling().id);
    BOOST_CHECK(2 == base.field().id);
  }

  /* ctor copy */
  {
    int c = 1;
    int f = 2;

    auto base = Base<int, int>(c, f);

    BOOST_CHECK(1 == base.coupling());
    BOOST_CHECK(2 == base.field());
  }
}

struct DummyVectorField {
  Utils::Vector3d operator()(const Utils::Vector3d &x, double t) const {
    return t * x;
  }
};

BOOST_AUTO_TEST_CASE(ForceField_test) {
  auto ff =
      ForceField<Id<true>, DummyVectorField>(Id<true>{}, DummyVectorField{});
  const Utils::Vector3d x{1., 2., 3.};
  const int p = 5;

  BOOST_CHECK((9. * x) == ff.force(5, x, 9.));
  BOOST_CHECK(1 == ff.coupling().count);
}

struct DummyScalarField {
  double operator()(const Utils::Vector3d &x, double t) const {
    return t * x.norm();
  }
  Utils::Vector3d jacobian(const Utils::Vector3d &x, double = {}) const {
    return 3. * x;
  }
};

BOOST_AUTO_TEST_CASE(PotentialField_test) {
  auto pf = PotentialField<Id<true>, DummyScalarField>(Id<true>{},
                                                       DummyScalarField{});
  const Utils::Vector3d x{1., 2., 3.};

  BOOST_CHECK((2. * x.norm()) == pf.energy(5, x, 2.));
  BOOST_CHECK(1 == pf.coupling().count);

  BOOST_CHECK(-(3. * x) == pf.force(5, x, 0));
  BOOST_CHECK(2 == pf.coupling().count);
}

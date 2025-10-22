/*
 * Copyright (C) 2010-2022 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define BOOST_TEST_MODULE Field coupling test for force fields
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "field_coupling/ForceField.hpp"
#include "field_coupling/PotentialField.hpp"

#include "field_coupling/detail/Base.hpp"
#include "field_coupling/detail/BindCoupling.hpp"

#include <utils/Vector.hpp>

#include <type_traits>

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
    static_assert(BindCoupling<Id<true>, Particle>::is_linear);
    static_assert(!BindCoupling<Id<false>, Particle>::is_linear);
  }

  /* make_bind_coupling */
  {
    /* Return type */
    using return_type = decltype(make_bind_coupling(int{}, float{}));
    static_assert(std::is_same_v<BindCoupling<int, float>, return_type>);
  }

  /* Call */
  {
    auto const id = Id<true>{};
    auto const bc = make_bind_coupling(id, Particle{});
    auto const x = 5;

    BOOST_CHECK_EQUAL(id.count, 0);
    BOOST_CHECK_EQUAL(bc(5), 5);
    BOOST_CHECK_EQUAL(id.count, 1);
    BOOST_CHECK_EQUAL(bc(x), x);
    BOOST_CHECK_EQUAL(id.count, 2);
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
    auto const c = 1;
    auto const f = 2;

    auto base = Base<int, int>(c, f);

    BOOST_CHECK_EQUAL(base.coupling(), c);
    BOOST_CHECK_EQUAL(base.field(), f);
  }
}

struct DummyVectorField {
  Utils::Vector3d operator()(const Utils::Vector3d &x, double t) const {
    return t * x;
  }
};

BOOST_AUTO_TEST_CASE(ForceField_test) {
  using FieldCoupling::ForceField;
  auto ff =
      ForceField<Id<true>, DummyVectorField>(Id<true>{}, DummyVectorField{});
  auto const x = Utils::Vector3d{1., 2., 3.};

  BOOST_CHECK_EQUAL(ff.coupling().count, 0);
  BOOST_CHECK_EQUAL(ff.force(5, x, 9.), 9. * x);
  BOOST_CHECK_EQUAL(ff.coupling().count, 1);
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
  using FieldCoupling::PotentialField;
  auto pf = PotentialField<Id<true>, DummyScalarField>(Id<true>{},
                                                       DummyScalarField{});
  auto const x = Utils::Vector3d{1., 2., 3.};
  BOOST_CHECK_EQUAL(pf.coupling().count, 0);

  BOOST_CHECK_EQUAL(pf.energy(5, x, 2.), 2. * x.norm());
  BOOST_CHECK_EQUAL(pf.coupling().count, 1);

  BOOST_CHECK_EQUAL(pf.force(5, x, 0), -3. * x);
  BOOST_CHECK_EQUAL(pf.coupling().count, 2);
}

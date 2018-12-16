#include <iostream>
#include <vector>

#include "device_matrix.hpp"
#include "multi_array.hpp"

#ifdef __CUDACC__
#define DEVICE_FUNC __host__ __device__
#else
#define DEVICE_FUNC
#endif

namespace sd {

// Compute distance between particle pairs and check for overlaps
template <typename Policy, typename T>
struct check_dist {
    device_vector_view<T, Policy> const x;
    device_matrix_view<T, Policy> pd;
    device_matrix_view<std::size_t, Policy> const part_id;
    DEVICE_FUNC void operator()(std::size_t i) {
        std::size_t k = 6 * part_id(0, i);
        std::size_t j = 6 * part_id(1, i);

        T dx = x(j + 0) - x(k + 0);
        T dy = x(j + 1) - x(k + 1);
        T dz = x(j + 2) - x(k + 2);
        T dr = std::sqrt(dx * dx + dy * dy + dz * dz);
        T dr_inv = 1 / dr;

#ifdef NDEBUG
#undef NDEBUG
        assert(dr > 2 && "Particles overlapped!");
#define NDEBUG
#else
        assert(dr > 2 && "Particles overlapped!");
#endif

        pd(0, i) = dx * dr_inv;
        pd(1, i) = dy * dr_inv;
        pd(2, i) = dz * dr_inv;
        pd(3, i) = dr_inv;
    }
};

// The expression for the mobility matrix is given in Appendix A of
//
//   Durlofsky, L. and Brady, J.F. and Bossis G., J . Fluid Mech. 180, 21-49
//   (1987) https://authors.library.caltech.edu/32068/1/DURjfm87.pdf
//
// The submatrices mob_a, mob_b, mob_c, mob_gt, mob_ht, mob_m are
// given by equation (A 2).
//
// Alternatively, the same expressions are given in
//
//   Kim, S. and Mifflin, R.T., Physics of Fluids 28, 2033 (1985)
//   https://doi.org/10.1063/1.865384
template <typename Policy, typename T, bool isself>
struct mobility;

template <typename Policy, typename T>
struct mobility<Policy, T, true> {
    device_matrix_view<T, Policy> zmuf, zmus, zmes;

    // Determine the self contribution
    // This is independent of dr_inv, dx, dy, dz
    DEVICE_FUNC void operator()(std::size_t part_id) {
        static constexpr multi_array<T, 3, 3> const mob_a = {
            // clang-format off
            1, 0, 0, 
            0, 1, 0, 
            0, 0, 1
            // clang-format on
        };
        static constexpr multi_array<T, 3, 3> const mob_c = {
            // clang-format off
            3./4., 0    , 0    ,
            0    , 3./4., 0    , 
            0    , 0    , 3./4.
            // clang-format on
        };
        static constexpr multi_array<T, 5, 5> const mob_m = {
            // clang-format off
            9./5. , 0    , 0    , 0    , 9./10., 
            0     , 9./5., 0    , 0    , 0     , 
            0     , 0    , 9./5., 0    , 0     , 
            0     , 0    , 0    , 9./5., 0     , 
            9./10., 0    , 0    , 0    , 9./5.
            // clang-format on
        };

        // Fill the self mobility terms
        std::size_t ph1 = 6 * part_id;
        std::size_t ph2 = ph1 + 3;
        std::size_t ph3 = 5 * part_id;

        for (std::size_t i = 0; i < 3; ++i) {
            zmuf(ph1 + i, ph1 + i) = mob_a(i, i);
            zmuf(ph2 + i, ph2 + i) = mob_c(i, i);
        }

        for (std::size_t i = 0; i < 5; ++i) {
            for (std::size_t j = 0; j < 5; ++j) {
                zmes(ph3 + i, ph3 + j) = mob_m(i, j);
            }
        }
    }
};

template <typename Policy, typename T>
struct mobility<Policy, T, false> {
    device_matrix_view<T, Policy> zmuf, zmus, zmes;
    device_matrix_view<T, Policy> const pd;
    device_matrix_view<std::size_t, Policy> const part_id;

    // Determine the pair contribution
    DEVICE_FUNC void operator()(std::size_t pair_id) {
        static constexpr multi_array<T, 3, 3> const delta = {
            // clang-format off
            1, 0, 0,
            0, 1, 0,
            0, 0, 1
            // clang-format on
        };

        static constexpr multi_array<T, 3, 3, 3> const eps = {
            // clang-format off
            0, 0, 0,   0, 0, 1,   0,-1, 0,
            0, 0,-1,   0, 0, 0,   1, 0, 0,
            0, 1, 0,  -1, 0, 0,   0, 0, 0
            // clang-format on
        };

        static constexpr multi_array<std::size_t, 2, 5> const mesid = {
            // clang-format off
            0, 0, 0, 1, 1,
            2, 1, 2, 2, 2
            // clang-format on
        };

        T dx = pd(0, pair_id);
        T dy = pd(1, pair_id);
        T dz = pd(2, pair_id);
        T dr_inv = pd(3, pair_id);

        multi_array<T, 3> e = {dx, dy, dz};
        auto ee = outer(e, e);

        T dr_inv2 = dr_inv * dr_inv;
        T dr_inv3 = dr_inv2 * dr_inv;
        T dr_inv4 = dr_inv3 * dr_inv;
        T dr_inv5 = dr_inv4 * dr_inv;

        // The following scalar mobility functions can be found in
        // equation (A 3).
        T x12a = T{3. / 2} * dr_inv - dr_inv3;
        T y12a = T{3. / 4.} * dr_inv + T{1. / 2.} * dr_inv3;

        T y12b = T{-3. / 4.} * dr_inv2;

        T x12c = T{3. / 4.} * dr_inv3;
        T y12c = T{-3. / 8.} * dr_inv3;

        T x12g = T{9. / 4.} * dr_inv2 - T{18. / 5.} * dr_inv4;
        T y12g = T{6. / 5.} * dr_inv4;

        T y12h = T{-1.0 * 9. / 8.} * dr_inv3;

        T x12m = T{-1.0 * 9. / 2.} * dr_inv3 + T{54. / 5.} * dr_inv5;
        T y12m = T{9. / 4.} * dr_inv3 - T{36. / 5.} * dr_inv5;
        T z12m = T{9. / 5.} * dr_inv5;

        // Equation (A 2) first, second, and third line
        multi_array<T, 3, 3> mob_a;
        multi_array<T, 3, 3> mob_b;
        multi_array<T, 3, 3> mob_c;
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                std::size_t k = (3 - i - j) % 3;
                mob_a(i, j) = x12a * ee(i, j) + y12a * (delta(i, j) - ee(i, j));
                mob_b(i, j) = y12b * eps(i, j, k) * e(k);
                mob_c(i, j) = x12c * ee(i, j) + y12c * (delta(i, j) - ee(i, j));
            }
        }

        // Equation (A 2) fourth and fifth line
        multi_array<T, 3, 3, 3> gt;
        multi_array<T, 3, 3, 3> ht;
        for (std::size_t k = 0; k < 3; ++k) {
            for (std::size_t i = 0; i < 3; ++i) {
                for (std::size_t j = 0; j < 3; ++j) {
                    gt(k, i, j) =
                        -(x12g * (ee(i, j) - T{1. / 3.} * delta(i, j)) * e(k) +
                          y12g * (e(i) * delta(j, k) + e(j) * delta(i, k) -
                                  T{2.0} * ee(i, j) * e(k)));

                    std::size_t l = (3 - j - k) % 3;
                    std::size_t m = (3 - i - k) % 3;
                    ht(k, i, j) = y12h * (ee(i, l) * eps(j, k, l) +
                                          ee(j, m) * eps(i, k, m));
                }
            }
        }

        // Equation (A 2) sixth line
        multi_array<T, 3, 3, 3, 3> m;
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                for (std::size_t k = 0; k < 3; ++k) {
                    for (std::size_t l = 0; l < 3; ++l) {
                        m(i, j, k, l) =
                            T{3. / 2.} * x12m *
                                (ee(i, j) - T{1. / 3.} * delta(i, j)) *
                                (ee(k, l) - T{1. / 3.} * delta(k, l)) +
                            T{1. / 2.} * y12m *
                                (ee(i, k) * delta(j, l) +
                                 ee(j, k) * delta(i, l) +
                                 ee(i, l) * delta(j, k) +
                                 ee(j, l) * delta(i, k) -
                                 T{4.0} * ee(i, j) * ee(k, l)) +
                            T{1. / 2.} * z12m *
                                (delta(i, k) * delta(j, l) +
                                 delta(j, k) * delta(i, l) -
                                 delta(i, j) * delta(k, l) +
                                 ee(i, j) * delta(k, l) +
                                 ee(k, l) * delta(i, j) -
                                 ee(i, k) * delta(j, l) -
                                 ee(j, k) * delta(i, l) -
                                 ee(i, l) * delta(j, k) -
                                 ee(j, l) * delta(i, k) + ee(i, j) * ee(k, l));
                    }
                }
            }
        }

        // This segment of code converts the pair contributions to the grand
        // mobility tensor into a symmetric matrix.  Using the following
        // conversions of the shear rate E and the stresslet S into the vectors
        // EV and SV respectively. EV_1 = E_11 - E_33, EV_2 = 2 E_12, EV_3 = 2
        // E_13, EV_4 = 2 E_23, EV_5 = E_22 - E_33 SV_1 = S_11, SV_2 = S_12 =
        // S_21, SV_3 = S_13 = S_31, SV_4 = S_23 = S_32, SV_5 = S_22
        multi_array<T, 3, 5> mob_gt;
        multi_array<T, 3, 5> mob_ht;
        for (std::size_t i = 0; i < 3; ++i) {
            mob_gt(i, 0) = gt(i, 0, 0) - gt(i, 2, 2);
            mob_gt(i, 1) = T{2.0} * gt(i, 0, 1);
            mob_gt(i, 2) = T{2.0} * gt(i, 0, 2);
            mob_gt(i, 3) = T{2.0} * gt(i, 1, 2);
            mob_gt(i, 4) = gt(i, 1, 1) - gt(i, 2, 2);

            mob_ht(i, 0) = ht(i, 0, 0) - ht(i, 2, 2);
            mob_ht(i, 1) = T{2.0} * ht(i, 0, 1);
            mob_ht(i, 2) = T{2.0} * ht(i, 0, 2);
            mob_ht(i, 3) = T{2.0} * ht(i, 1, 2);
            mob_ht(i, 4) = ht(i, 1, 1) - ht(i, 2, 2);
        }

        multi_array<T, 5, 5> mob_m;
        for (std::size_t i = 0; i < 5; ++i) {
            if (i == 0 || i == 4) {
                mob_m(i, 0) = m(mesid(0, i), mesid(0, i), 0, 0) -
                              m(mesid(0, i), mesid(0, i), 2, 2) -
                              (m(mesid(1, i), mesid(1, i), 0, 0) -
                               m(mesid(1, i), mesid(1, i), 2, 2));
                mob_m(i, 1) = T{2.0} * (m(mesid(0, i), mesid(0, i), 0, 1) -
                                        m(mesid(1, i), mesid(1, i), 0, 1));
                mob_m(i, 2) = T{2.0} * (m(mesid(0, i), mesid(0, i), 0, 2) -
                                        m(mesid(1, i), mesid(1, i), 0, 2));
                mob_m(i, 3) = T{2.0} * (m(mesid(0, i), mesid(0, i), 1, 2) -
                                        m(mesid(1, i), mesid(1, i), 1, 2));
                mob_m(i, 4) = m(mesid(0, i), mesid(0, i), 1, 1) -
                              m(mesid(0, i), mesid(0, i), 2, 2) -
                              (m(mesid(1, i), mesid(1, i), 1, 1) -
                               m(mesid(1, i), mesid(1, i), 2, 2));
            } else {
                mob_m(i, 0) = T{2.0} * (m(mesid(0, i), mesid(1, i), 0, 0) -
                                        m(mesid(0, i), mesid(1, i), 2, 2));
                mob_m(i, 1) = T{4.0} * m(mesid(0, i), mesid(1, i), 0, 1);
                mob_m(i, 2) = T{4.0} * m(mesid(0, i), mesid(1, i), 0, 2);
                mob_m(i, 3) = T{4.0} * m(mesid(0, i), mesid(1, i), 1, 2);
                mob_m(i, 4) = T{2.0} * (m(mesid(0, i), mesid(1, i), 1, 1) -
                                        m(mesid(0, i), mesid(1, i), 2, 2));
            }
        }

        // Fill the pair mobility terms
        std::size_t ph1 = part_id(0, pair_id);
        std::size_t ph2 = part_id(1, pair_id);

        std::size_t ph5 = 5 * ph1;
        std::size_t ph6 = 5 * ph2;

        ph1 = 6 * ph1;
        ph2 = 6 * ph2;

        std::size_t ph3 = ph1 + 3;
        std::size_t ph4 = ph2 + 3;

        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                zmuf(ph1 + i, ph2 + j) = mob_a(i, j);
                zmuf(ph3 + i, ph2 + j) = mob_b(i, j);
                zmuf(ph1 + i, ph4 + j) = -mob_b(j, i); // mob_b transpose
                zmuf(ph3 + i, ph4 + j) = mob_c(i, j);

                zmuf(ph2 + i, ph1 + j) = mob_a(j, i);
                zmuf(ph4 + i, ph1 + j) = mob_b(j, i);
                zmuf(ph2 + i, ph3 + j) = -mob_b(i, j); // mob_b transpose
                zmuf(ph4 + i, ph3 + j) = mob_c(j, i);
            }

            for (std::size_t j = 0; j < 5; ++j) {
                zmus(ph1 + i, ph6 + j) = mob_gt(i, j);
                zmus(ph2 + i, ph5 + j) = T{-1.0} * mob_gt(i, j);

                zmus(ph3 + i, ph6 + j) = mob_ht(i, j);
                zmus(ph4 + i, ph5 + j) = mob_ht(i, j);
            }
        }

        for (std::size_t i = 0; i < 5; ++i) {
            for (std::size_t j = 0; j < 5; ++j) {
                zmes(ph5 + i, ph6 + j) = mob_m(i, j);
                zmes(ph6 + i, ph5 + j) = mob_m(j, i);
            }
        }
    }
};

template <typename Policy, typename T>
struct solver {
    static_assert(policy::is_policy<Policy>::value,
                  "The execution policy must meet the requirements");
    static_assert(std::is_arithmetic<T>::value,
                  "Data type of device_matrix must be arithmetic for "
                  "arithmetic operations");

    template <typename U>
    using vector_type = typename Policy::template vector<U>;

    std::size_t const n_part;
    std::size_t const n_pair;

    device_matrix<T, Policy> zmuf;
    device_matrix<T, Policy> zmus;
    device_matrix<T, Policy> zmes;

    solver(std::size_t const n_part)
        : n_part(n_part), n_pair(n_part * (n_part - 1) / 2),
          zmuf(n_part * 6, n_part * 6), zmus(n_part * 6, n_part * 5),
          zmes(n_part * 5, n_part * 5) {}

    std::vector<T> calc_vel(std::vector<T> const &x_host,
                            std::vector<T> const &f_host) {
        assert(x_host.size() == 6 * n_part);
        vector_type<T> x(x_host.begin(), x_host.end());

        // TODO: More efficient!
        device_matrix<std::size_t, Policy> part_id(2, n_pair);
        std::size_t k = 0;
        for (std::size_t i = 0; i < n_part; ++i) {
            for (std::size_t j = i + 1; j < n_part; ++j) {
                part_id(0, k) = i;
                part_id(1, k) = j;
                k += 1;
            }
        }

        device_matrix<T, Policy> pd(4, n_pair);
        thrust::counting_iterator<std::size_t> begin(0UL);

        // check_dist
        thrust::for_each(Policy::par(), begin, begin + n_pair,
                         check_dist<Policy, T>{x, pd, part_id});

        // self mobility term
        thrust::for_each(Policy::par(), begin, begin + n_part,
                         mobility<Policy, T, true>{zmuf, zmus, zmes});

        // pair mobility term
        thrust::for_each(
            Policy::par(), begin, begin + n_pair,
            mobility<Policy, T, false>{zmuf, zmus, zmes, pd, part_id});

        // Invert the grand-mobility tensor.  This is done in several steps
        // which minimize the computation time.

        // Invert R1 = Muf ^ -1 => zmuf = zmuf ^ -1
        zmuf = zmuf.inverse();

        // if (mode == FTS) {

        // Compute R2 = Mus(t) * R1 => rsu = zmus(t) * zmuf
        device_matrix<T, Policy> const &rsu = zmus.transpose() * zmuf;

        // Compute R3 = R2 * Mus - Mes => zmes = rsu * zmus - zmes
        zmes = zmes - rsu * zmus;

        // Invert  R4 = R3 ^ -1 => zmes = zmes ^ -1
        zmes = zmes.inverse();

        // Compute R5 = -R3 * R4 => zmus = -rsu(t) * zmes
        zmus = -(rsu.transpose() * zmes);

        // Compute R6 = R1 - R5  => zmuf = zmuf - zmus * rsu
        zmuf = zmuf - zmus * rsu;

        //}

        device_matrix<T, Policy> const &rfu = zmuf;
        device_matrix<T, Policy> const &rfe = zmus;

        // TODO: COMPUTE LUBRICATION CORRECTION HERE

        assert(f_host.size() == 6 * n_part);
        vector_type<T> fext(f_host.begin(), f_host.end());
        vector_type<T> uinf(n_part * 6, T{0.0});
        vector_type<T> einf(n_part * 5, T{0.0});
        vector_type<T> u = rfu.inverse() * (fext + rfe * einf) + uinf;

        // return the change in velocity due to HI
        std::vector<T> out(u.size());
        thrust::copy(u.begin(), u.end(), out.begin());
        return out;
    }
};

} // namespace sd

/*
 * Copyright (C) 2018-2020 The ESPResSo project
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

/** @file
 *  In this file, various types are provided which are useful for matrix
 *  arithmetics on a parallel computing "device". Various BLAS and LAPACK
 *  routines are accessed for efficient matrix calculations.
 */
 
#ifndef SD_DEVICE_MATRIX_HPP
#define SD_DEVICE_MATRIX_HPP

#include <iomanip>
#include <limits>
#include <type_traits>

#include "thrust_wrapper.hpp"

#ifdef __CUDACC__
#include <cublas_v2.h>
#include <cusolverDn.h>
#endif

#ifdef __CUDACC__
#define DEVICE_FUNC __host__ __device__
#else
#define DEVICE_FUNC
#endif

#ifdef __GNUC__
#define MAYBE_UNUSED __attribute__((unused))
#else
#define MAYBE_UNUSED
#endif

/** The `host` and `device` structs are used to distinguish between routines
 *  that are executed on the Thrust "host" and on the Thrust "device". These
 *  types are passed to certain routines within this file as a template
 *  parameter. Within those routines, its `vector` type denotes either
 *  `host_vector` or `device_vector` depending on whether the routine is called
 *  for host or device. The `par()` function returns the according Thrust
 *  execution policy that is needed when Thrust routines are called.
 */
namespace policy {

struct host {
    template <typename T>
    using vector = thrust_wrapper::host_vector<T>;
    // the remove_const is necessary, because things cannot be static and const
    // at the same time
    static std::remove_const<decltype(thrust_wrapper::host)>::type par()
      { return thrust_wrapper::host; }
};

struct device {
    template <typename T>
    using vector = thrust_wrapper::device_vector<T>;
    static std::remove_const<decltype(thrust_wrapper::device)>::type par()
      { return thrust_wrapper::device; }
};

template <typename...>
using void_t = void;

template <typename Policy, typename = void>
struct is_policy : std::false_type {};

template <typename Policy>
struct is_policy<Policy, void_t<typename Policy::template vector<class T>,
                                decltype(Policy::par())>> : std::true_type {};

} // namespace policy

DEVICE_FUNC inline thrust_wrapper::tuple<std::size_t, std::size_t>
unravel_index(std::size_t index, std::size_t lda) {
    std::size_t i, j;
    i = index % lda;
    index /= lda;
    j = index;
    return thrust_wrapper::make_tuple(i, j);
    ;
}

/// \cond

// LAPACK prototypes
#ifndef __CUDACC__
extern "C" {
int dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha,
           double *a, int *lda, double *b, int *ldb, double *beta, double *c,
           int *ldc);

int dgemv_(char *trans, int *m, int *n, double *alpha, double *a, int *lda,
           double *x, int *incx, double *beta, double *y, int *incy);

int dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);

int dpotrs_(char *uplo, int *n, int *nrhs, double *a, int *lda, double *b,
            int *ldb, int *info);
}
#endif

namespace internal {

/** The `cublas` struct channels access to efficient basic matrix operations
 *  provided by either the cuBLAS library (after which it is named), which
 *  executes on the GPU, or by BLAS (Basic Linear Algebra Subprograms), which
 *  executes on the CPU.
 *  \tparam Policy Specifies whether the routines are executed on host or on
 *                 device, by either passing \ref policy::host or \ref
 *                 policy::device
 *  \tparam T      Specifies the data type (only `double` is available).
 */
template <typename Policy, typename T>
struct cublas {};

#ifdef __CUDACC__
/** Basic matrix operations on device (GPU)
 */
template <>
struct cublas<policy::device, double> {
    /** Transpose matrix
     *  \param A buffer on device for the matrix to be transposed.
     *  \param C buffer on device to store the result. May be the same as A.
     *  \param m number of rows of A
     *  \param n number of columns of A
     */
    static void geam(double const *A, double *C, int m, int n) {
        double const alpha = 1;
        double const beta = 0;

        MAYBE_UNUSED cublasStatus_t stat;
        cublasHandle_t handle;
        stat = cublasCreate(&handle);
        assert(CUBLAS_STATUS_SUCCESS == stat);

        stat = cublasDgeam(handle, CUBLAS_OP_T, CUBLAS_OP_T, n, m, &alpha, A, m,
                           &beta, A, m, C, n);
        assert(CUBLAS_STATUS_SUCCESS == stat);

        cublasDestroy(handle);
    }

    /** Matrix matrix multiplication
     *  \param A buffer on device for the matrix: first factor
     *  \param B buffer on device for the matrix: second factor
     *  \param C buffer on device to store the result
     *  \param m number of rows of A
     *  \param n number of columns of B
     *  \param k number of columns of A, rows of B
     */
    static void gemm(const double *A, const double *B, double *C, int m, int k,
                     int n) {
        int lda = m, ldb = k, ldc = m;
        double alpha = 1;
        double beta = 0;

        MAYBE_UNUSED cublasStatus_t stat;
        cublasHandle_t handle;
        stat = cublasCreate(&handle);
        assert(CUBLAS_STATUS_SUCCESS == stat);

        stat = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, A,
                           lda, B, ldb, &beta, C, ldc);
        assert(CUBLAS_STATUS_SUCCESS == stat);

        cublasDestroy(handle);
    }

    /** Matrix vector multiplication
     *  \param A buffer on device for the matrix
     *  \param x buffer on device for the vector
     *  \param y buffer on device for the result
     *  \param m number of rows of A
     *  \param n number of columns of A
     */
    static void gemv(const double *A, const double *x, double *y, int m,
                     int n) {
        int lda = m;
        double alpha = 1;
        double beta = 0;
        int incx = 1;
        int incy = 1;

        MAYBE_UNUSED cublasStatus_t stat;
        cublasHandle_t handle;
        stat = cublasCreate(&handle);
        assert(CUBLAS_STATUS_SUCCESS == stat);

        stat = cublasDgemv(handle, CUBLAS_OP_N, m, n, &alpha, A, lda, x, incx,
                           &beta, y, incy);
        assert(CUBLAS_STATUS_SUCCESS == stat);

        cublasDestroy(handle);
    }
};
#else
template <>
struct cublas<policy::host, double> {
    /** Transpose matrix
     *  \param A buffer on host for the matrix to be transposed.
     *  \param C buffer on host to store the result. May be the same as A.
     *  \param m number of rows of A
     *  \param n number of columns of A
     */
    static void geam(double const *A, double *C, int m, int n) {
        // m = m_rows, n = m_cols
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                // (row,col) = row + col * m_rows
                C[j + i * n] = A[i + j * m];
            }
        }
    }

    /** Matrix matrix multiplication
     *  \param A buffer on host for the matrix: first factor
     *  \param B buffer on host for the matrix: second factor
     *  \param C buffer on host to store the result
     *  \param m number of rows of A
     *  \param n number of columns of B
     *  \param k number of columns of A, rows of B
     */
    static void gemm(const double *A, const double *B, double *C, int m, int k,
                     int n) {
        int lda = m, ldb = k, ldc = m;
        double alpha = 1;
        double beta = 0;

        char N = 'N';
        dgemm_(&N, &N, &m, &n, &k, &alpha, const_cast<double *>(A), &lda,
               const_cast<double *>(B), &ldb, &beta, C, &ldc);
    }

    /** Matrix vector multiplication
     *  \param A buffer on host for the matrix
     *  \param x buffer on host for the vector
     *  \param y buffer on host for the result
     *  \param m number of rows of A
     *  \param n number of columns of A
     */
    static void gemv(const double *A, const double *x, double *y, int m,
                     int n) {
        int lda = m;
        double alpha = 1;
        double beta = 0;
        int incx = 1;
        int incy = 1;

        char N = 'N';
        dgemv_(&N, &m, &n, &alpha, const_cast<double *>(A), &lda,
               const_cast<double *>(x), &incx, &beta, y, &incy);
    }
};
#endif

/** The `cusolver` struct channels access to efficient linear algebra solvers
 *  provided by either the cuSOLVER library (after which it is named), which
 *  executes on the GPU, or by LAPACK (Linear Algebra PACKage), which executes
 *  on the CPU.
 *  The first template parameter specifies whether the routines are executed on
 *  host or on device, by either passing `policy::host` or `policy::device`.
 *  The second template parameter specifies the data type (only `double` is 
 *  available).
 */
template <typename, typename>
struct cusolver;

#ifdef __CUDACC__
template <>
struct cusolver<policy::device, double> {
    /** Computes the Cholesky factorization and the inverse of a real symmetric
     *  positive definite matrix, looking only in the top half of the symmetric
     *  matrix.
     *  
     *  \param A buffer on device for the symmetric input matrix,
     *           serves as output for the Cholesky decomposition
     *  \param B buffer on device expects identity matrix,
     *           serves as output for the inverse
     *  \param N size of the matrix
     */
    static void potrf(double *A, double *B, int N) {
        MAYBE_UNUSED cusolverStatus_t stat;
        cusolverDnHandle_t handle;
        stat = cusolverDnCreate(&handle);
        assert(CUSOLVER_STATUS_SUCCESS == stat);

        int lwork = -1;
        stat = cusolverDnDpotrf_bufferSize(handle, CUBLAS_FILL_MODE_UPPER, N, A,
                                           N, &lwork);
        assert(CUSOLVER_STATUS_SUCCESS == stat);

        assert(lwork != -1);

        thrust_wrapper::device_vector<double> workspace(lwork);
        thrust_wrapper::device_vector<int> info(1);
        stat = cusolverDnDpotrf(handle, CUBLAS_FILL_MODE_UPPER, N, A, N,
                                thrust_wrapper::raw_pointer_cast(workspace.data()),
                                lwork, thrust_wrapper::raw_pointer_cast(info.data()));
        assert(CUSOLVER_STATUS_SUCCESS == stat);
        assert(info[0] == 0);

        stat = cusolverDnDpotrs(handle, CUBLAS_FILL_MODE_UPPER, N, N, A, N, B,
                                N, thrust_wrapper::raw_pointer_cast(info.data()));
        assert(CUSOLVER_STATUS_SUCCESS == stat);
        assert(info[0] == 0);

        cusolverDnDestroy(handle);
    }
};
#else
template <>
struct cusolver<policy::host, double> {
    /** Computes the Cholesky factorization and the inverse of a real symmetric
     *  positive definite matrix, looking only in the top half of the symmetric
     *  matrix.
     *  
     *  \param A buffer on host for the symmetric input matrix,
     *           serves as output for the Cholesky decomposition
     *  \param B buffer on host expects identity matrix,
     *           serves as output for the inverse
     *  \param N size of the matrix
     */
    static void potrf(double *A, double *B, int N) {
        char uplo = 'U';
        int info;
        dpotrf_(&uplo, &N, A, &N, &info);
        assert(info == 0);

        dpotrs_(&uplo, &N, &N, A, &N, B, &N, &info);
        assert(info == 0);
    }
};
#endif

template <typename T>
struct almostEqual {
    T epsilon;
    DEVICE_FUNC bool operator()(T x, T y) const {
        return std::abs(x - y) < epsilon;
    }
};

template <typename T>
struct IdentityGenerator {
    std::size_t lda;

    DEVICE_FUNC T operator()(std::size_t index) {
        std::size_t i, j;

        thrust_wrapper::tie(i, j) = unravel_index(index, lda);
        return T(i == j ? 1.0 : 0.0);
    }
};

} // namespace internal

/// \endcond


/** Matrix datatype with arithemtic operations that can run on a THRUST DEVICE.
 *  Storage is column-major order.
 */
template <typename T, typename Policy = policy::host>
class device_matrix {
    static_assert(policy::is_policy<Policy>::value,
                  "The execution policy must meet the requirements");

public:
    using storage_type = typename Policy::template vector<T>;
    using value_type = typename storage_type::value_type;
    using size_type = typename storage_type::size_type;
    using difference_type = typename storage_type::difference_type;
    using reference = typename storage_type::reference;
    using const_reference = typename storage_type::const_reference;
    using pointer = typename storage_type::pointer;
    using const_pointer = typename storage_type::const_pointer;
    using iterator = typename storage_type::iterator;

private:
    // Storage
    size_type m_rows;
    size_type m_cols;
    storage_type m_data;

public:
    device_matrix() : m_rows(0), m_cols(0), m_data(m_rows * m_cols) {}

    explicit device_matrix(size_type rows, size_type cols)
        : m_rows(rows), m_cols(cols), m_data(m_rows * m_cols) {}

    reference operator()(size_type row, size_type col) noexcept {
        return m_data[row + col * m_rows];
    }

    const_reference operator()(size_type row, size_type col) const noexcept {
        return m_data[row + col * m_rows];
    }

    void fill(value_type const &value) noexcept(
        std::is_nothrow_copy_assignable<value_type>::value) {
        thrust_wrapper::fill(Policy::par(), m_data.begin(), m_data.end(), value);
    }

    void
    swap(device_matrix &other) noexcept(noexcept(m_data.swap(other.m_data))) {
        std::swap(m_rows, other.m_rows);
        std::swap(m_cols, other.m_cols);
        m_data.swap(other.m_data);
    }

    pointer data() noexcept { return m_data.data(); }
    const_pointer data() const noexcept { return m_data.data(); }
    size_type rows() const noexcept { return m_rows; }
    size_type cols() const noexcept { return m_cols; }
    size_type size() const noexcept { return m_rows * m_cols; }

    /// \defgroup arithmetic Arithmetic operators
    /// \{

    /// Matrix-matrix multiplication
    device_matrix operator*(device_matrix const &B) const {
        static_assert(std::is_arithmetic<T>::value,
                      "Data type of device_matrix must be arithmetic for "
                      "arithmetic operations");
        assert(m_cols == B.m_rows);
        device_matrix C(m_rows, B.m_cols);
        internal::cublas<Policy, value_type>::gemm(
            thrust_wrapper::raw_pointer_cast(data()),
            thrust_wrapper::raw_pointer_cast(B.data()),
            thrust_wrapper::raw_pointer_cast(C.data()), m_rows, m_cols, B.m_cols);
        return C;
    }

    /// Matrix-vector multiplication
    storage_type operator*(storage_type const &x) const {
        static_assert(std::is_arithmetic<T>::value,
                      "Data type of device_matrix must be arithmetic for "
                      "arithmetic operations");
        assert(m_cols == x.size());
        storage_type y(m_rows);
        internal::cublas<Policy, value_type>::gemv(
            thrust_wrapper::raw_pointer_cast(data()),
            thrust_wrapper::raw_pointer_cast(x.data()),
            thrust_wrapper::raw_pointer_cast(y.data()), m_rows, m_cols);
        return y;
    }

    /// Add element-wise
    device_matrix operator+(device_matrix const &B) const {
        static_assert(std::is_arithmetic<T>::value,
                      "Data type of device_matrix must be arithmetic for "
                      "arithmetic operations");
        assert(m_rows == B.m_rows);
        assert(m_cols == B.m_cols);
        device_matrix C(m_rows, m_cols);
        thrust_wrapper::transform(Policy::par(), data(), data() + size(), 
                                  B.data(), C.data(), 
                                  thrust_wrapper::plus<value_type>{});
        return C;
    }

    /// Subtract element-wise
    device_matrix operator-(device_matrix const &B) const {
        static_assert(std::is_arithmetic<T>::value,
                      "Data type of device_matrix must be arithmetic for "
                      "arithmetic operations");
        assert(m_rows == B.m_rows);
        assert(m_cols == B.m_cols);
        device_matrix C(m_rows, m_cols);
        thrust_wrapper::transform(Policy::par(), data(), data() + size(), 
                                  B.data(), C.data(), 
                                  thrust_wrapper::minus<value_type>{});
        return C;
    }

    /// Negate element-wise
    device_matrix &operator-() {
        static_assert(std::is_arithmetic<T>::value,
                      "Data type of device_matrix must be arithmetic for "
                      "arithmetic operations");
        thrust_wrapper::transform(Policy::par(), data(), data() + size(), 
                                  data(), 
                                  thrust_wrapper::negate<value_type>{});
        return *this;
    }

    /// \}

    /// \defgroup compare Comparison
    /// \{

    /// Compare for exact bit-wise equality
    bool operator==(device_matrix const &B) const {
        return m_rows == B.m_rows && m_cols == B.m_cols &&
               thrust_wrapper::equal(Policy::par(), data(), 
                                     data() + size(), B.data());
    }

    /// Compare for equality with threshold \p epsilon.
    ///
    /// \param epsilon maximum allowed difference
    bool almostEqual(device_matrix const &B,
                     value_type const epsilon =
                         std::numeric_limits<value_type>::epsilon()) const {
        static_assert(std::is_arithmetic<T>::value,
                      "Data type of device_matrix must be arithmetic for "
                      "arithmetic operations");
        return m_rows == B.m_rows && m_cols == B.m_cols &&
               thrust_wrapper::equal(Policy::par(), data(), data() + size(), 
                      B.data(), internal::almostEqual<value_type>{epsilon});
    }

    /// \}

    /// \defgroup linalg Linear algebra
    /// \{

    /// Compute the transpose.
    device_matrix transpose() const {
        static_assert(std::is_same<T, double>::value,
                      "Data type of device_matrix must be floating point for "
                      "BLAS/LAPACK operations");
        device_matrix C(m_cols, m_rows);
        internal::cublas<Policy, value_type>::geam(
            thrust_wrapper::raw_pointer_cast(data()),
            thrust_wrapper::raw_pointer_cast(C.data()), m_rows, m_cols);
        return C;
    }

    /// Compute the inverse and the Cholesky decomposition.
    thrust_wrapper::tuple<device_matrix, device_matrix> inverse_and_cholesky() const {
        static_assert(std::is_same<T, double>::value,
                      "Data type of device_matrix must be floating point for "
                      "BLAS/LAPACK operations");
        assert(m_rows == m_cols);
        device_matrix A = *this;
        device_matrix B = device_matrix::Identity(m_rows, m_cols);

        internal::cusolver<Policy, double>::potrf(
            thrust_wrapper::raw_pointer_cast(A.data()),
            thrust_wrapper::raw_pointer_cast(B.data()), m_rows);

        return thrust_wrapper::make_tuple(B, A);
    }

    /// Compute the inverse.
    device_matrix inverse() const {
        return thrust_wrapper::get<0>(inverse_and_cholesky());
    }

    /// \}

    /// \defgroup generators Convenience generator functions
    /// \{

    /// Generate an identity matrix with size \p rows x \p cols.
    static device_matrix Identity(size_type rows, size_type cols) {
        static_assert(std::is_same<T, double>::value,
                      "Data type of device_matrix must be floating point for "
                      "BLAS/LAPACK operations");
        device_matrix I(rows, cols);
        thrust_wrapper::tabulate(Policy::par(), I.data(), I.data() + I.size(),
                                 internal::IdentityGenerator<value_type>{rows});
        return I;
    }

    /// \}
};


/** Read-only reference to another device_matrix object. More accurately, it
 *  is a separate object that has different routines by which the same region
 *  of memory can be viewed.
 *
 *  Interestingly, by clever abuse of the rows and cols parameters and the data
 *  pointer, certain areas (e.g. blocks) of the referenced matrix could be
 *  accessed conveniently through a device_matrix_view object. 
 *  (This is not used in SD though?)
 */
template <typename T, typename Policy>
class device_matrix_view {
public:
    using storage_type = typename Policy::template vector<T>;
    using value_type = T;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using reference = value_type &;
    using const_reference = value_type const &;
    using pointer = T *;
    using const_pointer = T const *;
    using iterator = pointer;

private:
    // Storage
    size_type m_rows;
    size_type m_cols;
    T *m_data;

public:
    device_matrix_view(pointer data, size_type rows, size_type cols)
        : m_rows(rows), m_cols(cols), m_data(data) {}

    device_matrix_view(device_matrix<T, Policy> &v)
        : m_rows(v.rows()), m_cols(v.cols()),
          m_data(thrust_wrapper::raw_pointer_cast(v.data())) {}

    DEVICE_FUNC reference operator()(size_type row, size_type col) noexcept {
        return m_data[row + col * m_rows];
    }

    DEVICE_FUNC const_reference operator()(size_type row, size_type col) const
        noexcept {
        return m_data[row + col * m_rows];
    }

    DEVICE_FUNC pointer data() noexcept { return m_data; }
    DEVICE_FUNC const_pointer data() const noexcept { return m_data; }
    DEVICE_FUNC size_type rows() const noexcept { return m_rows; }
    DEVICE_FUNC size_type cols() const noexcept { return m_cols; }
    DEVICE_FUNC size_type size() const noexcept { return m_rows * m_cols; }
};

template <typename T, typename Policy = policy::host>
class device_vector_view {
public:
    using storage_type = typename Policy::template vector<T>;
    using value_type = T;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using reference = value_type &;
    using const_reference = value_type const &;
    using pointer = T *;
    using const_pointer = T const *;
    using iterator = pointer;

private:
    // Storage
    size_type m_size;
    T *m_data;

public:
    device_vector_view(pointer data, size_type size)
        : m_size(size), m_data(data) {}

    device_vector_view(storage_type &v)
        : m_size(v.size()), m_data(thrust_wrapper::raw_pointer_cast(v.data())) {}

    DEVICE_FUNC reference operator()(size_type i) noexcept { return m_data[i]; }

    DEVICE_FUNC const_reference operator()(size_type i) const noexcept {
        return m_data[i];
    }

    DEVICE_FUNC pointer data() noexcept { return m_data; }
    DEVICE_FUNC const_pointer data() const noexcept { return m_data; }
    DEVICE_FUNC size_type size() const noexcept { return m_size; }
};


/*
#ifdef SD_THRUST
// These cannot be moved to thrust_wrapper, because device_matrix appears
// Hence, the ifdef
// Just for debugging anyway?
template <typename T>
std::ostream &operator<<(std::ostream &os, 
                         thrust_wrapper::device_vector<T> const &v) {
    os << '{';
    for (std::size_t i = 0; i < v.size(); ++i) {
        os << std::setw(12) << v[i];
        if (i < v.size() - 1) {
            os << ',';
        }
    }
    os << '}';
    return os;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, 
                         thrust_wrapper::host_vector<T> const &v) {
    os << '{';
    for (std::size_t i = 0; i < v.size(); ++i) {
        os << std::setw(12) << v[i];
        if (i < v.size() - 1) {
            os << ',';
        }
    }
    os << '}';
    return os;
}

template <typename T, typename Policy>
std::ostream &operator<<(std::ostream &os, device_matrix<T, Policy> const &m) {
    os << '{';
    for (std::size_t i = 0; i < m.rows(); ++i) {
        if (i > 0) {
            os << ' ';
        }
        os << '{';
        for (std::size_t j = 0; j < m.cols(); ++j) {
            os << std::setw(12) << m(i, j);
            if (j < m.cols() - 1) {
                os << ',';
            }
        }
        os << '}';
        if (i < m.rows() - 1) {
            os << ",\n";
        }
    }
    os << '}';
    return os;
}
#endif
*/
#undef DEVICE_FUNC
#undef MAYBE_UNUSED
#endif

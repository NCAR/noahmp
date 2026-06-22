#ifndef NOAHMPARRAY_H
#define NOAHMPARRAY_H

#include <cstddef>
#include <initializer_list>
#include <array>
#include <NoahmpFatal.H>

#ifndef NDEBUG
// Debug-only bounds checking needs snprintf for the diagnostic; release
// (-DNDEBUG) builds compile the check out entirely.
#include <cstdio>
#endif

// Bounds check for the NoahmpArray access operators. In debug builds (NDEBUG not
// defined) an out-of-range index prints the offending index/range and calls
// noahmp::fatal() -- which routes through the host's abort handler so the failure
// reaches all ranks. In release builds this is a no-op (zero cost). `idx` is a
// long so a signed index below a (possibly negative) lower bound is still caught.
inline void noahmparray_check([[maybe_unused]] long idx, [[maybe_unused]] int lo,
                            [[maybe_unused]] int hi, [[maybe_unused]] int dim) {
#ifndef NDEBUG
    if (idx < lo || idx > hi) {
        char buf[160];
        std::snprintf(buf, sizeof(buf),
                      "NoahmpArray: index %ld out of bounds [%d, %d] on dimension %d",
                      idx, lo, hi, dim);
        noahmp::fatal(buf);
    }
#endif
}

// Floating-point precision of the ERF<->Noah-MP coupling boundary. ERF (the
// client) owns this: the build defines DOUBLE_PREC iff ERF/AMReX is double
// precision, so noahmp_real always matches amrex::Real. This mirrors the
// Fortran c_kind_noahmp parameter in utility/Machine.F90 (C_DOUBLE/C_FLOAT),
// and the two must agree for the C<->Fortran ABI to be correct.
#ifdef DOUBLE_PREC
using noahmp_real = double;
#else
using noahmp_real = float;
#endif

template <typename T>
struct NoahmpArray1D {

    int begin=1;
    int end=0;
    T* data = nullptr;

    constexpr NoahmpArray1D () noexcept = default;

    constexpr NoahmpArray1D (NoahmpArray1D<T> const& rhs) noexcept
        : begin(rhs.begin), end(rhs.end), data(rhs.data) {}

    constexpr NoahmpArray1D(T* fortranArray, int begin, int end) noexcept
        : begin(begin), end(end), data(fortranArray) {}

    // Access operator
    T& operator()(int i) {
        noahmparray_check(i, begin, end, 0);
        return data[i-begin];
    }

    // Const access operator
    const T& operator()(int i) const {
        noahmparray_check(i, begin, end, 0);
        return data[i-begin];
    }

};


template <typename T>
struct NoahmpArray2D {

    std::array<int,2> begin={1,1};
    std::array<int,2> end={0,0};
    T* data = nullptr;

    constexpr NoahmpArray2D () noexcept = default;

    constexpr NoahmpArray2D (NoahmpArray2D<T> const& rhs) noexcept
        : begin(rhs.begin), end(rhs.end), data(rhs.data) {}

    constexpr NoahmpArray2D(T* fortranArray, std::array<int,2> begin, std::array<int,2> end) noexcept
        : begin(begin), end(end), data(fortranArray) {}

    T& operator()(int i, int j) {
        noahmparray_check(i, begin[0], end[0], 0);
        noahmparray_check(j, begin[1], end[1], 1);
        const std::ptrdiff_t n0 = end[0]-begin[0]+1;
        return data[(i-begin[0]) + n0*(j-begin[1])];
    }

    const T& operator()(int i, int j) const {
        noahmparray_check(i, begin[0], end[0], 0);
        noahmparray_check(j, begin[1], end[1], 1);
        const std::ptrdiff_t n0 = end[0]-begin[0]+1;
        return data[(i-begin[0]) + n0*(j-begin[1])];
    }

};


template <typename T>
struct NoahmpArray3D {
    
    std::array<int,3> begin = {1,1,1};
    std::array<int,3> end = {0,0,0};
    T* data = nullptr;

    constexpr NoahmpArray3D() noexcept = default;

    constexpr NoahmpArray3D(NoahmpArray3D<T> const& rhs) noexcept
        : begin(rhs.begin), end(rhs.end), data(rhs.data) {}

    constexpr NoahmpArray3D(T* fortranArray, std::array<int,3> begin, std::array<int,3> end) noexcept
        : begin(begin), end(end), data(fortranArray) {}

    T& operator()(int i, int j, int k) {
        noahmparray_check(i, begin[0], end[0], 0);
        noahmparray_check(j, begin[1], end[1], 1);
        noahmparray_check(k, begin[2], end[2], 2);
        const std::ptrdiff_t n0 = end[0]-begin[0]+1;
        const std::ptrdiff_t n1 = end[1]-begin[1]+1;
        return data[(i-begin[0]) + n0*((j-begin[1]) + n1*(k-begin[2]))];
    }

    const T& operator()(int i, int j, int k) const {
        noahmparray_check(i, begin[0], end[0], 0);
        noahmparray_check(j, begin[1], end[1], 1);
        noahmparray_check(k, begin[2], end[2], 2);
        const std::ptrdiff_t n0 = end[0]-begin[0]+1;
        const std::ptrdiff_t n1 = end[1]-begin[1]+1;
        return data[(i-begin[0]) + n0*((j-begin[1]) + n1*(k-begin[2]))];
    }
};

#endif

#ifndef NOAHARRAY_H
#define NOAHARRAY_H

#include <cstddef>
#include <initializer_list>
#include <array>
 
template <typename T>
struct NoahArray1D {

    int begin=1;
    int end=0;
    T* data;

    constexpr NoahArray1D () noexcept = default;

    constexpr NoahArray1D (NoahArray1D<T> const& rhs) noexcept
        : data(rhs.data), begin(rhs.begin), end(rhs.end) {}

    constexpr NoahArray1D(T* fortranArray, int begin, int end) noexcept
        : data(fortranArray), begin(begin), end(end) {}

    // Access operator
    T& operator()(int i) {
        return data[i-begin];
    }

    // Const access operator
    const T& operator()(int i) const {
        return data[i-begin];
    }

};


template <typename T>
struct NoahArray2D {

    std::array<int,2> begin={1,1};
    std::array<int,2> end={0,0};
    T* data;

    constexpr NoahArray2D () noexcept = default;

    constexpr NoahArray2D (NoahArray2D<T> const& rhs) noexcept
        : data(rhs.data), begin(rhs.begin), end(rhs.end) {}

    constexpr NoahArray2D(T* fortranArray, std::array<int,2> begin, std::array<int,2> end) noexcept
        : data(fortranArray), begin(begin), end(end) {}

    T& operator()(size_t i, size_t j) {
        return data[i-begin[0] + (end[0]-begin[0]+1)*(j-begin[1])];
    }

    const T& operator()(size_t i, size_t j) const {
        return data[i-begin[0] + (end[0]-begin[0]+1)*(j-begin[1])];
    }

};

#endif

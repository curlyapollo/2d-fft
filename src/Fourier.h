#pragma once

#include <complex>

#include "Util.h"

// Abstract class for Fourier transforms
template <typename ComplT = std::complex<double>>
class Fourier {
public:
    // Forward Fourier transform of 1d complex data written in a C-style array
    virtual void Forward1d(int n, const ComplT* in, ComplT* out, unsigned options = 0);

    // Inverse Fourier transform of 1d complex data written in a C-style array
    virtual void Inverse1d(int n, const ComplT* in, ComplT* out, unsigned options = 0);

    // Forward Fourier transform of 2d complex data written in a C-style array
    virtual void Forward2d(int rows, int cols, const ComplT* in, ComplT* out, unsigned options = 0);

    // Inverse Fourier transform of 2d complex data written in a C-style array
    virtual void Inverse2d(int rows, int cols, const ComplT* in, ComplT* out, unsigned options = 0);

    // general Fourier transform, where sign refers to sign in complex exponent i.e. negative for
    // forward transforms and positive for inverse. Doesn't do any normalizing
    virtual void Transform1d(int n, int sign, const ComplT* in, ComplT* out, unsigned options) = 0;

    // same as 1d general Fourier transform but 2d
    virtual void Transform2d(int rows, int cols, int sign, const ComplT* in, ComplT* out,
                             unsigned options);
};

template <typename ComplT>
void Fourier<ComplT>::Forward1d(int n, const ComplT* in, ComplT* out, unsigned int options) {
    Transform1d(n, -1, in, out, options);
}

template <typename ComplT>
void Fourier<ComplT>::Inverse1d(int n, const ComplT* in, ComplT* out, unsigned int options) {
    Transform1d(n, 1, in, out, options);
    for (int i = 0; i < n; ++i) {
        out[i] /= n;
    }
}

template <typename ComplT>
void Fourier<ComplT>::Forward2d(int rows, int cols, const ComplT* in, ComplT* out,
                                unsigned int options) {
    Transform2d(rows, cols, -1, in, out, options);
}

template <typename ComplT>
void Fourier<ComplT>::Inverse2d(int rows, int cols, const ComplT* in, ComplT* out,
                                unsigned int options) {
    Transform2d(rows, cols, 1, in, out, options);
    for (int i = 0; i < rows * cols; ++i) {
        out[i] /= rows * cols;
    }
}

template <typename ComplT>
void Fourier<ComplT>::Transform2d(int rows, int cols, int sign, const ComplT* in, ComplT* out,
                                  unsigned options) {
    // 2d transform computed using only 1d transforms
    // 2d Fourier transform written in linear algebra form would look like F_n X F_m^T
    // So what we do can be rewritten as F_n X F_m^T = (F_n X) F_m^T = (F_m (F_n X)^T)^T

    // in -> temp1
    // untouched -> transformed rows
    ComplT* temp1 = new ComplT[rows * cols];
    for (int row = 0; row < rows; ++row) {
        Transform1d(cols, sign, in + row * cols, temp1 + row * cols, options);
    }

    // temp1 -> temp2
    // transformed rows -> transformed rows, transposed
    ComplT* temp2 = new ComplT[rows * cols];
    Transpose(rows, cols, temp1, temp2);

    // temp2 -> temp1
    // transformed rows, transposed -> transformed rows, transposed, transformed cols
    for (int col = 0; col < cols; ++col) {
        Transform1d(rows, sign, temp2 + col * rows, temp1 + col * rows, options);
    }

    // temp1 -> out
    // transformed rows, transposed, transformed cols -> fully transformed
    Transpose(cols, rows, temp1, out);

    delete[] temp1;
    delete[] temp2;
}

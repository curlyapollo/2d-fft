#pragma once

#include <math.h>

#include <complex>

#include "BitReverse.h"
#include "Fourier.h"
#include "FourierBasic.h"
#include "FourierFast.h"

// class for Fourier transforms where all transforms are done using Cooley-Tukey FFT algo
// if size of an array is a power of 2, algorithm avoids both memory allocation and recursion
// using technic known as bit-reversal permutation
// 2d transforms are done using several applications of 1d transforms
template <typename ComplT = std::complex<double>>
class FourierFastBitReverse : public Fourier<ComplT> {
    template <typename>
    friend class FourierAdvanced;

public:
    void Transform1d(int n, int sign, const ComplT* in, ComplT* out, unsigned options = 0) override;
    virtual ~FourierFastBitReverse(){};

private:
    static void StaticTransform1d(int n, int sign, const ComplT* in, ComplT* out,
                                  unsigned options = 0);
};

template <typename ComplT>
void FourierFastBitReverse<ComplT>::Transform1d(int n, int sign, const ComplT* in, ComplT* out,
                                                unsigned int options) {
    FourierFastBitReverse<ComplT>::StaticTransform1d(n, sign, in, out);
}

template <typename ComplT>
void FourierFastBitReverse<ComplT>::StaticTransform1d(int n, int sign, const ComplT* in,
                                                      ComplT* out, unsigned options) {

    if (!IsPowerOf2(n)) {
        FourierFast<ComplT>::StaticTransform1d(n, sign, in, out);
        return;
    }
    Copy(n, in, out);
    BitReverse* rev = BitReverse::GetInstance(BitLog(n));
    for (int i = 0; i < n; ++i) {
        if (i < rev->Get(i)) {
            std::swap(out[i], out[rev->Get(i)]);
        }
    }
    for (int len = 2; len <= n; len *= 2) {  // from last layer of virtual recursion to first
        double base_power = sign * 2 * M_PI / len;
        for (int i = 0; i < n; i += len) {  // corresponding to different calls on single layer
            for (int j = 0; j < len / 2; ++j) {
                ComplT& c0 = out[i + j];
                ComplT& c1 = out[i + j + len / 2];

                c1 *= std::polar(1., base_power * j);

                ButterflyTransform(c0, c1);
            }
        }
    }
}
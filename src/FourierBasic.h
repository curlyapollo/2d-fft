#pragma once

#include <math.h>

#include <complex>

#include "Fourier.h"

template <typename ComplT = std::complex<double>>
class FourierBasic : public Fourier<ComplT> {
    template <typename>
    friend class FourierFast;

public:
    void Transform1d(int n, int sign, const ComplT* in, ComplT* out, unsigned options = 0) override;
    virtual ~FourierBasic(){};

private:
    static void StaticTransform1d(int n, int sign, const ComplT* in, ComplT* out,
                                  unsigned options = 0);
};

template <typename ComplT>
void FourierBasic<ComplT>::Transform1d(int n, int sign, const ComplT* in, ComplT* out,
                                       unsigned int options) {
    FourierBasic<ComplT>::StaticTransform1d(n, sign, in, out, options);
}

template <typename ComplT>
void FourierBasic<ComplT>::StaticTransform1d(int n, int sign, const ComplT* in, ComplT* out,
                                             unsigned int options) {
    for (int i = 0; i < n; ++i) {
        double base_power = sign * 2 * M_PI * i / n;
        out[i] = 0;
        for (int j = 0; j < n; ++j) {
            out[i] += in[j] * std::polar(1., base_power * j);
        }
    }
}

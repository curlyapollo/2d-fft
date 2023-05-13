#pragma once

#include <math.h>

#include <complex>

#include "Fourier.h"

template <typename ComplT = std::complex<double>>
class FourierBasic : public Fourier<ComplT> {
public:
    void Transform1d(int n, int sign, const ComplT* in, ComplT* out, unsigned options = 0) override;
};

namespace fourier_transforms {

template <typename ComplT = std::complex<double>>
void Transform1dBasic(int n, int sign, const ComplT* in, ComplT* out) {
    for (int i = 0; i < n; ++i) {
        double base_power = sign * 2 * M_PI * i / n;
        out[i] = 0;
        for (int j = 0; j < n; ++j) {
            out[i] += in[j] * std::polar(1., base_power * j);
        }
    }
}

}  // namespace fourier_transforms

template <typename ComplT>
void FourierBasic<ComplT>::Transform1d(int n, int sign, const ComplT* in, ComplT* out,
                                       unsigned int options) {
    fourier_transforms::Transform1dBasic(n, sign, in, out);
}

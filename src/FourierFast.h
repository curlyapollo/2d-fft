#pragma once

#include <math.h>

#include <complex>

#include "Fourier.h"
#include "FourierBasic.h"

// class for Fourier transforms where all transforms are done using classic Cooley-Tukey FFT algo.
// 2d transforms are done using several applications of 1d transforms
template <typename ComplT = std::complex<double>>
class FourierFast : public Fourier<ComplT> {
private:
    void Transform1d(int n, int sign, const ComplT* in, ComplT* out, unsigned options = 0) override;
};

namespace fourier_transforms {

template <typename ComplT = std::complex<double>>
void Transform1dFast(int n, int sign, const ComplT* in, ComplT* out) {
    if (n % 2 != 0) {
        fourier_transforms::Transform1dBasic(n, sign, in, out);
        return;
    }
    ComplT* even_odd = new ComplT[n];

    for (int i = 0; 2 * i < n; ++i) {
        even_odd[i] = in[2 * i];
        even_odd[i + n / 2] = in[2 * i + 1];
    }
    fourier_transforms::Transform1dFast(n / 2, sign, even_odd, out);
    fourier_transforms::Transform1dFast(n / 2, sign, even_odd + n / 2, out + n / 2);

    double base_power = sign * 2 * M_PI / n;
    for (int i = 0; 2 * i < n; ++i) {
        ComplT temp = out[i] + out[i + n / 2] * std::polar(1., base_power * i);
        out[i + n / 2] = out[i] - out[i + n / 2] * std::polar(1., base_power * i);
        out[i] = temp;
    }

    delete[] even_odd;
}

}  // namespace fourier_transforms

template <typename ComplT>
void FourierFast<ComplT>::Transform1d(int n, int sign, const ComplT* in, ComplT* out,
                                      unsigned int options) {
    fourier_transforms::Transform1dFast(n, sign, in, out);
}

#pragma once

#include <math.h>

#include <complex>

#include "Fourier.h"
#include "FourierBasic.h"
#include "Util.h"

// class for Fourier transforms where all transforms are done using classic Cooley-Tukey FFT algo.
// 2d transforms are done using several applications of 1d transforms
template <typename ComplT = std::complex<double>>
class FourierFast : public Fourier<ComplT> {
    template <typename>
    friend class FourierAdvanced;

    template <typename>
    friend class FourierFastBitReverse;

public:
    void Transform1d(int n, int sign, const ComplT* in, ComplT* out, unsigned options = 0) override;
    virtual ~FourierFast(){};

private:
    static void StaticTransform1d(int n, int sign, const ComplT* in, ComplT* out,
                                  unsigned options = 0);
};

template <typename ComplT>
void FourierFast<ComplT>::Transform1d(int n, int sign, const ComplT* in, ComplT* out,
                                      unsigned int options) {
    FourierFast<ComplT>::StaticTransform1d(n, sign, in, out);
}

template <typename ComplT>
void FourierFast<ComplT>::StaticTransform1d(int n, int sign, const ComplT* in, ComplT* out,
                                            unsigned options) {

    if (n % 2 != 0) {
        FourierBasic<ComplT>::StaticTransform1d(n, sign, in, out);
        return;
    }
    ComplT* even_odd = new ComplT[n];

    for (int i = 0; 2 * i < n; ++i) {
        even_odd[i] = in[2 * i];
        even_odd[i + n / 2] = in[2 * i + 1];
    }
    FourierFast<ComplT>::StaticTransform1d(n / 2, sign, even_odd, out);
    FourierFast<ComplT>::StaticTransform1d(n / 2, sign, even_odd + n / 2, out + n / 2);

    double base_power = sign * 2 * M_PI / n;
    for (int i = 0; 2 * i < n; ++i) {
        ComplT& c0 = out[i];
        ComplT& c1 = out[i + n / 2];

        c1 *= std::polar(1., base_power * i);

        ButterflyTransform(c0, c1);
    }

    delete[] even_odd;
}

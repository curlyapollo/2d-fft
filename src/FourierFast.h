#pragma once

#include <math.h>

#include <complex>

#include "Fourier.h"

// class for Fourier transforms where all transforms are done using classic Cooley-Tukey FFT algo.
// 2d transforms are done using several applications of 1d transforms
template <typename ComplT = std::complex<double>>
class FourierFast : public Fourier<ComplT> {
private:
    void Transform1d(int n, int sign, const ComplT* in, ComplT* out, unsigned options = 0) override;
};

template <typename ComplT>
void FourierFast<ComplT>::Transform1d(int n, int sign, const ComplT* in, ComplT* out,
                                      unsigned int options) {
    // TODO
}

#pragma once

#include <math.h>

#include <complex>

#include "Fourier.h"
#include "FourierFast.h"

// class for Fourier transforms where all 1d transforms are done using classic Cooley-Tukey FFT
// algo and 2d transforms are done using algos from Tutatchikov, Starovoytov and Noskov 's papers
template <typename ComplT = std::complex<double>>
class FourierAdvanced : public FourierFast<ComplT> {
private:
    void Transform2d(int rows, int cols, int sign, const ComplT* in, ComplT* out,
                     unsigned options = 0) override;
};

template <typename ComplT>
void FourierAdvanced<ComplT>::Transform2d(int rows, int cols, int sign, const ComplT* in,
                                          ComplT* out, unsigned options) {
    // TODO
}

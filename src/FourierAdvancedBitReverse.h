#pragma once

#include <math.h>

#include <complex>

#include "Fourier.h"
#include "FourierAdvanced.h"
#include "FourierFastBitReverse.h"

// class for Fourier transforms where all 1d transforms are done using classic Cooley-Tukey FFT
// algo and 2d transforms are done using algos from Tutatchikov, Starovoytov and Noskov 's papers
template <typename ComplT = std::complex<double>>
class FourierAdvancedBitReverse : public FourierFastBitReverse<ComplT> {
public:
    void Transform2d(int rows, int cols, int sign, const ComplT* in, ComplT* out,
                     unsigned options = 0) override;
    virtual ~FourierAdvancedBitReverse(){};

private:
    static void StaticTransform2d(int rows, int cols, int sign, const ComplT* in, ComplT* out,
                                  unsigned options = 0);
};

template <typename ComplT>
void FourierAdvancedBitReverse<ComplT>::Transform2d(int rows, int cols, int sign, const ComplT* in,
                                                    ComplT* out, unsigned options) {
    FourierAdvancedBitReverse<ComplT>::StaticTransform2d(rows, cols, sign, in, out);
}

template <typename ComplT>
void FourierAdvancedBitReverse<ComplT>::StaticTransform2d(int rows, int cols, int sign,
                                                          const ComplT* in, ComplT* out,
                                                          unsigned int options) {
    if (!IsPowerOf2(rows) || !IsPowerOf2(cols)) {
        FourierAdvanced<ComplT>::StaticTransform2d(rows, cols, sign, in, out);
        return;
    }
    Copy(rows * cols, in, out);
    BitReverse* rev_rows = BitReverse::GetInstance(BitLog(rows));
    BitReverse* rev_cols = BitReverse::GetInstance(BitLog(cols));
    for (int i = 0; i < rows; ++i) {
        int rev_i = rev_rows->Get(i);
        if (i < rev_i) {
            for (int j = 0; j < cols; ++j) {
                std::swap(out[IND_2D(rows, cols, i, j)],
                          out[IND_2D(rows, cols, rev_i, rev_cols->Get(j))]);
            }
        } else if (i == rev_i) {
            for (int j = 0; j < cols; ++j) {
                if (j < rev_cols->Get(j)) {
                    std::swap(out[IND_2D(rows, cols, i, j)],
                              out[IND_2D(rows, cols, rev_i, rev_cols->Get(j))]);
                }
            }
        }
    }
    for (int len = 2; len <= std::min(rows, cols);
         len *= 2) {  // from last layer of virtual recursion to first
        double base_power = sign * 2 * M_PI / len;
        for (int i = 0; i < rows; i += len) {  // corresponding to different calls on single layer
            for (int j = 0; j < cols; j += len) {
                for (int u = 0; u < len / 2; ++u) {
                    for (int v = 0; v < len / 2; ++v) {
                        ComplT& c00 = out[IND_2D(rows, cols, i + u, j + v)];
                        ComplT& c01 = out[IND_2D(rows, cols, i + u, j + v + len / 2)];
                        ComplT& c10 = out[IND_2D(rows, cols, i + u + len / 2, j + v)];
                        ComplT& c11 = out[IND_2D(rows, cols, i + u + len / 2, j + v + len / 2)];

                        ComplT temp00 = c00 + c01 * std::polar(1., base_power * v) +
                                        c10 * std::polar(1., base_power * u) +
                                        c11 * std::polar(1., base_power * (v + u));
                        ComplT temp01 = c00 - c01 * std::polar(1., base_power * v) +
                                        c10 * std::polar(1., base_power * u) -
                                        c11 * std::polar(1., base_power * (v + u));
                        ComplT temp10 = c00 + c01 * std::polar(1., base_power * v) -
                                        c10 * std::polar(1., base_power * u) -
                                        c11 * std::polar(1., base_power * (v + u));
                        ComplT temp11 = c00 - c01 * std::polar(1., base_power * v) -
                                        c10 * std::polar(1., base_power * u) +
                                        c11 * std::polar(1., base_power * (v + u));
                        c00 = temp00;
                        c01 = temp01;
                        c10 = temp10;
                        c11 = temp11;
                    }
                }
            }
        }
    }
    if (cols > rows) {
        for (int len = rows * 2; len <= cols; len *= 2) {
            double base_power = sign * 2 * M_PI / len;
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j += len) {
                    for (int v = 0; v < len / 2; ++v) {
                        ComplT& c0 = out[IND_2D(rows, cols, i, j + v)];
                        ComplT& c1 = out[IND_2D(rows, cols, i, j + v + len / 2)];
                        ComplT temp = c0 + c1 * std::polar(1., base_power * v);
                        c1 = c0 - c1 * std::polar(1., base_power * v);
                        c0 = temp;
                    }
                }
            }
        }
    } else if (rows > cols) {
        for (int len = cols * 2; len <= rows; len *= 2) {
            double base_power = sign * 2 * M_PI / len;
            for (int i = 0; i < rows; i += len) {
                for (int u = 0; u < len / 2; ++u) {
                    for (int j = 0; j < cols; j++) {
                        ComplT& c0 = out[IND_2D(rows, cols, i + u, j)];
                        ComplT& c1 = out[IND_2D(rows, cols, i + u + len / 2, j)];
                        ComplT temp = c0 + c1 * std::polar(1., base_power * u);
                        c1 = c0 - c1 * std::polar(1., base_power * u);
                        c0 = temp;
                    }
                }
            }
        }
    }
}

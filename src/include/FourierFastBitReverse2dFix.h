#pragma once

#include <math.h>

#include <complex>

#include "BitReverse.h"
#include "Fourier.h"
#include "FourierBasic.h"
#include "FourierFastBitReverse.h"

// TODO
template <typename ComplT = std::complex<double>>
class FourierFastBitReverse2dFix : public FourierFastBitReverse<ComplT> {
    template <typename>
    friend class FourierAdvanced;

public:
    void Transform2d(int rows, int cols, int sign, const ComplT* in, ComplT* out,
                     unsigned options = 0) override;
    virtual ~FourierFastBitReverse2dFix(){};

private:
    static void StaticTransform2d(int rows, int cols, int sign, const ComplT* in, ComplT* out,
                                  unsigned options = 0);
};

template <typename ComplT>
void FourierFastBitReverse2dFix<ComplT>::Transform2d(int rows, int cols, int sign, const ComplT* in,
                                                     ComplT* out, unsigned options) {
    FourierFastBitReverse2dFix<ComplT>::StaticTransform2d(rows, cols, sign, in, out);
}

template <typename ComplT>
void FourierFastBitReverse2dFix<ComplT>::StaticTransform2d(int rows, int cols, int sign,
                                                           const ComplT* in, ComplT* out,
                                                           unsigned options) {

    if (!IsPowerOf2(rows) || !IsPowerOf2(cols)) {
        FourierFastBitReverse<ComplT>().Transform2d(rows, cols, sign, in, out);
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
    for (int len = 2; len <= cols; len *= 2) {
        double base_power = sign * 2 * M_PI / len;
        for (int i = 0; i < rows; ++i) {
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

    for (int len = 2; len <= rows; len *= 2) {
        double base_power = sign * 2 * M_PI / len;
        for (int i = 0; i < rows; i += len) {
            for (int u = 0; u < len / 2; ++u) {
                for (int j = 0; j < cols; ++j) {
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
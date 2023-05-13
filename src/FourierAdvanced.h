#pragma once

#include <math.h>

#include <complex>

#include "Fourier.h"
#include "FourierFast.h"

// class for Fourier transforms where all 1d transforms are done using classic Cooley-Tukey FFT
// algo and 2d transforms are done using algos from Tutatchikov, Starovoytov and Noskov 's papers
template <typename ComplT = std::complex<double>>
class FourierAdvanced : public FourierFast<ComplT> {
public:
    void Transform2d(int rows, int cols, int sign, const ComplT* in, ComplT* out,
                     unsigned options = 0) override;
    virtual ~FourierAdvanced(){};
};

namespace fourier_transforms {

template <typename ComplT = std::complex<double>>
void Transform2dFast(int rows, int cols, int sign, const ComplT* in, ComplT* out) {
    if (rows < cols) {
        ComplT* tmp = new ComplT[rows * cols];
        Transpose(rows, cols, in, out);
        Transform2dFast(cols, rows, sign, out, tmp);
        Transpose(cols, rows, tmp, out);
        delete[] tmp;
        return;
    }
    if (cols < rows && cols % 2 == 0) {
        ComplT* tmp = new ComplT[rows * cols];
        for (int i = 0; i * 2 < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                tmp[IND_2D(rows, cols, i, j)] = in[IND_2D(rows, cols, 2 * i, j)];
                tmp[IND_2D(rows, cols, i + rows / 2, j)] = in[IND_2D(rows, cols, 2 * i + 1, j)];
            }
        }
        fourier_transforms::Transform2dFast(rows / 2, cols, sign, tmp, out);
        fourier_transforms::Transform2dFast(rows / 2, cols, sign, tmp + (rows / 2) * cols,
                                            out + (rows / 2) * cols);

        double base_power = sign * 2 * M_PI / rows;
        for (int i = 0; i * 2 < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                ComplT temp =
                    out[IND_2D(rows, cols, i, j)] +
                    out[IND_2D(rows, cols, i + rows / 2, j)] * std::polar(1., base_power * i);
                out[IND_2D(rows, cols, i + rows / 2, j)] =
                    out[IND_2D(rows, cols, i, j)] -
                    out[IND_2D(rows, cols, i + rows / 2, j)] * std::polar(1., base_power * i);
                out[IND_2D(rows, cols, i, j)] = temp;
            }
        }

        delete[] tmp;
        return;
    }
    if (cols == rows && cols % 2 == 0) {
        ComplT* tmp_in = new ComplT[rows * cols];
        ComplT* tmp_out = new ComplT[rows * cols];
        ComplT* g00 = tmp_in;
        ComplT* g01 = tmp_in + rows * cols / 4;
        ComplT* g10 = tmp_in + rows * cols / 2;
        ComplT* g11 = tmp_in + rows * cols * 3 / 4;

        ComplT* res00 = tmp_out;
        ComplT* res01 = tmp_out + rows * cols / 4;
        ComplT* res10 = tmp_out + rows * cols / 2;
        ComplT* res11 = tmp_out + rows * cols * 3 / 4;

        for (int i = 0; i * 2 < rows; ++i) {
            for (int j = 0; j * 2 < rows; ++j) {
                g00[IND_2D(rows / 2, cols / 2, i, j)] = in[IND_2D(row, cols, 2 * i, 2 * j)];
                g01[IND_2D(rows / 2, cols / 2, i, j)] = in[IND_2D(row, cols, 2 * i, 2 * j + 1)];
                g10[IND_2D(rows / 2, cols / 2, i, j)] = in[IND_2D(row, cols, 2 * i + 1, 2 * j)];
                g11[IND_2D(rows / 2, cols / 2, i, j)] = in[IND_2D(row, cols, 2 * i + 1, 2 * j + 1)];
            }
        }

        fourier_transforms::Transform2dFast(rows / 2, cols / 2, sign, g00, res00);
        fourier_transforms::Transform2dFast(rows / 2, cols / 2, sign, g01, res01);
        fourier_transforms::Transform2dFast(rows / 2, cols / 2, sign, g10, res10);
        fourier_transforms::Transform2dFast(rows / 2, cols / 2, sign, g11, res11);
        double base_power = sign * 2 * M_PI / rows;

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                int i1 = (i < rows / 2) ? i : i - rows / 2;
                int j1 = (j < rows / 2) ? j : j - rows / 2;
                out[IND_2D(rows, cols, i, j)] =
                    res00[IND_2D(rows / 2, cols / 2, i1, j1)] +
                    res01[IND_2D(rows / 2, cols / 2, i1, j1)] * std::polar(1., base_power * j) +
                    res10[IND_2D(rows / 2, cols / 2, i1, j1)] * std::polar(1., base_power * i) +
                    res11[IND_2D(rows / 2, cols / 2, i1, j1)] *
                        std::polar(1., base_power * (i + j));
            }
        }
        delete[] tmp_in;
        delete[] tmp_out;
        return;
    } else {
        FourierFast().Transform2d(rows, cols, sign, in, out, 0);
    }
}

}  // namespace fourier_transforms

template <typename ComplT>
void FourierAdvanced<ComplT>::Transform2d(int rows, int cols, int sign, const ComplT* in,
                                          ComplT* out, unsigned options) {
    fourier_transforms::Transform2dFast(rows, cols, sign, in, out);
}

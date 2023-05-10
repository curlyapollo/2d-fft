#pragma once

#define IND_2D(d1, d2, ind1, ind2) (d2) * (ind1) + (ind2)

template <typename T>
void Transpose(int n1, int n2, const T* in, T* out) {
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            out[IND_2D(n2, n1, j, i)] = in[IND_2D(n1, n2, i, j)];
        }
    }
}

template <typename T>
void TransposeInplace(int n, T* in) {
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            swap(in[IND_2D(n, n, i, j)], in[IND_2D(n, n, j, i)]);
        }
    }
}
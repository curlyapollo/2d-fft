#include "BitReverse.h"

static const int MAX_LOG = 30;

BitReverse::BitReverse(int log) {
    data_ = std::make_unique<int[]>(1 << log);
    for (int i = 0; i < (1 << log); ++i) {
        int ans = 0;
        for (int j = 0; j < log; ++j) {
            if ((i >> j) & 1) {
                ans ^= (1 << (log - 1 - j));
            }
        }
        data_[i] = ans;
    }
}

BitReverse* BitReverse::GetInstance(int log) {
    static BitReverse* instances[MAX_LOG + 1] = {0};
    if (instances[log] == nullptr) {
        instances[log] = new BitReverse(log);
    }
    return instances[log];
}

int BitReverse::Get(int num) {
    return data_[num];
}

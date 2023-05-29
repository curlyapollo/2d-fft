#pragma once

#include <memory>

class BitReverse {
    BitReverse(int log);

public:
    static BitReverse* GetInstance(int log);
    int Get(int num);

private:
    std::unique_ptr<int[]> data_;
};
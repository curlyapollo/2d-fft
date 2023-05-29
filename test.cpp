#include <catch2/catch_test_macros.hpp>
#include <vector>

#include "_includes.h"

using namespace std::complex_literals;

bool AllClose(size_t n, const std::complex<double> vec1[], const std::complex<double> vec2[]) {
    for (size_t i = 0; i < n; ++i) {
        if (abs(vec1[i] - vec2[i]) > 1e-6 * std::max(1.0, abs(vec1[i]))) {
            return false;
        }
    }
    return true;
}

TEST_CASE("Precomputed 1d #1") {
    int n = 6;
    std::complex<double> in[] = {2, 2, 8, 3, 3, 2};
    std::complex<double> out[] = {20. + 0.00000000e+00i,  -4.5 - 4.33012702e+00i,
                                  -2.5 + 4.33012702e+00i, 6. + 8.88178420e-16i,
                                  -2.5 - 4.33012702e+00i, -4.5 + 4.33012702e+00i};
    std::complex<double> res_forward[n];
    std::complex<double> res_backward[n];

    std::unique_ptr<Fourier<>> transforms[] = {
        std::make_unique<FourierBasic<>>(),
        std::make_unique<FourierFast<>>(),
        std::make_unique<FourierFastBitReverse<>>(),
        std::make_unique<FourierAdvanced<>>(),
        std::make_unique<FourierAdvancedBitReverse<>>(),
        std::make_unique<FourierFastBitReverse2dFix<>>(),
    };

    for (auto& transform : transforms) {
        transform->Forward1d(n, in, res_forward);
        REQUIRE(AllClose(n, res_forward, out));
        transform->Inverse1d(n, out, res_backward);
        REQUIRE(AllClose(n, res_backward, in));
    }
}

TEST_CASE("Precomputed 1d #2") {
    int n = 6;
    std::complex<double> in[] = {1, 2, 3, 4, 5, 6};
    std::complex<double> out[] = {21. + 0.00000000e+00i, -3. + 5.19615242e+00i,
                                  -3. + 1.73205081e+00i, -3. - 4.44089210e-16i,
                                  -3. - 1.73205081e+00i, -3. - 5.19615242e+00i};
    std::complex<double> res_forward[n];
    std::complex<double> res_backward[n];

    std::unique_ptr<Fourier<>> transforms[] = {
        std::make_unique<FourierBasic<>>(),
        std::make_unique<FourierFast<>>(),
        std::make_unique<FourierFastBitReverse<>>(),
        std::make_unique<FourierAdvanced<>>(),
        std::make_unique<FourierAdvancedBitReverse<>>(),
        std::make_unique<FourierFastBitReverse2dFix<>>(),
    };

    for (auto& transform : transforms) {
        transform->Forward1d(n, in, res_forward);
        REQUIRE(AllClose(n, res_forward, out));
        transform->Inverse1d(n, out, res_backward);
        REQUIRE(AllClose(n, res_backward, in));
    }
}

TEST_CASE("Precomputed 2d #1") {
    int n1 = 2, n2 = 4;
    // clang-format off
    std::complex<double> in[] = {2, 2, 8,         3,
                                 3, 2, 1. + 0.5i, 2};
    std::complex<double> out[] = {23. + 0.5i, -4. + 0.5i, 5. + 0.5i, -4. - 1.5i,
                                  7. - 0.5i,  -8. + 1.5i, 5. - 0.5i, -8. - 0.5i};
    // clang-format on
    std::complex<double> res_forward[n1 * n2];
    std::complex<double> res_backward[n1 * n2];

    std::unique_ptr<Fourier<>> transforms[] = {
        std::make_unique<FourierBasic<>>(),
        std::make_unique<FourierFast<>>(),
        std::make_unique<FourierFastBitReverse<>>(),
        std::make_unique<FourierAdvanced<>>(),
        std::make_unique<FourierAdvancedBitReverse<>>(),
        std::make_unique<FourierFastBitReverse2dFix<>>(),
    };

    for (auto& transform : transforms) {
        transform->Forward2d(n1, n2, in, res_forward);
        REQUIRE(AllClose(n1 * n2, res_forward, out));
        transform->Inverse2d(n1, n2, out, res_backward);
        REQUIRE(AllClose(n1 * n2, res_backward, in));
    }
}

TEST_CASE("Precomputed 2d #2") {
    int n1 = 4, n2 = 4;
    // clang-format off
    std::complex<double> in[] = {2, 2,         8,         3,
                                 3, 2,         1. + 0.5i, 2,
                                 2, 8,         3,         3,
                                 2, 1. + 0.5i, 2,         2};
    std::complex<double> out[] = {46. + 1.i, -4.5 - 3.5i, 0. + 0.i,  -5.5 + 2.5i,
                                  -1. - 1.i, -6.5 + 4.5i, 12. + 1.i, -4.5 - 8.5i,
                                  16. - 1.i, -9.5 - 4.5i, -2. + 0.i, -8.5 + 5.5i,
                                  -1. + 1.i, -3.5 + 7.5i, 10. - 1.i, -5.5 - 3.5i};

    // clang-format on
    std::complex<double> res_forward[n1 * n2];
    std::complex<double> res_backward[n1 * n2];

    std::unique_ptr<Fourier<>> transforms[] = {
        std::make_unique<FourierBasic<>>(),
        std::make_unique<FourierFast<>>(),
        std::make_unique<FourierFastBitReverse<>>(),
        std::make_unique<FourierAdvanced<>>(),
        std::make_unique<FourierAdvancedBitReverse<>>(),
        std::make_unique<FourierFastBitReverse2dFix<>>(),
    };

    for (auto& transform : transforms) {
        transform->Forward2d(n1, n2, in, res_forward);
        REQUIRE(AllClose(n1 * n2, res_forward, out));
        transform->Inverse2d(n1, n2, out, res_backward);
        REQUIRE(AllClose(n1 * n2, res_backward, in));
    }
}

#define TRIES 1000

TEST_CASE("Backward of forward is id, 1d #1") {
    int n = 16;
    std::complex<double> in[n];
    std::complex<double> res_forward[n];
    std::complex<double> res_backward[n];

    std::unique_ptr<Fourier<>> transforms[] = {
        std::make_unique<FourierBasic<>>(),
        std::make_unique<FourierFast<>>(),
        std::make_unique<FourierFastBitReverse<>>(),
        std::make_unique<FourierAdvanced<>>(),
        std::make_unique<FourierAdvancedBitReverse<>>(),
        std::make_unique<FourierFastBitReverse2dFix<>>(),
    };

    for (int i = 0; i < TRIES; ++i) {
        for (int j = 0; j < n; ++j) {
            in[j] = rand();
        }

        for (auto& transform : transforms) {
            transform->Forward1d(n, in, res_forward);
            transform->Inverse1d(n, res_forward, res_backward);
            REQUIRE(AllClose(n, in, res_backward));
        }
    }
}

TEST_CASE("Backward of forward is id, 1d #2") {
    int n = 120;
    std::complex<double> in[n];
    std::complex<double> res_forward[n];
    std::complex<double> res_backward[n];

    std::unique_ptr<Fourier<>> transforms[] = {
        std::make_unique<FourierBasic<>>(),
        std::make_unique<FourierFast<>>(),
        std::make_unique<FourierFastBitReverse<>>(),
        std::make_unique<FourierAdvanced<>>(),
        std::make_unique<FourierAdvancedBitReverse<>>(),
        std::make_unique<FourierFastBitReverse2dFix<>>(),
    };

    for (int i = 0; i < TRIES; ++i) {
        for (int j = 0; j < n; ++j) {
            in[j] = rand();
        }

        for (auto& transform : transforms) {
            transform->Forward1d(n, in, res_forward);
            transform->Inverse1d(n, res_forward, res_backward);
            REQUIRE(AllClose(n, in, res_backward));
        }
    }
}

TEST_CASE("Backward of forward is id, 2d #1") {
    int n1 = 8, n2 = 8;
    std::complex<double> in[n1 * n2];
    std::complex<double> res_forward[n1 * n2];
    std::complex<double> res_backward[n1 * n2];

    std::unique_ptr<Fourier<>> transforms[] = {
        std::make_unique<FourierBasic<>>(),
        std::make_unique<FourierFast<>>(),
        std::make_unique<FourierFastBitReverse<>>(),
        std::make_unique<FourierAdvanced<>>(),
        std::make_unique<FourierAdvancedBitReverse<>>(),
        std::make_unique<FourierFastBitReverse2dFix<>>(),
    };

    for (int i = 0; i < TRIES; ++i) {
        for (int j = 0; j < n1 * n2; ++j) {
            in[j] = rand();
        }

        for (auto& transform : transforms) {
            transform->Forward2d(n1, n2, in, res_forward);
            transform->Inverse2d(n1, n2, res_forward, res_backward);
            REQUIRE(AllClose(n1 * n2, in, res_backward));
        }
    }
}

TEST_CASE("Backward of forward is id, 2d #2") {
    int n1 = 16, n2 = 4;
    std::complex<double> in[n1 * n2];
    std::complex<double> res_forward[n1 * n2];
    std::complex<double> res_backward[n1 * n2];

    std::unique_ptr<Fourier<>> transforms[] = {
        std::make_unique<FourierBasic<>>(),
        std::make_unique<FourierFast<>>(),
        std::make_unique<FourierFastBitReverse<>>(),
        std::make_unique<FourierAdvanced<>>(),
        std::make_unique<FourierAdvancedBitReverse<>>(),
        std::make_unique<FourierFastBitReverse2dFix<>>(),
    };

    for (int i = 0; i < TRIES; ++i) {
        for (int j = 0; j < n1 * n2; ++j) {
            in[j] = rand();
        }

        for (auto& transform : transforms) {
            transform->Forward2d(n1, n2, in, res_forward);
            transform->Inverse2d(n1, n2, res_forward, res_backward);
            REQUIRE(AllClose(n1 * n2, in, res_backward));
        }
    }
}

TEST_CASE("Rand vs Basic 1") {
    int n1 = 16, n2 = 2;
    std::complex<double> in[n1 * n2];
    std::complex<double> out[n1 * n2];
    std::complex<double> res_forward[n1 * n2];
    std::complex<double> res_backward[n1 * n2];

    std::unique_ptr<Fourier<>> transforms[] = {
        std::make_unique<FourierFast<>>(),
        std::make_unique<FourierFastBitReverse<>>(),
        std::make_unique<FourierAdvanced<>>(),
        std::make_unique<FourierAdvancedBitReverse<>>(),
        std::make_unique<FourierFastBitReverse2dFix<>>(),
    };

    for (int i = 0; i < TRIES; ++i) {
        for (int j = 0; j < n1 * n2; ++j) {
            in[j] = rand();
        }
        FourierBasic<>().Forward2d(n1, n2, in, out);

        for (auto& transform : transforms) {
            transform->Forward2d(n1, n2, in, res_forward);
            REQUIRE(AllClose(n1 * n2, res_forward, out));
            transform->Inverse2d(n1, n2, out, res_backward);
            REQUIRE(AllClose(n1 * n2, res_backward, in));
        }
    }
}

TEST_CASE("Rand vs Basic 2") {
    int n1 = 4, n2 = 8;
    std::complex<double> in[n1 * n2];
    std::complex<double> out[n1 * n2];
    std::complex<double> res_forward[n1 * n2];
    std::complex<double> res_backward[n1 * n2];

    std::unique_ptr<Fourier<>> transforms[] = {
        std::make_unique<FourierFast<>>(),
        std::make_unique<FourierFastBitReverse<>>(),
        std::make_unique<FourierAdvanced<>>(),
        std::make_unique<FourierAdvancedBitReverse<>>(),
        std::make_unique<FourierFastBitReverse2dFix<>>(),
    };

    for (int i = 0; i < TRIES; ++i) {
        for (int j = 0; j < n1 * n2; ++j) {
            in[j] = rand();
        }
        FourierBasic<>().Forward2d(n1, n2, in, out);

        for (auto& transform : transforms) {
            transform->Forward2d(n1, n2, in, res_forward);
            REQUIRE(AllClose(n1 * n2, res_forward, out));
            transform->Inverse2d(n1, n2, out, res_backward);
            REQUIRE(AllClose(n1 * n2, res_backward, in));
        }
    }
}

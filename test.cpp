#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <vector>

#include "Fourier.h"
#include "FourierAdvanced.h"
#include "FourierAdvancedBitReverse.h"
#include "FourierBasic.h"
#include "FourierFast.h"
#include "FourierFastBitReverse.h"
#include "FourierFastBitReverse2dFix.h"

using namespace std::complex_literals;

bool is_identity(size_t n, const std::complex<double> vec1[], const std::complex<double> vec2[]) {
    for (size_t i = 0; i < n; ++i) {
        if (abs(vec1[i] - vec2[i]) > 1e-6) {
            return false;
        }
    }
    return true;
}

TEST_CASE("test1") {
    int n = 6;
    std::complex<double> ans[] = {20. + 0.00000000e+00i,  -4.5 - 4.33012702e+00i,
                                  -2.5 + 4.33012702e+00i, 6. + 8.88178420e-16i,
                                  -2.5 - 4.33012702e+00i, -4.5 + 4.33012702e+00i};
    std::complex<double> vec1[] = {2, 2, 8, 3, 3, 2};
    std::complex<double> transformed[n];

    std::unique_ptr<Fourier<>> transforms[] = {
        std::make_unique<FourierBasic<>>(),
        std::make_unique<FourierFast<>>(),
        std::make_unique<FourierFastBitReverse<>>(),
        std::make_unique<FourierAdvanced<>>(),
        std::make_unique<FourierAdvancedBitReverse<>>(),
        std::make_unique<FourierFastBitReverse2dFix<>>(),
    };

    for (auto& transform : transforms) {
        transform->Forward1d(n, vec1, transformed);
        REQUIRE(is_identity(n, transformed, ans));
    }
}

TEST_CASE("test2") {
    int n = 6;
    std::complex<double> ans[] = {21. + 0.00000000e+00i, -3. + 5.19615242e+00i,
                                  -3. + 1.73205081e+00i, -3. - 4.44089210e-16i,
                                  -3. - 1.73205081e+00i, -3. - 5.19615242e+00i};
    std::complex<double> vec1[] = {1, 2, 3, 4, 5, 6};
    std::complex<double> transformed[n];

    std::unique_ptr<Fourier<>> transforms[] = {
        std::make_unique<FourierBasic<>>(),
        std::make_unique<FourierFast<>>(),
        std::make_unique<FourierFastBitReverse<>>(),
        std::make_unique<FourierAdvanced<>>(),
        std::make_unique<FourierAdvancedBitReverse<>>(),
        std::make_unique<FourierFastBitReverse2dFix<>>(),
    };

    for (auto& transform : transforms) {
        transform->Forward1d(n, vec1, transformed);
        REQUIRE(is_identity(n, transformed, ans));
    }
}

TEST_CASE("test3") {
    int n = 6;
    std::complex<double> ans[6];
    std::complex<double> vec1[] = {2, 2, 8, 3, 3, 2};
    std::complex<double> transformed[n];

    std::unique_ptr<Fourier<>> transforms[] = {
        std::make_unique<FourierBasic<>>(),
        std::make_unique<FourierFast<>>(),
        std::make_unique<FourierFastBitReverse<>>(),
        std::make_unique<FourierAdvanced<>>(),
        std::make_unique<FourierAdvancedBitReverse<>>(),
        std::make_unique<FourierFastBitReverse2dFix<>>(),
    };

    for (auto& transform : transforms) {
        transform->Forward1d(n, vec1, transformed);
        transform->Inverse1d(n, transformed, ans);
        REQUIRE(is_identity(n, ans, vec1));
    }
}

TEST_CASE("test4") {
    int n = 4;
    std::complex<double> ans[] = {46. + 1.i, -4.5 - 3.5i, 0. + 0.i,  -5.5 + 2.5i,
                                  -1. - 1.i, -6.5 + 4.5i, 12. + 1.i, -4.5 - 8.5i,
                                  16. - 1.i, -9.5 - 4.5i, -2. + 0.i, -8.5 + 5.5i,
                                  -1. + 1.i, -3.5 + 7.5i, 10. - 1.i, -5.5 - 3.5i};
    // clang-format off
    std::complex<double> vec1[4 * 4] = {2, 2, 8, 3, 3, 2, 1.+0.5i, 2, 2, 8, 3, 3, 2, 1.+0.5i,2, 2};
    // clang-format on
    std::complex<double> transformed[n * n];

    std::unique_ptr<Fourier<>> transforms[] = {
        std::make_unique<FourierBasic<>>(),
        std::make_unique<FourierFast<>>(),
        std::make_unique<FourierFastBitReverse<>>(),
        std::make_unique<FourierAdvanced<>>(),
        std::make_unique<FourierAdvancedBitReverse<>>(),
        std::make_unique<FourierFastBitReverse2dFix<>>(),
    };

    for (auto& transform : transforms) {
        transform->Forward2d(n, n, vec1, transformed);
        REQUIRE(is_identity(n * n, ans, transformed));
    }
}

TEST_CASE("test5") {
    int n = 4;
    std::complex<double> ans[4];
    std::complex<double> vec1[] = {2, 2, 8, 3};
    std::complex<double> transformed[n];

    std::unique_ptr<Fourier<>> transforms[] = {
        std::make_unique<FourierBasic<>>(),
        std::make_unique<FourierFast<>>(),
        std::make_unique<FourierFastBitReverse<>>(),
        std::make_unique<FourierAdvanced<>>(),
        std::make_unique<FourierAdvancedBitReverse<>>(),
        std::make_unique<FourierFastBitReverse2dFix<>>(),
    };

    for (auto& transform : transforms) {
        transform->Forward1d(n, vec1, transformed);
        transform->Inverse1d(n, transformed, ans);
        REQUIRE(is_identity(n, ans, vec1));
    }
}

TEST_CASE("test6") {
    int n1 = 2, n2 = 4;
    std::complex<double> ans[] = {23. + 0.5i, -4. + 0.5i, 5. + 0.5i, -4. - 1.5i,
                                  7. - 0.5i,  -8. + 1.5i, 5. - 0.5i, -8. - 0.5i};
    // clang-format off
    std::complex<double> vec1[] = {2, 2, 8,       3,
                                        3, 2, 1.+0.5i, 2};
    // clang-format on
    std::complex<double> transformed[n1 * n2];

    std::unique_ptr<Fourier<>> transforms[] = {
        std::make_unique<FourierBasic<>>(),
        std::make_unique<FourierFast<>>(),
        std::make_unique<FourierFastBitReverse<>>(),
        std::make_unique<FourierAdvanced<>>(),
        std::make_unique<FourierAdvancedBitReverse<>>(),
        std::make_unique<FourierFastBitReverse2dFix<>>(),
    };

    for (auto& transform : transforms) {
        transform->Forward2d(n1, n2, vec1, transformed);
        REQUIRE(is_identity(n1 * n2, ans, transformed));
    }
}

TEST_CASE("test7") {
    int n1 = 16, n2 = 2;
    std::complex<double> in[] = {
        0.9351048932114933, 0.734202124261979,   0.742206893695806,   0.4341030567854304,
        0.5725264269224833, 0.29584613658283165, 0.7328781224696177,  0.9400944014471981,
        0.9898853046059555, 0.2795066745708422,  0.45477857008489386, 0.9049397147415742,
        0.7438572360405848, 0.05631233021341864, 0.531495448777522,   0.5547780025653606,
        0.597249410063219,  0.5925202058586136,  0.801860738633819,   0.3041629885122463,
        0.4062551857251263, 0.5656660140126529,  0.6719279080264052,  0.6608523329991923,
        0.6280796293914833, 0.4543024360341411,  0.3821523664179747,  0.30814307574593325,
        0.9864763454520385, 0.6077610849833506,  0.17783658443208716, 0.21842993829326207};
    std::complex<double> out[n1 * n2];
    std::complex<double> transformed[n1 * n2];

    FourierBasic<>().Forward2d(n1, n2, in, out);

    std::unique_ptr<Fourier<>> transforms[] = {
        std::make_unique<FourierFast<>>(),
        std::make_unique<FourierFastBitReverse<>>(),
        std::make_unique<FourierAdvanced<>>(),
        std::make_unique<FourierAdvancedBitReverse<>>(),
        std::make_unique<FourierFastBitReverse2dFix<>>(),
    };

    for (auto& transform : transforms) {
        transform->Forward2d(n1, n2, in, transformed);
        REQUIRE(is_identity(n1 * n2, out, transformed));
    }
}
#include <benchmark/benchmark.h>

#include <complex>
#include <iostream>
#include <memory>

#include "_includes.h"

template <class Transform, int N>
static void BMTransform(benchmark::State& state) {
    std::unique_ptr<Fourier<>> transform = std::make_unique<Transform>();
    std::unique_ptr<std::complex<double>[]> in = std::make_unique<std::complex<double>[]>(N);
    std::unique_ptr<std::complex<double>[]> out = std::make_unique<std::complex<double>[]>(N);
    for (size_t i = 0; i < N; ++i) {
        in[i] = rand();
    }

    // precompute Bit-reverse permutations to not include it running time of algorithms
    BitReverse::GetInstance(BitLog(N));

    for (auto _ : state) {
        transform->Forward1d(N, in.get(), out.get());
    }
}

template <class Transform, int N1, int N2>
static void BMTransform2d(benchmark::State& state) {
    std::unique_ptr<Fourier<>> transform = std::make_unique<Transform>();
    std::unique_ptr<std::complex<double>[]> in = std::make_unique<std::complex<double>[]>(N1 * N2);
    std::unique_ptr<std::complex<double>[]> out = std::make_unique<std::complex<double>[]>(N1 * N2);
    for (size_t i = 0; i < N1 * N2; ++i) {
        in[i] = rand();
    }

    // precompute Bit-reverse permutations to not include it running time of algorithms
    BitReverse::GetInstance(BitLog(N1));
    BitReverse::GetInstance(BitLog(N2));

    for (auto _ : state) {
        transform->Forward2d(N1, N2, in.get(), out.get());
    }
}

#define SINGLE_TIME 1

#define RACE2d(Class1, Class2, n1, n2)                              \
    BENCHMARK(BMTransform2d<Class1, n1, n2>)->MinTime(SINGLE_TIME); \
    BENCHMARK(BMTransform2d<Class2, n1, n2>)->MinTime(SINGLE_TIME);

RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 128, 128);
RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 256, 256);
RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 512, 512);
RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 1024, 1024);

RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 128, 256);
RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 256, 512);
RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 512, 1024);
RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 1024, 2048);

RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 64, 512);
RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 128, 1024);
RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 256, 2048);
RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 512, 4096);

RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 512, 64);
RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 1024, 128);
RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 2048, 256);
RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 4096, 512);

// Big tests, takes long time.
#undef SINGLE_TIME
#define SINGLE_TIME 10

RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 2048, 2048);
RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 4096, 4096);
RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 8192, 8192);

RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 2048, 4096);
RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 4096, 8192);
RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 8192, 16384);

RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 1024, 8192);
RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 2048, 16384);
RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 4096, 32768);

RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 8192, 1024);
RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 16384, 2048);
RACE2d(FourierFastBitReverse2dFix<>, FourierAdvancedBitReverse<>, 32768, 4096);

BENCHMARK_MAIN();

Реализация алгоритма быстрого двумерного преобразования Фурье Тутатчикова и пр. [1], [2] (и др.) на языке С++ и сравнение времени их работы с классическим алгоритмов последовательных применений 1d БПФ сначала по одному измерению, потом по другому.

Также к обоим алгоритмам применяется оптимизация Bit-reverse permutation и в классическом алгоритме соптимизировано два транспонирования матрицы входных данных.

В пределе алгоритм Тутатчикова должен быть на 14% быстрее и бенчмарк примерно это отражает.

Результаты (в конце): https://disk.yandex.ru/i/igQPbXvcjQo8ag

/src - реализация
test.cpp - тесты с использованием библиотеки Catch2
bench.cpp - бенчмарк с использованием библиотеки google/benchmark

[1] V. S. Tutatchikov. “Two-dimensional fast Fourier transform: Batterfly in analog of CooleyTukey algorithm”. В: 2016 11th International Forum on Strategic Technology (IFOST). 2016, с. 495—498. doi: 10.1109/IFOST.2016.7884163.

[2] Mikhail Noskov и Valeriy Tutatchikov. “Modification of a two-dimensional fast Fourier transform algorithm by the analog of the Cooley-Tukey algorithm for a rectangular signal”. В: Pattern Recognition and Image Analysis 25 (янв. 2015), с. 81—83. doi: 10.1134/S1054661815010137.

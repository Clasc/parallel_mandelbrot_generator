#include <time.h>
#include <chrono>
#include <functional>
#include<vector>
#include "lib/a2/a2-helpers.hpp"

using namespace std;

double measure(function<void()> function) {
    auto t1 = chrono::high_resolution_clock::now();
    function();
    auto t2 = chrono::high_resolution_clock::now();
    return chrono::duration<double>(t2 - t1).count();
}

// A set of random gradients, adjusted for this mandelbrot algorithm
std::vector<gradient> const gradients = {
    gradient({0, 0, 0}, {76, 57, 125}, 0.0, 0.010, 2000),
    gradient({76, 57, 125}, {255, 255, 255}, 0.010, 0.020, 2000),
    gradient({255, 255, 255}, {0, 0, 0}, 0.020, 0.050, 2000),
    gradient({0, 0, 0}, {0, 0, 0}, 0.050, 1.0, 2000) };






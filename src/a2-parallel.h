#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <time.h>
#include <cmath>
#include <complex>
#include <chrono>
#include <omp.h>
#include <functional>
#include "lib/a2/a2-helpers.hpp"
using namespace std;

double measure(function<void()> function) {
    auto t1 = chrono::high_resolution_clock::now();
    function();
    auto t2 = chrono::high_resolution_clock::now();
    return chrono::duration<double>(t2 - t1).count();
}

void safeImage(Image& image, const char* filename) {
    ofstream ofs(filename, std::ofstream::out);
    ofs << "P3" << std::endl;
    ofs << image.width << " " << image.height << std::endl;
    ofs << 255 << std::endl;

    for (int j = 0; j < image.height; j++) {
        for (int i = 0; i < image.width; i++) {
            ofs << " " << image(0, j, i) << " " << image(1, j, i) << " " << image(2, j, i) << std::endl;
        }
    }
    ofs.close();
}

// A set of random gradients, adjusted for this mandelbrot algorithm
vector<gradient> const gradients = {
    gradient({0, 0, 0}, {76, 57, 125}, 0.0, 0.010, 2000),
    gradient({76, 57, 125}, {255, 255, 255}, 0.010, 0.020, 2000),
    gradient({255, 255, 255}, {0, 0, 0}, 0.020, 0.050, 2000),
    gradient({0, 0, 0}, {0, 0, 0}, 0.050, 1.0, 2000) };

// Test if point c belongs to the Mandelbrot set
bool mandelbrot_kernel(complex<double> c, vector<int>& pixel) {
    int max_iterations = 2048, iteration = 0;

    complex<double> z(0, 0);
    for (; iteration < max_iterations && abs(z) <= 4; iteration++) {
        z = z * z + c;
    }
    // now the computation of the color gradient and interpolation
    double log_len = log(sqrt(z.real() * z.real() + z.imag() * z.imag()));
    double log_2 = log(2.0);
    long double m = (iteration + 1 - log_len / log_2);
    double q = m / (double)max_iterations;

    q = iteration + 1 - log(log_len) / log_2;
    q /= max_iterations;

    colorize(pixel, q, iteration, gradients);

    return (iteration < max_iterations);
}

/**
 * Compute the Mandelbrot set for each pixel of a given image.
 * Image is the Image data structure for storing RGB image
 * The default value for ratio is 0.15.
 *
 * @param[inout] image
 * @param[in] ratio
 *
*/
int mandelbrot(Image& image, double ratio = 0.15) {
    // reduction: gives each thread a private pixels_inside variable that is summed at the end
    int pixels_inside = 0;
    ratio /= 10.0;
    int channels = image.channels;
    int w = image.width;
    int h = image.height;
    vector<int> pixel = { 0, 0, 0 }; // red, green, blue (each range 0-255)
    complex<double> c;

    #pragma omp parallel shared(image, ratio, pixels_inside) private (pixel, c)
    {
        #pragma omp single
        {
            #pragma omp taskgroup
            {
                for (int j = 0; j < h; j++) {

                    #pragma omp task shared(w, h, channels, image, ratio, pixels_inside)
                    {
                        for (int i = 0; i < w; i++) {

                            double dx = (double)i / (w)*ratio - 1.10;
                            double dy = (double)j / (h) * 0.1 - 0.35;

                            c = complex<double>(dx, dy);

                            // the actual mandelbrot kernel
                            if (mandelbrot_kernel(c, pixel)) {
                                #pragma omp atomic
                                pixels_inside++;
                            }

                            // apply to the image
                            for (int ch = 0; ch < channels; ch++)
                                image(ch, j, i) = pixel[ch];
                        }
                    }
                }
            }
        }
    }
    return pixels_inside;
}


/**
 * 2D Convolution
 * src is the source Image to which we apply the filter.
 * Resulting image is saved in dst. The size of the kernel is
 * given with kernel_width (must be odd number). Sigma represents
 * the standard deviation of the filter. The number of iterations
 * is given with the nstep (default=1)
 *
 * @param[in] src
 * @param[out] dst
 * @param[in] kernel_width
 * @param[in] sigma
 * @param[in] nsteps
 *
*/
void convolution_2d(Image& src, Image& dst, int kernel_width, double sigma, int nsteps = 1) {
    int h = src.height;
    int w = src.width;
    int channels = src.channels;

    std::vector<std::vector<double>> kernel = get_2d_kernel(kernel_width, kernel_width, sigma);

    int displ = (kernel.size() / 2); // height==width!
    for (int step = 0; step < nsteps; step++) {
        for (int ch = 0; ch < channels; ch++) {
            for (int i = 0; i < h; i++) {
                for (int j = 0; j < w; j++) {
                    double val = 0.0;

                    for (int k = -displ; k <= displ; k++) {
                        for (int l = -displ; l <= displ; l++) {
                            int cy = i + k;
                            int cx = j + l;
                            int src_val = 0;

                            // if it goes outside we disregard that value
                            if (cx < 0 || cx > w - 1 || cy < 0 || cy > h - 1) {
                                continue;
                            }
                            else {
                                src_val = src(ch, cy, cx);
                            }

                            val += kernel[k + displ][l + displ] * src_val;
                        }
                    }
                    dst(ch, i, j) = (int)(val > 255 ? 255 : (val < 0 ? 0 : val));
                }
            }
        }

        if (step < nsteps - 1) {
            // swap references
            // we can reuse the src buffer for this example
            Image tmp = src;
            src = dst;
            dst = tmp;
        }
    }
}

void generate(Image& image, int n_steps = 20, const char* outputName = "Mandelbrot_parallel.ppm") {
    Image filtered_image(image.channels, image.height, image.width);
    auto ratio = image.width / (double)image.height;
    int pixels_inside = 0;
    auto duration = measure([&]() {
        pixels_inside = mandelbrot(image, ratio);
        });

    // TODO Use OpenMP tasking to implement a parallel version
    cout << "Mandelbrot time: " << duration << endl;
    cout << "Total Mandelbrot pixels: " << pixels_inside << endl;

    // Actual 2D convolution part
    // Use OpenMP tasking to implement a parallel version

    auto duration2 = measure([&]() { convolution_2d(image, filtered_image, 5, 0.37, n_steps); });

    cout << "Convolution time: " << duration2 << endl;
    cout << "Total time: " << duration + duration2 << endl;
    safeImage(filtered_image, outputName);
}

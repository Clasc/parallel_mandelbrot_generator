#include <vector>
#include <tuple>
#include <cmath>
#include <omp.h>
#include "ParallelGenerator.h"

using namespace std;

/**
 * @brief Use only Omp tasks for parallelization
 *
 */
class TasksGenerator : public ParallelGenerator {
protected:
    char const* getName() override { return "TasksGenerator"; }

public:
    TasksGenerator(int max_tasks) :ParallelGenerator(max_tasks) {}

    bool mandelbrot_kernel(complex<double> c, vector<int>& pixel) override {
        int max_iterations = 2048, iteration = 0;

        complex<double> z(0, 0);

        for (; iteration < max_iterations && abs(z) <= 4; iteration++) {
            z = z * z + c;
        }
        // now the computation of the color gradient and interpolation
        double log_len = log(sqrt(z.real() * z.real() + z.imag() * z.imag()));
        double log_2 = log(2.0);
        long double m = (iteration + 1 - log_len / log_2);

        auto q = iteration + 1 - log(log_len) / log_2;
        q /= max_iterations;

        colorize(pixel, q, iteration, gradients);

        return (iteration < max_iterations);
    }

    int mandelbrot(Image& image, double ratio = 0.15)override {
        int pixels_inside = 0;
        ratio /= 10.0;
        vector<int> pixel = { 0, 0, 0 }; // red, green, blue (each range 0-255)
        complex<double> c;

        #pragma omp parallel 
        {
            #pragma omp single nowait
            {
                #pragma omp taskgroup
                {
                    for (int j = 0; j < image.height; j++) {
                        double dy = (double)j / (image.height) * 0.1 - 0.35;
                        #pragma omp task shared(image, ratio, pixels_inside) firstprivate(dy) private(pixel, c) untied final(j > max_tasks)
                        {
                            for (int i = 0; i < image.width; i++) {
                                double dx = (double)i / (image.width) * ratio - 1.10;

                                c = complex<double>(dx, dy);

                                // the actual mandelbrot kernel
                                if (mandelbrot_kernel(c, pixel)) {
                                    #pragma omp atomic
                                    pixels_inside++;
                                }

                                // apply to the image
                                for (int ch = 0; ch < image.channels; ch++)
                                    image(ch, j, i) = pixel[ch];
                            }
                        }
                    }
                }
            }
        }
        return pixels_inside;
    }

    void convolution_2d(Image& src, Image& dst, int kernel_width, double sigma, int nsteps = 1) override {
        std::vector<std::vector<double>> kernel = get_2d_kernel(kernel_width, kernel_width, sigma);
        int displ = (kernel.size() / 2); // height==width!

        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                for (int step = 0; step < nsteps; step++) {
                    #pragma omp taskgroup
                    {
                        for (int ch = 0; ch < src.channels; ch++) {
                            for (int i = 0; i < src.height; i++) {
                                #pragma omp task shared(src, dst, displ, ch, i, step, kernel) untied final(i > max_tasks)
                                {
                                    for (int j = 0; j < src.width; j++) {
                                        double val = 0.0;

                                        for (int k = -displ; k <= displ; k++) {
                                            for (int l = -displ; l <= displ; l++) {
                                                int cy = i + k;
                                                int cx = j + l;
                                                int src_val = 0;

                                                // if it goes outside we disregard that value
                                                if (cx < 0 || cx > src.width - 1 || cy < 0 || cy > src.height - 1) {
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
        }
    }
};

/**
 * @brief Use Omp taskloops for parallelization
 *
 */
class TaskloopGenerator : public ParallelGenerator {
protected:
    char const* getName() override { return "TaskloopGenerator"; }

public:
    TaskloopGenerator(int max_tasks) :ParallelGenerator(max_tasks) {}

    bool mandelbrot_kernel(complex<double> c, vector<int>& pixel) override {
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

    int mandelbrot(Image& image, double ratio = 0.15) override {
        int pixels_inside = 0;
        ratio /= 10.0;
        vector<int> pixel = { 0, 0, 0 }; // red, green, blue (each range 0-255)
        complex<double> c;

        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                #pragma omp taskloop shared(image, ratio) private(pixel, c) num_tasks(max_tasks) reduction(+:pixels_inside) untied
                for (int j = 0; j < image.height; j++) {
                    #pragma omp taskloop shared(image, ratio) private(pixel, c) num_tasks(max_tasks) reduction(+:pixels_inside) untied
                    for (int i = 0; i < image.width; i++) {

                        double dx = (double)i / (image.width) * ratio - 1.10;
                        double dy = (double)j / (image.height) * 0.1 - 0.35;

                        c = complex<double>(dx, dy);

                        // the actual mandelbrot kernel
                        if (mandelbrot_kernel(c, pixel)) {
                            pixels_inside++;
                        }

                        // apply to the image
                        for (int ch = 0; ch < image.channels; ch++)
                            image(ch, j, i) = pixel[ch];
                    }
                }
            }
        }
        return pixels_inside;
    }

    void convolution_2d(Image& src, Image& dst, int kernel_width, double sigma, int nsteps = 1) override {
        std::vector<std::vector<double>> kernel = get_2d_kernel(kernel_width, kernel_width, sigma);
        int displ = (kernel.size() / 2); // height==width!
        #pragma omp parallel shared(src, dst, displ)
        {
            #pragma omp single nowait
            {
                for (int step = 0; step < nsteps; step++) {
                    for (int ch = 0; ch < src.channels; ch++) {
                        #pragma omp taskloop num_tasks(max_tasks) shared(src, dst, displ, ch, step, kernel) untied
                        for (int i = 0; i < src.height; i++) {
                            for (int j = 0; j < src.width; j++) {
                                double val = 0.0;
                                for (int k = -displ; k <= displ; k++) {
                                    for (int l = -displ; l <= displ; l++) {
                                        int cy = i + k;
                                        int cx = j + l;
                                        int src_val = 0;

                                        // if it goes outside we disregard that value
                                        if (cx < 0 || cx > src.width - 1 || cy < 0 || cy > src.height - 1) {
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
        }
    }
};
#include <complex>
#include <iostream>
#include <fstream>
#include "a2-utils.h"
/**
 * Superclass for parallelgenerators of the mandelbrot set
 * */
class ParallelGenerator {
private:
public:

    /**
     * @brief  prepares the algorithm, makes the 2d convolution, measures the performance and saves the image to a file
     *
     * @param image The input image
     * @param outputName name of the outputfile
     * @param n_steps number of steps for the convolution
     */
    void generate(Image image, const char* outputName = "Mandelbrot_parallel.ppm", int n_steps = 20) {
        Image filtered_image(image.channels, image.height, image.width);
        auto ratio = image.width / (double)image.height;
        int pixels_inside = 0;
        auto duration = measure([&]() {
            pixels_inside = mandelbrot(image, ratio);
            });

        cout << "Mandelbrot time: " << duration << endl;
        cout << "Total Mandelbrot pixels: " << pixels_inside << endl;

        // Actual 2D convolution part
        auto duration2 = measure([&]() { convolution_2d(image, filtered_image, 5, 0.37, n_steps); });

        cout << "Convolution time: " << duration2 << endl;
        cout << "Total time: " << duration + duration2 << endl;
        safeImage(filtered_image, outputName);
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


    /**
     * @brief  Test if point c belongs to the Mandelbrot set
     */
    virtual bool mandelbrot_kernel(std::complex<double> c, std::vector<int>& pixel) = 0;

    /**
     * Compute the Mandelbrot set for each pixel of a given image.
     * Image is the Image data structure for storing RGB image
     * The default value for ratio is 0.15.
     *
     * @param[inout] image
     * @param[in] ratio
     *
    */
    virtual int mandelbrot(Image& image, double ratio = 0.15) = 0;

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
    virtual void convolution_2d(Image& src, Image& dst, int kernel_width, double sigma, int nsteps = 1) = 0;
};
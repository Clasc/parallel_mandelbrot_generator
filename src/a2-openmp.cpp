#include "TasksGenerator.h"

int main(int argc, char const* argv[]) {
    auto image = Image(3, 1024, 1536);
    omp_set_num_threads(16);
    #pragma omp parallel
    {
        #pragma omp single
        {
            cout << "threads:" << omp_get_num_threads() << endl;
        }
    }
    auto tasks = TasksGenerator();
    auto taskloops = TaskloopGenerator();

    tasks.generate(image);
    taskloops.generate(image, "Mandelbrot-taskloop.ppm");

    return 0;
}

#include "TasksGenerator.h"
#include <cstdlib>

int main(int argc, char const* argv[]) {
    auto image = Image(3, 1024, 1536);
    int threads;
    int tasks_per_thread;

    std::cout << "argc" << argc << std::endl;

    if (argc >= 3) {
        std::istringstream input_threads(argv[1]);
        std::istringstream input_max_t(argv[2]);

        if (!(input_threads >> threads)) {
            threads = 16;
        }

        if (!(input_max_t >> tasks_per_thread)) {
            tasks_per_thread = 4;
        }
    }

    omp_set_num_threads(threads);

    #pragma omp parallel
    {
        #pragma omp single
        {
            cout << "threads:" << omp_get_num_threads() << endl << "tasks_per_thread:" << tasks_per_thread << endl;
        }
    }

    auto tasks = TasksGenerator(threads * tasks_per_thread);
    auto taskloops = TaskloopGenerator(threads * tasks_per_thread);

    tasks.generate(image);
    taskloops.generate(image, "Mandelbrot-taskloop.ppm");

    return 0;
}

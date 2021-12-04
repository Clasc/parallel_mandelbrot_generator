#include "TasksGenerator.h"

int main(int argc, char const* argv[]) {
    auto image = Image(3, 1024, 1536);
    auto tasks = TasksGenerator();
    auto taskloops = TaskloopGenerator();

    tasks.generate(image);
    taskloops.generate(image, "Mandelbrot-taskloop.ppm");

    return 0;
}

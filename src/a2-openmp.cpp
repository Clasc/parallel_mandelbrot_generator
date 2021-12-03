#include "TasksGenerator.h"

int main(int argc, char const* argv[]) {
    auto image = Image(3, 1024, 1536);
    auto tasks = TasksGenerator();
    auto taskloops = TaskloopGenerator();
    std::cout << "Generating by Tasks..." << std::endl;
    tasks.generate(image);
    std::cout << "Generating by Taskloop..." << std::endl;
    taskloops.generate(image, "Mandelbrot-taskloop.ppm");
    return 0;
}

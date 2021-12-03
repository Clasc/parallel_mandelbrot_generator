#include "TasksGenerator.h"

int main(int argc, char const* argv[]) {
    auto image = Image(3, 1024, 1536);
    auto gen = TasksGenerator();
    gen.generate(image);
    return 0;
}

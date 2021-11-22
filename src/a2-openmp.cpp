#include "a2-parallel.h"

int main(int argc, char const* argv[]) {
    auto image = Image(3, 1024, 1536);
    generate(image);
    return 0;
}

#include "a2-parallel.h"

int main(int argc, char const *argv[])
{
    auto x = 5;
    auto y = 23;
    auto z = x + y;
    printf("%i \n", z);
    generate(Image(3, 1024, 1536));
    return 0;
}

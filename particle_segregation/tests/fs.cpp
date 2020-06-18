#include <iostream>
#include <vector>
#include "filesystem.hpp"

int main()
{
    util::IO io_object;
    io_object.scan_config();
    io_object.dump(std::vector<float>{1,2,3,4});

    return 0;
}

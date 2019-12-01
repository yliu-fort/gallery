#include <iostream>
#include <string>
#include "cmake_source_dir.h"

int main()
{
    printf("Test macro function.\n");
    printf("file path: %s \n", INFILE(file.dat));
    printf("file path: %s \n", std::string("/root/").append("file.dat").c_str());
    printf("file path: %s \n", FP("file.dat"));
    return 0;
}

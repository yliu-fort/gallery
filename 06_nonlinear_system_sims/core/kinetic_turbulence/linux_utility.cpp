#include <bits/stdc++.h>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include "linux_utility.h"
void checkDir(const char* path)
{
    // Creating a directory
    if (mkdir(path, 0777) == -1)
        std::cerr << "Mkdir::Error :  " << strerror(errno) << std::endl;

    else
        std::cout << "Directory " << path << " created\n";
}

#include <iostream>
#include <fstream>
#include "cmake_source_dir.h"

template <typename T>
void SwapEnd(T& var)
{
  char* varArray = reinterpret_cast<char*>(&var);
  for(long i = 0; i < static_cast<long>(sizeof(var)/2); i++)
    std::swap(varArray[sizeof(var) - 1 - i],varArray[i]);
}

int main()
{
    std::ifstream is (FP("test.bin.vtk"), std::ifstream::binary);
    if (is) {
      // get length of file:
      is.seekg (0, is.end);
      int length = is.tellg();
      is.seekg (0, is.beg);

      char * buffer = new char [length];

      //std::cout << "Reading " << length << " characters... ";
      // read data as a block:
      //is.read (buffer,length);

      //vtkstream<<"# vtk DataFile Version 2.0"<<"\n";
      //vtkstream<<"Exemple"<<"\n";
      //vtkstream<<"BINARY"<<"\n";
      //vtkstream<<"DATASET STRUCTURED_GRID"<<std::endl;
      //vtkstream<<"DIMENSIONS 6 4 1"<<std::endl;
      //vtkstream<<"POINTS 24 double"<<std::endl;

      char header[255];
      is.getline(header,255);
      std::cout << header << std::endl;

      is.getline(header,255);
      std::cout << header << std::endl;

      is.getline(header,255);
      std::cout << header << std::endl;

      is.getline(header,255);
      std::cout << header << std::endl;

      is.getline(header,255);
      char ss[80];int nx,ny,nz;
      sscanf(header, "%s %d %d %d", ss, &nx, &ny, &nz);
      std::cout << ss << " " << nx  << " " << ny << " " << nz << std::endl;

      is.getline(header,255);
      char ft[80];char ftype[80];int nelem;
      sscanf(header, "%s %d %s", ft, &nelem, ftype);
      std::cout << ft << " " << nelem << " " << ftype << std::endl;

      double myarray[72];
      for (unsigned int i = 0; i < 72; ++i) {
        is.read(reinterpret_cast<char*>(&myarray[i]), sizeof(double));
        SwapEnd(myarray[i]);
        std::cout << myarray[i] << std::endl;
      }

      if (is)
        std::cout << "all characters read successfully.";
      else
        std::cout << "error: only " << is.gcount() << " could be read";
      is.close();

      // ...buffer contains the entire file...

      delete[] buffer;
    }
    return 0;
}

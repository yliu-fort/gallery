#include <iostream>
#include <fstream>
#include "cmake_source_dir.h"

// Thanks to https://stackoverflow.com/questions/105252
template <typename T>
void SwapEnd(T& var)
{
  char* varArray = reinterpret_cast<char*>(&var);
  for(long i = 0; i < static_cast<long>(sizeof(var)/2); i++)
    std::swap(varArray[sizeof(var) - 1 - i],varArray[i]);
}

double myarray[72] = {
  0.001,0.002,0,1,0,0,2,0,0,3,0,0,4,0,0,
  5,0,0,0,1,0,1,1,0,2,1,0,3,1,0,
  4,1,0,5,1,0,0,2,0,1,2,0,2,2,0,
  3,2,0,4,2,0,5,2,0,0,3,0,1,3,0,
  2,3,0,3,3,0,4,3,0,5,3,0};

int main()
{
  std::ofstream vtkstream;
  vtkstream.open(FP("test.bin.vtk"), std::ios::out | std::ios::app | std::ios::binary);
  if (vtkstream) {
    vtkstream<<"# vtk DataFile Version 2.0"<<"\n";
    vtkstream<<"Exemple"<<"\n";
    vtkstream<<"BINARY"<<"\n";
    vtkstream<<"DATASET STRUCTURED_GRID"<<std::endl;
    vtkstream<<"DIMENSIONS 6 4 1"<<std::endl;
    vtkstream<<"POINTS 24 double"<<std::endl;
    for (unsigned int i = 0; i < 72; ++i) {
      SwapEnd(myarray[i]);
      vtkstream.write((char*)&myarray[i], sizeof(double));
    }
    vtkstream.close();
  } else {
    std::cout<<"ERROR"<<std::endl;
  }
  return 0;
}

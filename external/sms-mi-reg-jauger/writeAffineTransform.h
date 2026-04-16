#ifndef _WRITE_AFFINE_TRANSFORM_INCLUDE
#define _WRITE_AFFINE_TRANSFORM_INCLUDE 1

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <itkMacro.h>
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkAffineTransform.h>

#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>

inline bool file_exists_ifstream(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

inline bool file_exists(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

#endif


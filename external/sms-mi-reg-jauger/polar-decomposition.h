#ifndef _POLAR_DECOMPOSITION_INCLUDE
#define _POLAR_DECOMPOSITION_INCLUDE 1

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkAffineTransform.h>

#include <itkCenteredTransformInitializer.h>
#include <itkCenteredVersorTransformInitializer.h>
#include <itkVersor.h>
#include <itkVersorRigid3DTransform.h>

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


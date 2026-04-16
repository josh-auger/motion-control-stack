#include <cmath>
#include <itkMacro.h>
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkAffineTransform.h>
#include <itkVersor.h>
#include <itkMatrix.h>
#include <iomanip>
#include <iostream>
#include <vector>

// Include the Eigen template library for math
#include ITK_EIGEN(Core)
#include ITK_EIGEN(Dense)
#include ITK_EIGEN(LU)

typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> MatrixType33d;

#include "writeAffineTransform.h"

int main(int argc, char *argv[])
{

  using ScalarType = double;
  constexpr unsigned int Dimension = 3;

  //New Argument error check
  printf("The value of argc is %03d\n",argc);
  if (argc <=4) {
    std::cout << argv[0] << " affineTransformFile angle scale shear  " << std::endl;
    exit(1);
  }
 
  std::string affineTransformFile = argv[1];
  double angle = std::stod(argv[2]);
  double scale = std::stod(argv[3]);
  double shear = std::stod(argv[4]);

  using AffineTransformType = itk::AffineTransform<ScalarType, Dimension>;
  auto affineTransform = AffineTransformType::New();

  affineTransform->Print( std::cout );

  std::cout << "Rotate by " << angle << std::endl;
  affineTransform->Rotate(0, 1, angle);

  std::cout << "Scale by " << scale << std::endl;
  affineTransform->Scale(scale);

  std::cout << "Shear by " << shear << std::endl;
  affineTransform->Shear( 0, 1, shear );

  affineTransform->Print( std::cout );
  AffineTransformType::ParametersType parms = affineTransform->GetParameters();

  using TransformFileWriterType = itk::TransformFileWriterTemplate<ScalarType>;
  auto writer = TransformFileWriterType::New();
  writer->SetInput( affineTransform );
  writer->SetFileName( "affineTransform.tfm" );

  try {
    writer->Update();
  } catch (const itk::ExceptionObject & excp) {
    std::cerr << "Error while writing the affine transform." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }


  MatrixType33d M;
  M = MatrixType33d::Identity();
  // Now assign the contents of M based on the affine transform that has been read.
  for (unsigned int k = 0; k < 3; k++) {
    for (unsigned int j = 0; j < 3; j++) {
      M(k,j) = parms[k*3 + j];
    }
  }

  // At this point in the code the Eigen matrix M has been initialized based on the parameters from the affine transform.
  MatrixType33d PQ = M;
  MatrixType33d NQ = M;
  MatrixType33d PQNQDiff = MatrixType33d::Identity();
  MatrixType33d PQit = MatrixType33d::Identity();
  MatrixType33d PQi = MatrixType33d::Identity();

  const unsigned int maximumIterations = 100;
  const double eps = 1e-11;

  // Calculate part of the polar decomposition of the matrix M.
  // Compute the matrix NQ as the best orthogonal approximation to the input matrix M.
  for(unsigned int ni = 0; ni < maximumIterations; ni++ )
  {
    std::cout << "PQ at iteration " << ni << " is : " << std::endl << PQ << std::endl;
    // Average current Qi with its inverse transpose
    PQi = PQ.inverse();
    std::cout << "PQ inverse at iteration " << ni << " is : " << std::endl << PQi << std::endl;
    PQit = PQi.transpose();
    std::cout << "PQ inverse tranpose at iteration " << ni << " is : " << std::endl << PQit << std::endl;
    NQ = ( PQ + PQit) / 2.0;
    std::cout << "NQ at iteration " << ni << " is : " << std::endl << NQ << std::endl;
    PQNQDiff = NQ - PQ;
    if( PQNQDiff.norm() < eps ) {
      std::cout << "Polar decomposition used " << ni << " iterations " << std::endl;
      break;
    } else {
      PQ = NQ;
    }
  }

  MatrixType33d QMatrix;
  QMatrix = NQ;
  std::cout << "Initial Matrix = " << std::endl << M << std::endl;
  std::cout << "Q Matrix = " << std::endl << QMatrix << std::endl;

  // A versor can be initialized from an orthogonal matrix.
  // https://itk.org/Doxygen/html/classitk_1_1Versor.html
  itk::Versor<double>::MatrixType matrix;
  itk::Versor<double> rotation;

  // Initialize the itk::Matrix with the newly computed rotation matrix
  for (unsigned int k = 0; k < 3; k++) {
    for (unsigned int j = 0; j < 3; j++) {
      matrix(k,j) = QMatrix(k,j);
    }
  }
        
  // Initialize the versor from the matrix.
  rotation.Set( matrix );

  std::cout << " Versor transform is : " << rotation << std::endl;

  exit(0);
}


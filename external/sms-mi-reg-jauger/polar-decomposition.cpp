#include <cmath>
#include <itkMacro.h>
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkAffineTransform.h>
#include <itkVersorRigid3DTransform.h>
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

#include "polar-decomposition.h"

int main(int argc, char *argv[])
{

  using ScalarType = double;
  constexpr unsigned int Dimension = 3;
  using VersorRigid3DTransformType = itk::VersorRigid3DTransform< ScalarType >;
  using AffineTransformType = itk::AffineTransform<ScalarType, Dimension>;
  // auto affineTransform = AffineTransformType::New();


  //New Argument error check
  if (argc <=2) {
    std::cout << argv[0] << " affineTransform rotationTransform" << std::endl;;
    printf("The value of argc is %03d\n",argc);
    exit(1);
  }
 
  std::string affineTransformFile = argv[1];
  std::string rotationTransformFile = argv[2];

  itk::TransformFileReader::Pointer trsfReader;
  trsfReader = itk::TransformFileReader::New();
  trsfReader->SetFileName( affineTransformFile );
  
  if (!file_exists( affineTransformFile )) {
    std::cout << "There is no input file called " << affineTransformFile << std::endl;
  } else {

    try {
      trsfReader->Update();
    } catch ( itk::ExceptionObject & excp ) {
      // Display error from reading the transform file.
      std::cerr << "Error while reading the transform file " <<
         affineTransformFile << std::endl;
      std::cerr << excp << std::endl;
      std::cerr << "[FAILED]" << std::endl;
      return EXIT_FAILURE;
    }
  }

 // Now try to work out how many and what type of transforms were read.
 // We only want to get one transform.
  VersorRigid3DTransformType::Pointer vrTransform =
    VersorRigid3DTransformType::New(); // This is the transform to be read in.
  VersorRigid3DTransformType::Pointer trsf_read =
    VersorRigid3DTransformType::New(); // This is the transform to be read in.
  AffineTransformType::Pointer affine_read =
    AffineTransformType::New(); // This is the transform that may be read in.

  VersorRigid3DTransformType::Pointer finalVRTransform = VersorRigid3DTransformType::New(); // This is the transform that is constructed.
  VersorRigid3DTransformType::ParametersType parmArray( finalVRTransform->GetNumberOfParameters() );

  VersorRigid3DTransformType::Pointer identityVRTransform = VersorRigid3DTransformType::New(); // This is to hold an identity transform.

  typedef itk::TransformFileReader::TransformListType * TransformListType;
  TransformListType transforms = trsfReader->GetModifiableTransformList();
  std::cout << "Number of transforms = " << transforms->size() << std::endl;
  itk::TransformFileReader::TransformListType::const_iterator it =
      transforms->begin();
  if (transforms->size() <= 0 || transforms->size() > 1) {
      std::cerr << "Read " << transforms->size() << " transforms but want 1." << std::endl;
      return EXIT_FAILURE;
  }

  itk::TransformFileReader::TransformListType::const_iterator trsfit =
           transforms->begin();
  if (!strcmp((*trsfit)->GetNameOfClass(),"VersorRigid3DTransform")) {
   trsf_read = static_cast<VersorRigid3DTransformType*>((*trsfit).GetPointer());
    std::cout << "Just loaded a versor rigid 3d transform." << std::endl;
    trsf_read->Print(std::cout);

    // This is the recommended procedure for initializing when the
    // transforms have different mappings of actual parameter to index
    // in the array of parameters.
    vrTransform->SetCenter( trsf_read->GetCenter() );
    vrTransform->SetTranslation( trsf_read->GetTranslation() );
        // API BUG: itk::VersorRigid3DTransform<double>' has no member
        //   named 'GetRotation'
        // API BUG: ScaleSkewVersor3DTransform has no SetVersor()
    vrTransform->SetRotation( trsf_read->GetVersor() );
  } else if (!strcmp((*trsfit)->GetNameOfClass(),"AffineTransform")) {
        affine_read = static_cast<AffineTransformType*>((*trsfit).GetPointer());
        std::cout << "Just loaded an affine transform." << std::endl;
        affine_read->Print(std::cout);
       // This is the recommended procedure for initializing when the
        // transforms have different mappings of actual parameter to index
        // in the array of parameters, as is the case here.
        vrTransform->SetCenter( affine_read->GetCenter() );
        vrTransform->SetTranslation( affine_read->GetTranslation() );
        // There is no built in means to infer the rotation implied by the
        // transformation, so we have to compute it.
        // Shoemake in the paper "Matrix Animation and Polar Decomposition".
        // http://citeseer.ist.psu.edu/shoemake92matrix.html
        // https://research.cs.wisc.edu/graphics/Courses/838-s2002/Papers/polar-decomp.pdf

        // Switched to using Eigen3 for the matrix operations
        // ITK-v5.3.0/Modules/ThirdParty/Eigen3

        // Get the 3x3 set of parameters.
        // There is also a set of fixed parameters that represent the origin of the transformation.

        AffineTransformType::ParametersType parms = affine_read->GetParameters();

        MatrixType33d M;
        M = MatrixType33d::Identity();
        // Now assign M from the affine transform that has been read.
        for (unsigned int k = 0; k < 3; k++) {
          for (unsigned int j = 0; j < 3; j++) {
            M(k,j) = parms[k*3 + j];
          }
        }

        // At this point in the code the Eigen matrix M has been initialized based on the parameters from the affine transform.
        MatrixType33d PQ = M;
        MatrixType33d NQ = M;
        MatrixType33d PQNQDiff = MatrixType33d::Identity();
        MatrixType33d PQi = MatrixType33d::Identity();
        MatrixType33d PQit = MatrixType33d::Identity();

        const unsigned int maximumIterations = 100;
        const double eps = 1e-11;

        // Calculate part of the polar decomposition of the matrix M.
        // Compute the matrix NQ as the best orthogonal approximation to 
        // the input matrix M.
        for(unsigned int ni = 0; ni < maximumIterations; ni++ )
        {
          // Average current Qi with its inverse transpose
          PQi = PQ.inverse();
          PQit = PQi.transpose();
          NQ = ( PQ + PQit) / 2.0;
          PQNQDiff = NQ - PQ;
          if( PQNQDiff.norm() < eps )
          {
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

        // Set the versor part of the rigid body transform using the versor
        vrTransform->SetRotation( rotation );

  } else {
        std::cerr << "Can't initialize from transform of type " <<
                ( (*trsfit)->GetNameOfClass() ) << " ." << std::endl;
      return EXIT_FAILURE;
  }

  parmArray[0] = vrTransform->GetParameters()[0];
  parmArray[1] = vrTransform->GetParameters()[1];
  parmArray[2] = vrTransform->GetParameters()[2];
  parmArray[3] = vrTransform->GetParameters()[3];
  parmArray[4] = vrTransform->GetParameters()[4];
  parmArray[5] = vrTransform->GetParameters()[5];
  double maxf = 0.0;

  std::cout << "Now setting the parameters of the transform." << std::endl 
          << std::flush;
  finalVRTransform->SetParameters( parmArray );

  std::cout << "Assigned the parameter array to the finalVRTransform." 
      << std::endl << std::flush;

std::cout << "Now writing out the final transform." << std::endl << std::flush;

  // Now write out the final VR Transform
  itk::TransformFileWriter::Pointer trsfWriter;
  trsfWriter = itk::TransformFileWriter::New();
  trsfWriter->SetInput( finalVRTransform );
  trsfWriter->SetFileName( rotationTransformFile );

  try
  {
    trsfWriter->Update();
  }
  catch( itk::ExceptionObject & excp)
  {
    std::cerr << "Error while saving the transform to file "
              << std::endl;
    std::cerr << excp << std::endl;
    exit(1);
  }


  exit(0);
}


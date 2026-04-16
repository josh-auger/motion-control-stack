#include <cmath>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <iomanip>
#include <iostream>
#include <vector>

#include "happly.h"


int main(int argc, char *argv[])
{

  const unsigned int Dimension = 3;
  typedef double PixelType;
  typedef itk::Image< PixelType,  Dimension >   ImageType;
  typedef itk::ImageFileReader< ImageType  >  ReaderType;
  typedef itk::ImageFileWriter< ImageType  >  WriterType;

  //New Argument error check
  printf("The value of argc is %03d\n",argc);
  if (argc <= 2) {
    std::cout << "Usage: " << argv[0] << " inputFile.nrrd outputFile.ply" << std::endl;
    exit(1);
  }
 
  //adding new arguments
  std::string referencevolume = argv[1];
  std::string outputPlyFileName = argv[2];

  ReaderType::Pointer readerRefVol = ReaderType::New();

  readerRefVol->SetFileName( referencevolume.c_str() );

  try {
    readerRefVol->Update();
  } catch ( itk::ExceptionObject & excp )
  {
    // Display error from reading the reference volume file.
    std::cerr << "Error while reading the reference volume " <<
                     (referencevolume.c_str()) << std::endl;
    std::cerr << excp << std::endl;
    std::cerr << "[FAILED]" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Successfully read in the image : " << referencevolume << std::endl;


  // Iterate over the input image, finding the coordinates of every voxel.
  itk::ImageRegionIteratorWithIndex<ImageType> imageIterator(
      readerRefVol->GetOutput(), 
      readerRefVol->GetOutput()->GetLargestPossibleRegion() 
      );

  std::vector< std::array<double,3> > vertexPositions;
  std::vector<double> vertexColor;

  std::array<double,3> vertexPosition;

  imageIterator.GoToBegin();
  while (!imageIterator.IsAtEnd()) {
    ImageType::IndexType idx = imageIterator.GetIndex();
    ImageType::PointType pointInVolume;

    double voxelValue = readerRefVol->GetOutput()->GetPixel(idx);
    readerRefVol->GetOutput()->TransformIndexToPhysicalPoint( idx, pointInVolume );
    // std::cout << "The point " << pointInVolume << " has value : " <<  voxelValue << std::endl;

    vertexPosition[0] = pointInVolume[0];
    vertexPosition[1] = pointInVolume[1];
    vertexPosition[2] = pointInVolume[2];
    vertexPositions.push_back( vertexPosition );
    vertexColor.push_back( voxelValue );

    ++imageIterator;
  }

  // Create an empty object
  happly::PLYData plyOut;
  plyOut.addVertexPositions( vertexPositions );
  plyOut.getElement("vertex").addProperty<double>("intensity", vertexColor);
  

  // Write the object to file
  plyOut.write( outputPlyFileName, happly::DataFormat::ASCII);


  exit(0);
}


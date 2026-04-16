#ifndef _SMS_MI_REG_INCLUDE
#define _SMS_MI_REG_INCLUDE 1

#include <cmath>
#include <iostream>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <itkImage.h>
#include <itkOrientImageFilter.h>
#include <itkBSplineInterpolateImageFunction.h>
#include "itkLinearInterpolateImageFunction.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkAffineTransform.h>
#include <itkPoint.h>
#include <itkContinuousIndex.h>

#include <itkCenteredTransformInitializer.h>
#include <itkCenteredVersorTransformInitializer.h>
#include <itkVersor.h>
#include <itkVersorRigid3DTransform.h>

#include <itkPowellOptimizer.h>

// static allocation of arrays causes a compile time limit on the maximum
// number of threads that can be used.
const int MAX_NUM_THREADS = 64;
// const int MAX_NUM_THREADS = 1;

// Parameters associated with resampling volume to match a slice.
const int ARRAY_SIZE = 4096;
// pad to avoid false sharing of cache line.
const unsigned int padBest = 64/sizeof(float);// Should avoid false sharing.
const unsigned int padWorst = 0; // Generates false sharing of cache lines.
const unsigned int pad = padBest;
    // We actually eliminated the resampling to this
    // array n[][] in favor of computing the joint histograms, and then
    // computing MI.
float n[ARRAY_SIZE][pad]; 

// Parameters associated with computing mutual information
// We need one histogram for each thread as they can't be written to
// in parallel.
const int MaxNumBins = 66;
float jointHistograms[MAX_NUM_THREADS*MaxNumBins*MaxNumBins][pad];
float marginalFixed[MaxNumBins][pad];
float marginalMoving[MaxNumBins][pad];
float sumHistogram[MaxNumBins*MaxNumBins][pad];
float (*sumHistogramPtr)[MaxNumBins*MaxNumBins][pad] = &sumHistogram;

const bool useDirichletPrior = true;

const unsigned int Dimension = 3;
typedef float PixelType;
typedef itk::Image< PixelType,  Dimension >   ImageType;
typedef itk::ImageFileReader< ImageType  >  ReaderType;
typedef itk::ImageFileWriter< ImageType  >  WriterType;
typedef itk::AffineTransform< double, Dimension >  AffineTransformType;
typedef itk::VersorRigid3DTransform< double >  VersorRigid3DTransformType;

const unsigned int JointHistogramDimension = 2;
typedef itk::Image< PixelType,  JointHistogramDimension > HistogramImageType;

void updateJointHistogram(
    int refIndex,    // Moving image
    int fixedIndex,  // Fixed image
    int nbins,       // Size of quantization of histogram
    float weight,    // Amount to update histogram bin
    int threadIndex); // The thread doing the histogram update

double alignmentEvaluation(
  ImageType::Pointer inputRefVolImage, // The volume being resampled.
  const std::vector<ImageType::Pointer> &inputSliceImages,  // The fixed slices.
  VersorRigid3DTransformType::Pointer vrTransform,
  int nbins // The number of bins to use for the histogram.
);

HistogramImageType::Pointer allocImage( float **ptrHistogram );

ImageType::Pointer loadImage( std::string inputFileName, 
    float shift, float scale, int nbins );

#endif


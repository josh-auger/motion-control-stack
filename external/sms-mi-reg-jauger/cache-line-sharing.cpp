#include <iostream>
#include <cstring>

#include "omp.h"

#include <itkImage.h>
#include <itkOrientImageFilter.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkAffineTransform.h>
#include <itkPoint.h>
#include <itkContinuousIndex.h>

#include <itkVersorRigid3DTransform.h>
#include <itkPowellOptimizer.h>

// static allocation of arrays causes a compile time limit on the maximum
// number of threads that can be used.
const int MAX_NUM_THREADS = 64;

// Parameters associated with resampling volume to match a slice.
const int ARRAY_SIZE = 4096;
// pad to avoid false sharing of cache line.
const unsigned int padBest = 64/sizeof(float);// Should avoid false sharing.
const unsigned int padWorst = 0; // Generates false sharing of cache lines.
const unsigned int pad = padBest;
// int n[ARRAY_SIZE];
float n[ARRAY_SIZE][pad];

// Parameters associated with computing mutual information
// We need one histogram for each thread as they can't be written to
// in parallel.
const int nbins = 64;
//const int nbins = 32;
float jointHistograms[MAX_NUM_THREADS*nbins*nbins][pad];
float marginalFixed[nbins][pad];
float marginalMoving[nbins][pad];
float sumHistogram[nbins*nbins][pad];


int main(int argc, char *argv[])
{

  std::cout << "pad size of 64 bytes is " << 
                 pad << " unsigned int." << std::endl;

	double t0 = 0.0;
	double t1 = 0.0;

        // dynamic allocation of maximum number of threads this computer
        // can handle, capped by static limit.
	int n_threads = omp_get_max_threads();
        if (n_threads > MAX_NUM_THREADS) n_threads = MAX_NUM_THREADS;

        double threadStartTime[n_threads];
        double threadEndTime[n_threads];

        memset(threadStartTime,0,n_threads*sizeof(double));
        memset(threadEndTime,0,n_threads*sizeof(double));
	
	std::cout << "Maximum " << n_threads << " threads ..." << std::endl;

        //initialize n
        memset(n,0,pad*ARRAY_SIZE*sizeof(n[0][0]));
  
  // initialize jointHistograms
  memset(jointHistograms,0,pad*MAX_NUM_THREADS*nbins*nbins*
                    sizeof(jointHistograms[0][0]));
  memset(marginalFixed,0, pad*nbins*sizeof(marginalFixed[0][0]));
  memset(marginalMoving,0, pad*nbins*sizeof(marginalMoving[0][0]));
  memset(sumHistogram,0, pad*nbins*nbins*sizeof(sumHistogram[0][0]));

  const unsigned int Dimension = 3;
  typedef float PixelType;
  typedef itk::Image< PixelType,  Dimension >   ImageType;
  typedef itk::ImageFileReader< ImageType  >  ReaderType;
  typedef itk::ImageFileWriter< ImageType  >  WriterType;
  typedef itk::AffineTransform< double, Dimension >  AffineTransformType;

  std::string *inputRefImageFile = new std::string("refvol.nrrd");
  std::string *inputSliceImageFile = new std::string("inputSlice.nrrd");
  std::string *inputTransformFile = new std::string("transform.tfm");

  std::string *outputSliceImageFile = new std::string("outputSlice.nrrd");

  ReaderType::Pointer readerRefVol = ReaderType::New();
  ReaderType::Pointer readerSlice = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  readerRefVol->SetFileName( inputRefImageFile->c_str() );
  readerSlice->SetFileName( inputSliceImageFile->c_str() );
  writer->SetFileName( outputSliceImageFile->c_str() );

  try {
    readerRefVol->Update();
  } catch ( itk::ExceptionObject & excp )
  {
    // Display error from reading the reference volume file.
    std::cerr << "Error while reading the reference volume " <<
                     (*inputRefImageFile->c_str()) << std::endl;
    std::cerr << excp << std::endl;
    std::cerr << "[FAILED]" << std::endl;
    return EXIT_FAILURE;
  }

  ImageType::Pointer inputRefVolImage = readerRefVol->GetOutput();

  try {
    readerSlice->Update();
  } catch ( itk::ExceptionObject & excp )
  {
    // Display error from reading the reference volume file.
    std::cerr << "Error while reading the slice " <<
                     (*inputSliceImageFile->c_str()) << std::endl;
    std::cerr << excp << std::endl;
    std::cerr << "[FAILED]" << std::endl;
    return EXIT_FAILURE;
  }
  ImageType::Pointer inputSliceImage = readerSlice->GetOutput();

// NEED TO ALLOCATE THE OUTPUT IMAGE SLICE
  ImageType::Pointer outputSliceImage = ImageType::New();
  outputSliceImage->SetRegions(inputSliceImage->GetRequestedRegion());
  outputSliceImage->CopyInformation(inputSliceImage);
  outputSliceImage->Allocate();
  writer->SetInput(outputSliceImage);

  itk::TransformFileReader::Pointer trsfreader;
  trsfreader = itk::TransformFileReader::New();
  trsfreader->SetFileName( inputTransformFile->c_str() );

  try {
    trsfreader->Update();
  } catch ( itk::ExceptionObject & excp )
  {
    // Display error from reading the transform file.
       std::cerr << "Error while reading the transform file " <<
         (*inputTransformFile) << std::endl;
       std::cerr << excp << std::endl;
       std::cerr << "[FAILED]" << std::endl;
      return EXIT_FAILURE;
  }

 // Now try to work out how many and what type of transforms were read.
 // We only want to get one transform.
    typedef itk::TransformFileReader::TransformListType * TransformListType;
    TransformListType transforms = trsfreader->GetTransformList();
    std::cout << "Number of transforms = " << transforms->size() << std::endl;
    itk::TransformFileReader::TransformListType::const_iterator it =
      transforms->begin();
    if (transforms->size() <= 0 || transforms->size() > 1) {
      std::cerr << "Read " << transforms->size() << " transforms but want 1." << std::endl;
      return EXIT_FAILURE;
    }

    AffineTransformType::Pointer affine_read = 0;
    AffineTransformType::Pointer affine_inverse = 0;

    if (!strcmp((*it)->GetNameOfClass(), "AffineTransform")) {
        affine_read = 
                    static_cast<AffineTransformType*>((*it).GetPointer());
        affine_read->Print(std::cout);
        affine_inverse = AffineTransformType::New();
      if(!affine_read->GetInverse(affine_inverse)) {
        std::cerr << "Input transform is not invertible." << std::endl;
      }
    }

    

    typedef itk::BSplineInterpolateImageFunction<
                       ImageType, double >  BSplineInterpolatorType;
    BSplineInterpolatorType::Pointer bsplineInterpolator =
                       BSplineInterpolatorType::New();

    typedef itk::LinearInterpolateImageFunction<
                       ImageType, double >  LinearInterpolatorType;
    LinearInterpolatorType::Pointer linearInterpolator =
                       LinearInterpolatorType::New();

   linearInterpolator->SetInputImage(inputRefVolImage);

    ImageType::RegionType region = inputSliceImage->GetLargestPossibleRegion();

    // Construct an array of the input slice image
    unsigned int numVoxels = region.GetSize()[0]*
                             region.GetSize()[1]*
                             region.GetSize()[2];

    std::cout << "numVoxels is " << numVoxels << std::endl;

    // Remember that slices are not adjacent in the SMS acquisition.
    // The physical position needs to be accounted for.

  // for each slice
  for (int slice = 0; slice < region.GetSize()[2]; slice++) {

    for (unsigned int threadcount = 1; 
                      threadcount <= omp_get_max_threads(); threadcount++)
    {
      n_threads = threadcount;

      int job_size = ARRAY_SIZE/n_threads;
      int half_job_size = job_size/2;
      int quarter_job_size = half_job_size/2;
      int threadid = 0;

      // Initialize the jointHistograms back to zero.
      memset(jointHistograms,0,pad*MAX_NUM_THREADS*nbins*nbins*
                    sizeof(jointHistograms[0][0]));
      // Initialize marginal distributions back to zero.
      memset(marginalFixed,0, pad*nbins*sizeof(marginalFixed[0][0]));
      memset(marginalMoving,0, pad*nbins*sizeof(marginalMoving[0][0]));
      memset(sumHistogram,0, pad*nbins*nbins*sizeof(sumHistogram[0][0]));

       //start parallelization
       t0 = omp_get_wtime();
                
      #pragma omp parallel num_threads(n_threads)
      {
        threadStartTime[omp_get_thread_num()] = omp_get_wtime();
//#pragma omp for schedule(dynamic, quarter_job_size) private(index) nowait
#pragma omp for schedule(static, job_size) nowait
        for (unsigned int i = 0;i < region.GetSize()[0]*region.GetSize()[1];
                  i++) {
          ImageType::IndexType index;
          ImageType::PointType pointInSlice;
          ImageType::PointType pointTransformed;
          itk::ContinuousIndex<double, 3> cindex;

          float currentSliceValue = 0.0;
          float newVolValue = 0.0;
          int fixedIndex = 0;
          int refIndex = 0;

          index[2] = slice;
          index[1] = floor(i / region.GetSize()[0]);
          index[0] = i % region.GetSize()[0]; 

          currentSliceValue = inputSliceImage->GetPixel(index);

          inputSliceImage->TransformIndexToPhysicalPoint(index, pointInSlice);
          pointTransformed = affine_read->TransformPoint( pointInSlice );
          inputRefVolImage->TransformPhysicalPointToContinuousIndex( 
                pointTransformed, cindex );
          // Now interpolate the inputRefVolImage at the cindex. 
	  newVolValue = linearInterpolator->EvaluateAtContinuousIndex( cindex );
          // Since we turn this off and store the data in array n, 
          // the output slice is not being updated.
          // outputSliceImage->SetPixel(index, newVolValue);
          // Store in the data array in order to check moving image histogram
          n[i][0] = newVolValue; // Not a histogram, an image.
          // Joint histogram index is :
          // thread number , slice pixel value, resampled pixel value.
//convert to histogram indexes,
//check limits,
//rescale input images.
          fixedIndex = static_cast<int>( currentSliceValue );
          if (fixedIndex < 0 || fixedIndex >= nbins) {
            if (fixedIndex < 0) fixedIndex = 0;
            if (fixedIndex >= nbins) fixedIndex = nbins - 1;
            
            /* Don't display this warning for now.
            std::cerr << "Data value out of range: " << currentSliceValue
               << ", " <<fixedIndex << std::endl << std::flush;
            */
          }
          refIndex = static_cast<int>( newVolValue );
          if (refIndex < 0 || refIndex >= nbins) {
            if (refIndex < 0) refIndex = 0;
            if (refIndex >= nbins) refIndex = nbins - 1;
            
            std::cerr << "Resampled data value out of range: " 
               << newVolValue << ", " << refIndex
               << std::endl << std::flush;
          }
          int histogramIndex = 0;
          histogramIndex = fixedIndex * nbins + refIndex;
          if ( (histogramIndex < 0) || 
               (histogramIndex >= nbins*nbins) ) {
            std::cerr << "histogramIndex is out of range." << std::endl <<
                    std::flush;
          }
          // Now adjust to the location this thread is storing data in:
          histogramIndex += (omp_get_thread_num())*nbins*nbins;
          jointHistograms[histogramIndex][0] += 1;
        }
        threadEndTime[omp_get_thread_num()] = omp_get_wtime();
      }
      t1 = omp_get_wtime();
         
      // The following code is all running serially.
      float hitsCount = 0.0;
      // Copy the per-thread histograms to a unified histogram
      for (unsigned int i = 0; i < n_threads; i++) {
        for (unsigned int j = 0; j < nbins; j++) {
          for (unsigned int k = 0; k < nbins; k++) {
            hitsCount += jointHistograms[i*nbins*nbins + j*nbins + k][0];
            sumHistogram[j*nbins + k][0] += 
                       jointHistograms[i*nbins*nbins + j*nbins + k][0];
          }
        }
      }
      // Now normalize the summed histogram
      for (unsigned int j = 0; j < nbins; j++) {
        for (unsigned int k = 0; k < nbins; k++) {
          sumHistogram[j*nbins + k][0] /= hitsCount;
        }
      }

      // Compute the marginal fixed image and moving image histogram
      for (unsigned int j = 0; j < nbins; j++) {
        for (unsigned int k = 0; k < nbins; k++) {
          marginalFixed[j][0] += sumHistogram[j*nbins + k][0];
          marginalMoving[k][0] += sumHistogram[j*nbins + k][0];
        }
      }

      // Now we can compute the mutual information
      double mi = 0.0;
      double jointPDF = 0.0;
      double fixedPDF = 0.0;
      double movingPDF = 0.0;
      for (unsigned int j = 0; j < nbins; j++) {
        for (unsigned int k = 0; k < nbins; k++) {
          jointPDF = sumHistogram[j*nbins + k][0];
          if (jointPDF == 0.0) jointPDF = (1.0/hitsCount);
          fixedPDF = marginalFixed[j][0];
          if (fixedPDF == 0.0) fixedPDF = (1.0/hitsCount);
          movingPDF = marginalMoving[k][0];
          if (movingPDF == 0.0) movingPDF = (1.0/hitsCount);
          mi += (jointPDF * log2(jointPDF/(fixedPDF * movingPDF)));
        }
      }
      std::cout << "Mutual Information : " << mi << std::endl;
      std::cout << std::endl;

      std::cout << "job size = " << job_size << " with nthreads " <<
       n_threads << ", spent " << (t1 - t0) << " seconds." << std::endl;

      for (unsigned int i = 0; i < n_threads; i++) {
        std::cout << "Thread id " << i << " elapsed run time " <<
            (threadEndTime[i] - threadStartTime[i]) << std::endl;
        int pixelCount = 0;
        for (unsigned int j = 0; j < nbins*nbins; j++) {
          pixelCount += jointHistograms[i*nbins*nbins + j][0];
        }
        std::cout << "Thread " << i << " histogram hits count : " << 
          pixelCount << std::endl;
      }
    }

    writer->Update();
  }

  exit(0);
}


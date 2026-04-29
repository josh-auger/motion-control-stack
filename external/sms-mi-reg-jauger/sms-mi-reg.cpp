#include <cmath>
#include <vnl/vnl_inverse.h>
#include <vnl/algo/vnl_determinant.h>
#include <itkOptimizerParameters.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMacro.h>
#include <itkTransformFileWriter.h>
#include <itkStatisticsImageFilter.h>
#include "itkImageDuplicator.h"
#include "itkShiftScaleImageFilter.h"
#include "itkClampImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkOtsuMultipleThresholdsImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <iomanip>
#include <iostream>
#include <vector>
#include "nlopt.hpp"

// quill for logging
#include "quill/Backend.h"
#include "quill/Frontend.h"
#include "quill/LogMacros.h"
#include "quill/Logger.h"
#include "quill/sinks/FileSink.h"
#include "quill/sinks/JsonSink.h"
#include "quill/sinks/ConsoleSink.h"


#include <argparse/argparse.hpp>
#include <cassert>

// Include the Eigen template library for math
#include ITK_EIGEN(Core)
#include ITK_EIGEN(Dense)
#include ITK_EIGEN(LU)

typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> MatrixType33d;


#include "sms-mi-reg.h"


void updateJointHistogram(
    int refIndex,    // Moving image
    int fixedIndex,  // Fixed image
    int nbins,       // Size of quantization of histogram
    float weight,    // Amount to update histogram bin
    int threadIndex) // The thread doing the histogram update
// The joint histogram is a global variable
{
  assert(fixedIndex >= 0);
  assert(fixedIndex < nbins);
  assert(refIndex >= 0);
  assert(refIndex < nbins);

  int histogramIndex = fixedIndex * nbins + refIndex;
  // Now adjust to the location this thread is storing data in:
  histogramIndex += threadIndex*MaxNumBins*MaxNumBins;
  jointHistograms[histogramIndex][0] += weight;
};

// Carry out a partial volume based update to the joint histogram.
// We compute 8 weights and 8 bin updates.
// INPUT: fixed signal intensity
// INPUT: continuous index of point to weight
// computed: 2^Dimension = 8 moving image signal intensities
// computed: 8 weights
// OUTPUT: jointHistogram is directly updated.
void weightedUpdateJointHistogram(
    ImageType::Pointer inputRefVolImage, // reference volume
    itk::ContinuousIndex<double, 3> cindex, // ContinuousIndex location
    int fixedIndex,  // Fixed image
    int nbins,       // Size of quantization of histogram
    int threadIndex) // The thread doing the histogram update
{
  typedef float InternalComputationType;
  typedef ImageType::IndexType IndexType;
  typedef ImageType::IndexValueType IndexValueType;

  // Compute base index = closet index below point
  // Compute distance from point to base index
  ImageType::IndexType baseIndex;
  constexpr int ImageDimension = ImageType::ImageDimension; // JDA: evaluate at compile time rather than runtime
  InternalComputationType    distance[ImageDimension];

  /**
   * Interpolated value is the weighted sum of each of the surrounding
   * neighbors. The weight for each neighbor is the fraction overlap
   * of the neighbor pixel with respect to a pixel centered on point.
   */
  constexpr unsigned int m_Neighbors = 1 << ImageDimension;  // JDA: 2^ImageDimension using bitwise shifting, at compile time

  if (!inputRefVolImage->GetLargestPossibleRegion().IsInside(cindex)) {
    return;  // JDA: Early exit to avoid unnecessary processing of a given pixel
  }

  for ( unsigned int dim = 0; dim < ImageDimension; ++dim ) {
    baseIndex[dim] = itk::Math::Floor< IndexValueType >( cindex[dim] );
    distance[dim] = cindex[dim] - static_cast< InternalComputationType >( baseIndex[dim] );
  }

  // JDA: Get image bounds to check neighbor
  IndexType minIndex = inputRefVolImage->GetLargestPossibleRegion().GetIndex();
  IndexType maxIndex = minIndex + inputRefVolImage->GetLargestPossibleRegion().GetSize();

  // The bit pattern of the counter is used to index through the 2^dim neighbors
  // At each one a distance and overlap weight is calculated
  for ( unsigned int counter = 0; counter < m_Neighbors; ++counter ) {
    InternalComputationType overlap = 1.0;    // fraction overlap
    unsigned int upper = counter;  // each bit indicates upper/lower neighbour
    IndexType neighIndex( baseIndex ); // Variable for index of the neighbor

    // get neighbor index and overlap fraction
    for ( unsigned int dim = 0; dim < ImageDimension; ++dim ) {
      if ( upper & 1 ) {
        ++(neighIndex[dim]);
        // Take care of the case where the pixel is just
        // in the outer upper boundary of the image grid.
        if (neighIndex[dim] < minIndex[dim] || neighIndex[dim] >= maxIndex[dim]) {
          overlap = 0.0; // Set update weight to zero for points outside.
          break; // JDA: early break if current neighbor is out of bounds
        }
        overlap *= distance[dim];
      } else {
        if (neighIndex[dim] < minIndex[dim] || neighIndex[dim] >= maxIndex[dim]) {
          overlap = 0.0; // This could only happen if cindex was outside of the image...?
          break; // JDA: early break if current neighbor is out of bounds
        }
        overlap *= 1.0 - distance[dim];
      }
      upper >>= 1; // Shift all the bits right to choose upper or lower for the next dimension to process.
    }

    // overlap now holds the weight for the particular neighbor index.
    // This will fail if neighIndex is outside the image buffer.
    if (overlap > 0.0) {
      int referenceImageSample = inputRefVolImage->GetPixel( neighIndex );
      updateJointHistogram(referenceImageSample, fixedIndex, nbins, overlap, threadIndex);
    }
  }
};

// JDA: Define a struct to hold trace of all evaluated transform parameters and MI values
struct RegistrationInfo {
    std::vector<double> evaluatedTransformParams;
    std::vector<double> fixedTransformParams;
    double mutualInformation;
    double initTime;
    double time01;
    double time02;
    double time03;
    double parallelTime;
    double serialTime;
};
std::vector<RegistrationInfo> registrationInfos;

// JDA: Define a struct to hold joint histogram values
struct HistogramData {
    double xindex;
    double yindex;
    double bincount;
};
std::vector<HistogramData> histogramDatas;


double alignmentEvaluation(
  ImageType::Pointer inputRefVolImage, // The volume being resampled.
  const std::vector<ImageType::Pointer> &inputSliceImages,  // The fixed slices.
  VersorRigid3DTransformType::Pointer vrTransform,
  int nbins // The number of bins to use for the histogram.
)
{
  double mi = 0.0;

  // true when we want a single pair update, and 
  // false when we want a neighborhood weighted update of the histogram.
  // We believe that the neighbor weighted updates are the better estimator
  // for mutual information.
  bool singlePairUpdate = true;

  if (inputSliceImages.size() < 1) {
    std::cerr << "Only " << inputSliceImages.size() 
              << " provided. Need more." << std::endl;
    return 0.0;
  }
  // We are going to use the input image slices to determine the size of
  // the target geometry. The reference volume is resampled to match the slices.
  ImageType::RegionType region =inputSliceImages[0]->GetLargestPossibleRegion();

  double t00 = 0.0;
//  double t01 = 0.0;
//  double t02 = 0.0;
//  double t03 = 0.0;
  double t0 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;

  t00 = omp_get_wtime();

  // dynamic allocation of maximum number of threads this computer
  // can handle, capped by static limit.
  int n_threads = 1;
#ifdef _OPENMP
  n_threads = omp_get_max_threads();
#endif
  if (n_threads > MAX_NUM_THREADS) n_threads = MAX_NUM_THREADS;

  double threadStartTime[n_threads];
  double threadEndTime[n_threads];

//  t01 = omp_get_wtime();

  memset(threadStartTime,0,n_threads*sizeof(double));
  memset(threadEndTime,0,n_threads*sizeof(double));

//  t02 = omp_get_wtime();
	
  // initialize jointHistograms, sumHistogram and marginals to zero counts.
  memset(jointHistograms,0,pad*n_threads*MaxNumBins*MaxNumBins*sizeof(jointHistograms[0][0]));
//  #pragma omp parallel for
//  for (int i = 0; i < n_threads; i++) {
//    memset(jointHistograms[i],0,pad*MaxNumBins*MaxNumBins*sizeof(jointHistograms[0][0]));
//  }
  memset(marginalFixed,0,pad*MaxNumBins*sizeof(marginalFixed[0][0]));
  memset(marginalMoving,0,pad*MaxNumBins*sizeof(marginalMoving[0][0]));
  memset(sumHistogram,0,pad*MaxNumBins*MaxNumBins*sizeof(sumHistogram[0][0]));

//  t03 = omp_get_wtime();

  typedef itk::LinearInterpolateImageFunction< ImageType, double >  LinearInterpolatorType;
  LinearInterpolatorType::Pointer linearInterpolator = LinearInterpolatorType::New();
  linearInterpolator->SetInputImage(inputRefVolImage);

  // Construct an array of the input slice image
  unsigned int numVoxels = region.GetSize()[0]*
                             region.GetSize()[1]*
                             region.GetSize()[2]*
                             inputSliceImages.size();

  int colMult = 1;
  int rowMult = colMult*region.GetSize()[0];
  int slabMult = rowMult*region.GetSize()[1];
  int sliceMult = slabMult*region.GetSize()[2];

  // Remember that slices are not adjacent in the SMS acquisition.
  // The physical position needs to be accounted for.

  int job_size = numVoxels/n_threads;
  int half_job_size = job_size/2;
  int quarter_job_size = half_job_size/2;
  int threadid = 0;

//  std::cout << "Evaluated transform params : " << vrTransform->GetParameters() << std::endl;

  //start parallelization
  t0 = omp_get_wtime();

  #pragma omp parallel num_threads(n_threads)
  {
    int threadid = omp_get_thread_num();
    double localStartTime = omp_get_wtime();

    #pragma omp for schedule(dynamic, quarter_job_size) nowait
    //#pragma omp for schedule(static, job_size)
    for (unsigned int i = 0;i < numVoxels; i++) {   // JDA: Only search the number of pixels constituting the newly acquired slice(s)
      ImageType::IndexType index;
      ImageType::PointType pointInSlice;
      ImageType::PointType pointTransformed;
      itk::ContinuousIndex<double, 3> cindex;

      index[0] = static_cast<int>( floor(i/1) ) % static_cast<int>(rowMult);
      index[1] = static_cast<int>( floor(i/rowMult) ) % static_cast<int>(region.GetSize()[1]);
      index[2] = static_cast<int>( floor(i/slabMult)) % region.GetSize()[2];
      int slice = static_cast<int>( floor(i/sliceMult)) % inputSliceImages.size();

      ImageType* currentSlice = inputSliceImages[slice];
      float currentSliceValue = currentSlice->GetPixel(index);  // JDA: grab a specific pixel from the newly acquired slice(s)

      currentSlice->TransformIndexToPhysicalPoint( index, pointInSlice );
      pointTransformed = vrTransform->TransformPoint( pointInSlice );   // JDA: identify the new spatial location of the pixel after transformation
      inputRefVolImage->TransformPhysicalPointToContinuousIndex( pointTransformed, cindex ); // JDA: identify the equivalent spatial location within the reference volume
      /* LOGGER For debugging:
      if (threadid == 0) {
        std::cout << "Point: " << pointInSlice << " - " << "Point Transformed: " << pointTransformed << std::endl;
        std::cout << "cindex " << cindex << " - cindex " << cindex[0] << ", " << cindex[1] << ", " << cindex[2] << std::endl;
      }
      */
      // The recent ITK implementations do bounds checking on the neighbors.
      // They do not do bounds checking on the center point.
      // Values outside the buffer will be set to zero, where they are
      // accumulated in histogram bins on the edge,
      // where they are not included in the
      // mutual information calculation.
      float newVolValue = 0.0;
      if (linearInterpolator->IsInsideBuffer(cindex)) {
        newVolValue = linearInterpolator->EvaluateAtContinuousIndex(cindex);    // JDA: interpolate the pixel value at the transformed spatial location within the reference volume
      }

      // These represent the location in the histogram:
      int fixedIndex = static_cast<int>( currentSliceValue );
      int refIndex = static_cast<int>( newVolValue );

      // JDA: compile the new reference pixel value and newly acquired slice pixel into the joint histogram
      if (singlePairUpdate) {
        updateJointHistogram(
          refIndex,     // Moving image
          fixedIndex,   // Fixed image
          nbins,        // Size of quantization of histogram
          1.0,          // Amount to update histogram bin
          threadid      // The thread doing the histogram update
        );
      } else { // Eight pair update of the histogram:
        weightedUpdateJointHistogram(
          inputRefVolImage, // reference volume
          cindex,       // ContinuousIndex location
          fixedIndex,   // Fixed image
          nbins,        // Size of quantization of histogram
          threadid      // The thread doing the histogram update
        );
      }
    }   // JDA: repeat this per-pixel process for all pixel indices in the newly acquired slice(s)
    // JDA: after completing for all slice indices, we are effectively comparing a transformed-and-resampled slice(s) from the reference volume with the newly acquired slice(s)
    double localEndTime = omp_get_wtime();
    // Store time results in a thread-local way
    threadStartTime[threadid] = localStartTime;
    threadEndTime[threadid] = localEndTime;
  }
  t1 = omp_get_wtime();

  // The following code is all running serially.

  // The histogram uses the first and last bin to collect data <= min value and >= max value.
  // These bins are excluded from the PDF estimation to calculate mutual information.
  int minBinIndex = 1;
  int maxBinIndex = nbins - 1;
//  int minBinIndex = 0;
//  int maxBinIndex = nbins;

  // When estimating a pdf from a sample, an optimal estimator may use
  // prior information. We are choosing to encode this prior in the form
  // of a Dirichlet distribution.
  // This alters the values in the joint histogram so that every bin
  // has a count of at least one, which is the mathematically correct 
  // implementation of the Dirichlet prior.
  if (useDirichletPrior) {
    for (unsigned int j = minBinIndex; j < maxBinIndex; j++) {
      for (unsigned int k = minBinIndex; k < maxBinIndex; k++) {
        sumHistogram[j*nbins + k][0] = 1;
      }
    }
  }

  // Copy the per-thread histograms to a unified joint histogram
  for (unsigned int i = 0; i < n_threads; i++) {
    for (unsigned int j = minBinIndex; j < maxBinIndex; j++) {
      for (unsigned int k = minBinIndex; k < maxBinIndex; k++) {
        sumHistogram[j*nbins + k][0] += jointHistograms[i*MaxNumBins*MaxNumBins+ j*nbins + k][0];
      }
    }
  }

  // JDA: The sumHistogram is the joint histogram composed from the per-thread jobs, before normalization and before handling empty bins

  float hitsCount = 0.0;
  float totalBins = 0.0;
  for (unsigned int j = minBinIndex; j < maxBinIndex; j++) {
    for (unsigned int k = minBinIndex; k < maxBinIndex; k++) {
      // JDA: Populate struct with the joint histogram data
      HistogramData histoData;
      histoData.xindex = j;
      histoData.yindex = k;
      histoData.bincount = sumHistogram[j * nbins + k][0];
      histogramDatas.push_back(histoData);

      // Cumulative sum of all the counts in the joint histogram
      hitsCount += sumHistogram[j*nbins + k][0];
      totalBins += 1.0;
    }
  }
//  std::cout << "Total hitsCount : " << hitsCount << std::endl;
//  std::cout << "totalBins used : " << totalBins << std::endl;

  // Now normalize the joint histogram into probability density (PDF)
  for (unsigned int j = minBinIndex; j < maxBinIndex; j++) {
    for (unsigned int k = minBinIndex; k < maxBinIndex; k++) {
      sumHistogram[j*nbins + k][0] /= hitsCount;
    }
  }

  // Extract the marginal fixed image and moving image PDF histograms by summing the joint PDF along each axis
  for (unsigned int j = minBinIndex; j < maxBinIndex; j++) {
    for (unsigned int k = minBinIndex; k < maxBinIndex; k++) {
      marginalFixed[j][0] += sumHistogram[j*nbins + k][0]; // cumulative sum along j
      marginalMoving[k][0] += sumHistogram[j*nbins + k][0]; // cumulative sum along k
    }
  }


  // Now we can compute the mutual information.
  // See Maes et al. 2003 Equation 8.
  //
  mi = 0.0;
//  double mientropy = 0.0;
  double jointPDF = 0.0;
  double fixedPDF = 0.0;
  double movingPDF = 0.0;
//  double jointentropy = 0.0;
//  double fixedentropy = 0.0;
//  double movingentropy = 0.0;
  for (unsigned int j = minBinIndex; j < maxBinIndex; j++) {
    for (unsigned int k = minBinIndex; k < maxBinIndex; k++) {
      jointPDF = sumHistogram[j*nbins + k][0];
      fixedPDF = marginalFixed[j][0];
      movingPDF = marginalMoving[k][0];

      // Only non-zero PDF values contribute to MI
      if (jointPDF > 0.0) {
        // Calculate joint entropy
//        jointentropy += -(jointPDF * log2(jointPDF));

        // Calculate MI from PDFs
        mi += (jointPDF * log2(jointPDF/(fixedPDF * movingPDF)));
      }
    }
  }

//  // JDA: Cycle through indices of each marginal PDF separately to avoid multi-counting the same value towards entropy
//  for (unsigned int j = minBinIndex; j < maxBinIndex; j++) {
//    fixedPDF = marginalFixed[j][0];
//    if (fixedPDF > 0.0) {
//      fixedentropy += -(fixedPDF * log2(fixedPDF));
//    }
//
//    movingPDF = marginalMoving[j][0];
//    if (movingPDF > 0.0) {
//      movingentropy += -(movingPDF * log2(movingPDF));
//    }
//  }
//
//  // JDA: calculate mutual information from entropy measures
//  mientropy = fixedentropy + movingentropy - jointentropy;
//
//  std::cout << "Fixed entropy : " << fixedentropy << std::endl;
//  std::cout << "Moving entropy : " << movingentropy << std::endl;
//  std::cout << "Joint entropy : " << jointentropy << std::endl;
//  std::cout << "MI from entropy : " << mientropy << std::endl;

  // End of the serial part of the code.
  t2 = omp_get_wtime();

  /* LOGGER: TURN THIS ON TO SEE THREADING PERFORMANCE REPORTED.
  std::cout << "Initial pre-processing time (sec) : " << (t0 - t00) << std::endl;
  std::cout << "Parallel time (sec) : " << (t1 - t0) << std::endl;
  std::cout << "Serial time (sec) : " << (t2 - t1) << std::endl;
  std::cout << "Mutual Information is : " << mi << std::endl;
  */

  /* // JDA: Populate struct with most recent registration parameters and MI value */
  RegistrationInfo regInfo;
  regInfo.evaluatedTransformParams = std::vector<double>(vrTransform->GetParameters().begin(), vrTransform->GetParameters().end());
  regInfo.fixedTransformParams = std::vector<double>(vrTransform->GetFixedParameters().begin(), vrTransform->GetFixedParameters().end());
  regInfo.mutualInformation = mi;
  regInfo.initTime = (t0 - t00);
//  regInfo.time01 = (t01 - t00);
//  regInfo.time02 = (t02 - t01);
//  regInfo.time03 = (t03 - t02);
  regInfo.parallelTime = (t1 - t0);
  regInfo.serialTime = (t2 - t1);
  registrationInfos.push_back(regInfo);

  return mi;

}; // End of metric evaluation

class MICostFunction : public itk::SingleValuedCostFunction
{
  public:
    typedef MICostFunction Self;
    typedef itk::SingleValuedCostFunction Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;
    itkNewMacro( Self );
    // itkTypeMacro( MICostFunction, SingleValuedCostFunction );

    static const unsigned int Dimension = 3;
    typedef float PixelType;
    typedef itk::Image< PixelType,  Dimension >   ImageType;
    int nbins = MaxNumBins;

    ImageType::Pointer inputRefVolImage = nullptr; // The volume being resampled.
    // The pointer to slices to match to.
    std::vector<ImageType::Pointer> inputSliceImages; 
    VersorRigid3DTransformType::Pointer vrTransform = nullptr;

    typedef Superclass::ParametersType ParametersType;
    typedef Superclass::DerivativeType DerivativeType;
    typedef Superclass::MeasureType MeasureType;

    MICostFunction() 
    {
      vrTransform = VersorRigid3DTransformType::New();
      // Initialize the input array
      inputSliceImages.clear();
    }

    virtual unsigned int GetNumberOfParameters() const override {
      // There are three fixed parameters that are not optimized.
      // These are the x,y,z center of rotation coordinates.
      // There are then three translation and three rotation parameters.
      return Dimension*2;
    }

    void SetFixedParameters(
              VersorRigid3DTransformType::FixedParametersType fixedParameters)
    {
      if (vrTransform == nullptr) {
        std::cerr << "vrTransform should not be null." << std::endl;
        exit(1);
      }
      vrTransform->SetFixedParameters( fixedParameters );
    }

    void SetInputReferenceVolume(ImageType::Pointer input)
    {
      inputRefVolImage = input;
    }

    void SetInputSliceImage(ImageType::Pointer slice)
    {
      inputSliceImages.push_back(slice);
    }

    int GetNumberOfSlices()
    {
      return inputSliceImages.size();
    }

    void SetNumberOfBins( int numBins ) {
      nbins = numBins;
    };

    int GetNumberOfBins( ) {
      return nbins;
    };

    // Need some way to set the center of rotation.

    void GetDerivative( const ParametersType &, 
             DerivativeType & ) const ITK_OVERRIDE
    {
    }

    // This is the three translation and three rotation parameters only.
    void const ConvertParametersToTransform( const 
                    VersorRigid3DTransformType::ParametersType &parms )
    {
      if (vrTransform == nullptr) {
        std::cerr << "vrTransform should not be null." << std::endl;
        exit(1);
      }
      vrTransform->SetParameters( parms );
    }

    virtual MeasureType GetValue( const ParametersType &parameters )
        const ITK_OVERRIDE
    {
      // PRECONDITIONS: Center of rotation has been set.
      //   The center of rotation is set by the fixed parameters.
      //                input volume and slices have been set.

      double mi = 0.0;

      if (vrTransform == nullptr) {
        std::cerr << "vrTransform should not be null." << std::endl;
        exit(1);
      }
      // 1. convert the parameters to a transform
      // 2. evaluate the alignment implied by the transform
      vrTransform->SetParameters( parameters );

      mi = alignmentEvaluation(
             inputRefVolImage, // The volume being resampled.
             inputSliceImages,  // The fixed slices being targeted.
             vrTransform,
             nbins
            );

      return mi;
    }
};

struct my_func_data_struct
{
  MICostFunction::Pointer costFunctionPointer = nullptr;
};

// This requires all of the data structures have already been initialized, including the
// FixedParameters of the vrTransform.
//
// double my_cost_function(unsigned int n, const double *x, double *grad, void *my_func_data)
double my_cost_function(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
  double mi = 0.0;
  my_func_data_struct *localPtr = (my_func_data_struct *)my_func_data;
  MICostFunction::Pointer costFunctionPointer = localPtr->costFunctionPointer;
  itk::OptimizerParameters<double> transformParameters(6);

  if (x.size() != 6) {
    std::cerr << "Wrong number of parameters." << std::endl;
    return 2; // Wrong number of parameters.
  }

  if (!grad.empty()) {
    std::cerr << "grad is not yet implemented." << std::endl;
    std::cerr << "set grad to be empty." << std::endl;
    exit(3);
  }

  if (!my_func_data) {
    std::cerr << "Wrong call to my_cost_function." << std::endl;
    exit(1);
  }

  for (unsigned int i = 0; i < x.size(); i++) {
    transformParameters[i] = x[i];
  }

  mi = costFunctionPointer->GetValue( transformParameters );

  return mi;
}


bool writeImage( ImageType::Pointer inputImage, std::string outputFileName )
{
  using WriterType = itk::ImageFileWriter<ImageType>;
  auto writer = WriterType::New();
  writer->SetFileName( outputFileName );
  writer->SetInput( inputImage );

  try
  {
    writer->Update();
    // ITKv5 ? itk::WriteImage(clonedImage, outputFileName );
  }
  catch (const itk::ExceptionObject & excp)
  {
    std::cerr << "Error: " << excp.GetDescription() << std::endl;
    return false;
  }

  writer = 0;

  return true;
}

bool writeHistogramImage( HistogramImageType::Pointer inputImage, 
                          std::string outputFileName )
{
  using WriterType = itk::ImageFileWriter<HistogramImageType>;
  auto writer = WriterType::New();
  writer->SetFileName( outputFileName );
  writer->SetInput( inputImage );

  try
  {
    writer->Update();
    // ITKv5 ? itk::WriteImage(clonedImage, "otsuThresholded.nrrd" );
  }
  catch (const itk::ExceptionObject & excp)
  {
    std::cerr << "Error: " << excp.GetDescription() << std::endl;
    return false;
  }

  writer = 0;

  return true;
}

bool printImageStatistics( ImageType::Pointer inputImage )
{
  using StatisticsImageFilterType = itk::StatisticsImageFilter<ImageType>;
  auto statFilter = StatisticsImageFilterType::New();
  statFilter->SetInput( inputImage );
  statFilter->Update();
 
  std::cout << "Statistics for image:" << std::endl;
  std::cout << "Mean: " << statFilter->GetMean() << std::endl;
  std::cout << "Std.: " << statFilter->GetSigma() << std::endl;
  std::cout << "Min: " << statFilter->GetMinimum() << std::endl;
  std::cout << "Max: " << statFilter->GetMaximum() << std::endl;

  statFilter = 0;

  return true;
}

bool printHistogramImageStatistics( HistogramImageType::Pointer inputImage )
{
  using StatisticsImageFilterType = 
        itk::StatisticsImageFilter<HistogramImageType>;
  auto statFilter = StatisticsImageFilterType::New();
  statFilter->SetInput( inputImage );
  statFilter->Update();
 
  std::cout << "Statistics for histogram image:" << std::endl;
  std::cout << "Mean: " << statFilter->GetMean() << std::endl;
  std::cout << "Std.Dev: " << statFilter->GetSigma() << std::endl;
  std::cout << "Min: " << statFilter->GetMinimum() << std::endl;
  std::cout << "Max: " << statFilter->GetMaximum() << std::endl;

  statFilter = 0;

  return true;
}

// The purpose of this function is to allocate an image to be able to
// represent the 2D histograms we calculate as an ITK image.
HistogramImageType::Pointer allocImage( 
    float (*histogramPtr)[MaxNumBins*MaxNumBins][pad] )
{
  const unsigned int Dimension = 2;
  typedef float PixelType;
  typedef itk::Image< PixelType,  Dimension > HistogramImageType;

  auto image = HistogramImageType::New();

  HistogramImageType::IndexType start;
  start[0] = 0; // first index on X
  start[1] = 0; // first index on Y

  HistogramImageType::SizeType size;
  size[0] = MaxNumBins; // size along X
  size[1] = MaxNumBins; // size along Y

  HistogramImageType::RegionType region;
  region.SetSize(size);
  region.SetIndex(start);

  image->SetRegions(region);
  image->Allocate();
  image->FillBuffer(itk::NumericTraits<PixelType>::Zero);

  // Now copy the data from the histogram to the image
  HistogramImageType::IndexType idx;
  for (unsigned int j = 0; j < MaxNumBins; j++) {
    for (unsigned int k = 0; k < MaxNumBins; k++) {
      idx[0] = k;
      idx[1] = j;
      float voxelValue = (*histogramPtr)[j*MaxNumBins+ k][0];
      image->SetPixel(idx, voxelValue);
    }
  }

  return image;
}

ImageType::Pointer loadImage( std::string inputFileName, float shift, float scale , int nbins)
{
    /* INPUT: 
     * string: slice file name
     * float: shift, scale.
     * OUTPUT: slice image, rescaled
     */
  const unsigned int Dimension = 3;
  typedef float PixelType;
  typedef itk::Image< PixelType, Dimension >   ImageType;
  typedef itk::ImageFileReader< ImageType  >  ReaderType;
  typedef itk::ImageFileWriter< ImageType  >  WriterType;

  using StatisticsImageFilterType = itk::StatisticsImageFilter<ImageType>;
  using shiftScaleFilterType = itk::ShiftScaleImageFilter<ImageType, ImageType>;
  using ClampFilterType = itk::ClampImageFilter<ImageType, ImageType>;

  ReaderType::Pointer readerSlice = ReaderType::New();
  readerSlice->SetFileName( inputFileName.c_str() );

  try {
    readerSlice->Update();
  } catch ( itk::ExceptionObject & excp )
  {
    std::cerr << "Error while reading the input data from " <<
                     (inputFileName.c_str()) << std::endl;
    std::cerr << excp << std::endl;
    std::cerr << "[FAILED]" << std::endl;
    return nullptr;
  }

  auto filterSlice = shiftScaleFilterType::New();
  filterSlice->SetInput( readerSlice->GetOutput() );
  filterSlice->SetShift( shift );
  filterSlice->SetScale( scale );
  filterSlice->Update();
  std::cout << "Underflow: " <<  filterSlice->GetUnderflowCount() << std::endl;
  std::cout << "Overflow : " <<  filterSlice->GetOverflowCount() << std::endl;

  // We are also going to set the min and max range to ensure we don't exceed
  // the histogram bin range. The purpose of this is to accelerate the
  // later binning used in pdf estimation.
  //   The clamp filter sets values less than or equal to the lower bound
  // to the lower bound. 
  //    It sets values greater than or equal to the upper bound to 
  // the upper bound.
  //    The number of bins should be set based on whether or not the lower
  // and upper bins will be excluded.
  
  using ClampFilterType = itk::ClampImageFilter<ImageType, ImageType>;
  auto clampFilterSlice = ClampFilterType::New();
  clampFilterSlice->SetInput( filterSlice->GetOutput() );
  clampFilterSlice->SetBounds(0.0f, nbins - 1.0f);
  // clampFilterSlice->SetBounds(0.0f, nbins );
  clampFilterSlice->Update();

  ImageType::Pointer inputSliceImage = clampFilterSlice->GetOutput();

  // Check the remapped slice signal intensities range [0, nbins - 1]
  StatisticsImageFilterType::Pointer statisticsImageFilterSlice = 
      StatisticsImageFilterType::New();
  statisticsImageFilterSlice->SetInput( inputSliceImage );
  statisticsImageFilterSlice->Update();
 
  std::cout << "Statistics for slice : " + inputFileName  << std::endl;
  std::cout << "Mean: " << statisticsImageFilterSlice->GetMean() << std::endl;
  std::cout << "Std.: " << statisticsImageFilterSlice->GetSigma() << std::endl;
  std::cout << "Min: "<< statisticsImageFilterSlice->GetMinimum() << std::endl;
  std::cout << "Max: "<< statisticsImageFilterSlice->GetMaximum() << std::endl;

  printImageStatistics( inputSliceImage );

  return inputSliceImage;

}


// function to map input optimizer name to NLOPT algorithm
nlopt::algorithm getOptimizerAlgorithm(const std::string& name) {
    // Local derivative-free optimization schemes
    if (name == "LN_COBYLA") return nlopt::LN_COBYLA;
    if (name == "LN_BOBYQA") return nlopt::LN_BOBYQA;
    if (name == "LN_NELDERMEAD") return nlopt::LN_NELDERMEAD;
    if (name == "LN_SBPLX")  return nlopt::LN_SBPLX;
    // Add more here as needed...
    throw std::invalid_argument("Unknown optimizer: " + name);
}



int main(int argc, char *argv[])
{

  const unsigned int Dimension = 3;
  typedef float PixelType;
  typedef itk::Image< PixelType,  Dimension >   ImageType;
  typedef itk::ImageFileReader< ImageType  >  ReaderType;
  typedef itk::ImageFileWriter< ImageType  >  WriterType;
  typedef itk::AffineTransform< double, Dimension >  AffineTransformType;
  typedef itk::VersorRigid3DTransform< double >  VersorRigid3DTransformType;

  // Start the backend thread
  quill::BackendOptions backend_options;
  quill::Backend::start(backend_options);

  // The maximum number of bins allowed by the memory allocation is MaxNumBins
  // JDA: Here, set nbins to the actual number of bins used to populate the joint histogram
  int nbins = 64;
  std::cout << "nbins used for joint histogram : " << nbins << std::endl;

  if (nbins > MaxNumBins) {
    std::cout << "More bins requested than allocated." << std::endl;
    std::exit(1);
  }

  // Frontend
//  auto file_sink = quill::Frontend::create_or_get_sink<quill::FileSink>(
//    "example_file_logging.log",
//    []()
//    {
//      quill::FileSinkConfig cfg;
//      cfg.set_open_mode('w');
//      cfg.set_filename_append_option(quill::FilenameAppendOption::StartDateTime);
//      return cfg;
//    }(),
//    quill::FileEventNotifier{});
//
//  quill::Logger* logger = quill::Frontend::create_or_get_logger(
//    "root", std::move(file_sink),
//    quill::PatternFormatterOptions{"%(time) [%(thread_id)] %(short_source_location:<28) "
//                                          "LOG_%(log_level:<9) %(logger:<12) %(message)",
//                                   "%H:%M:%S.%Qns", quill::Timezone::GmtTime});
//
//  // set the log level of the logger to debug (default is info)
//  logger->set_log_level(quill::LogLevel::Debug);
//
//  LOG_INFO(logger, "Starting sms-mi-reg!");

  argparse::ArgumentParser program("sms-mi-reg");
  program.add_argument("referenceVolume")
      .help("The volume that is moving to be aligned to the slices.");
  program.add_argument("inputTransform")
      .help("The transform to initialize the alignment.");
  program.add_argument("outputTransformLabel")
      .help("Name phrase used in the construction of the output transform file name.");
  program.add_argument("inputSlices")
      .help("The list of file names of the fixed target slices.")
      .nargs(argparse::nargs_pattern::at_least_one); // "+". This accepts one or more file name arguments.
  program.add_argument("--optimizer")
    .help("Choice of optimizer (LN_COBYLA, LN_BOBYQA, LN_NELDERMEAD, LN_SBPLX). Default is LN_BOBYQA.")
    .default_value(std::string("LN_BOBYQA"));
  program.add_argument("--maxiter")
    .help("Maximum number of optimizer iterations. Default is 1000.")
    .default_value(1000)
    .scan<'i', int>();   // ensure it parses as int

  try {
    program.parse_args(argc, argv);
  }
  catch (const std::runtime_error &err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    std::exit(1);
  }
 
  //adding new arguments
  std::string referenceVolume = program.get<std::string>("referenceVolume");
//  LOG_INFO(logger, "Reference volume is {}." , referenceVolume);
  std::string transformfile = program.get<std::string>("inputTransform");
  std::string outputTransformLabel = program.get<std::string>("outputTransformLabel");
  std::vector<std::string> inputSliceFileNames = program.get<std::vector<std::string>>("inputSlices");
  std::string optimizerName = program.get<std::string>("--optimizer");
  int maxeval_value = program.get<int>("--maxiter");


  std::cout << "Processing for " << inputSliceFileNames.size() << " input slices." << std::endl;
//  LOG_INFO(logger, "Number of input file names : {}", inputSliceFileNames.size());
  for (auto i : inputSliceFileNames) {
    std::cout << "input slice file names : " << std::endl << i << std::endl;
  }

  std::string *inputTransformFile = new std::string( transformfile );
  std::string *inputRefImageFile = new std::string( referenceVolume );

  ReaderType::Pointer readerRefVol = ReaderType::New();
  readerRefVol->SetFileName( inputRefImageFile->c_str() );

  std::string * outputSliceImageFile = new std::string("scaled-outputSlice.nrrd");


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

  using StatisticsImageFilterType = itk::StatisticsImageFilter<ImageType>;
  StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New();
  statisticsImageFilter->SetInput( readerRefVol->GetOutput() );
  statisticsImageFilter->Update();
  std::cout << "Statistics for native reference volume:" << std::endl;
  printImageStatistics( readerRefVol->GetOutput() );

// This example program uses a heuristic to compute a mapping from the
// signal intensity range of the input to the histogram bins.
//  
// It uses a rescale that maps the input signal intensities to the desirable 
// bin range [0, nbins - 1], save the mapping function, apply it to all the 
// slices.

// 1. Otsu threshold to find background and foreground.
// 2. Measure the mean and variance of the foreground voxels.
// 3. Determine a linear remapping based on the foreground intensities.
// 4. Ensuring valid range of [0, nbins - 1] is preserved
// 5. Apply scaling to the slices.

  /* The Otsu threshold is a binary threshold and appears inadequate for fmri,
   * so instead we use a multiple threshold Otsu filter.
   *
   * using OtsuFilterType = itk::OtsuThresholdImageFilter<ImageType, ImageType>;
  */

//  LOG_INFO(logger, "Number of bins : {}", nbins);

  using MOtsuFilterType = itk::OtsuMultipleThresholdsImageFilter<ImageType, ImageType>;
  auto motsuFilter = MOtsuFilterType::New();
  motsuFilter->SetInput( readerRefVol->GetOutput() );
  motsuFilter->SetNumberOfHistogramBins( nbins );
  motsuFilter->SetNumberOfThresholds( 2 ); // 2 thresholds is suitable for fmri.
  motsuFilter->SetLabelOffset( 0 );
  motsuFilter->Update();

  MOtsuFilterType::ThresholdVectorType thresholds = motsuFilter->GetThresholds();

  /*
  // Now make a duplicate of the thresholded image to hold and save.
  using DuplicatorType = itk::ImageDuplicator<ImageType>;
  auto duplicator = DuplicatorType::New();
  duplicator->SetInputImage( motsuFilter->GetOutput() );
  duplicator->Update();
  ImageType::Pointer clonedImage = duplicator->GetOutput();

  // Writing out the otsu thresholded image.
  std::cout << "Writing out a copy of the otsu thresholded image " <<
    " otsuThresholded.nii " << std::endl;

  // Already defined: using WriterType = itk::ImageFileWriter<ImageType>;
  auto writerOtsu = WriterType::New();
  writerOtsu->SetFileName( "otsuThresholded.nii" );
  writerOtsu->SetInput( clonedImage );

  try
  {
    writerOtsu->Update();
    // ITKv5 ? itk::WriteImage(clonedImage, "otsuThresholded.nrrd" );
  }
  catch (const itk::ExceptionObject & excp)
  {
    std::cerr << "Error: " << excp.GetDescription() << std::endl;
    return EXIT_FAILURE;
  }
  */


  // Iterate over the input image, checking the otsu threshold image, 
  // accumulating the min, max, mean and variance.
  itk::ImageRegionIteratorWithIndex<ImageType> imageIterator(
      readerRefVol->GetOutput(), 
      readerRefVol->GetOutput()->GetLargestPossibleRegion() 
      );

  double voxelValue = 0.0;
  double mean = 0.0;
  double variance = 0.0;
  double oldmean = 0.0;
  double oldvariance = 0.0;
  double stddev = 0.0;
  unsigned long int voxelCount = 0;

  // Modified to use Welford's algorithm for variance
  imageIterator.GoToBegin();
  double minValue = std::numeric_limits<double>::max();
  double maxValue = std::numeric_limits<double>::lowest();
  double M2 = 0.0;
  while (!imageIterator.IsAtEnd()) {
    double delta = 0.0;
    double delta2 = 0.0;
    ImageType::IndexType idx = imageIterator.GetIndex();
    float otsuVoxel = motsuFilter->GetOutput()->GetPixel(idx);
    if (otsuVoxel > 0) {
      voxelCount += 1;
      voxelValue = imageIterator.Get();
      if (voxelValue < minValue) minValue = voxelValue;
      if (voxelValue > maxValue) maxValue = voxelValue;
      delta = voxelValue - mean;
      mean += delta / (double)voxelCount;
      delta2 = voxelValue - mean;
      M2 += delta * delta2;
      oldmean += voxelValue;
      oldvariance += voxelValue*voxelValue;
    }
    ++imageIterator;
  }

  if (voxelCount == 0) {
    std::cerr << "Otsu filter found no foreground voxels." << std::endl;
    std::cerr << "Exiting with a fatal error." << std::endl;
    exit(2);
  }
  oldmean /= voxelCount;
  oldvariance = (oldvariance/voxelCount - oldmean*oldmean);
  stddev = std::sqrt(oldvariance);
  std::cout << "Mean: " << oldmean << ", Variance: " << oldvariance 
        << " StdDev: " << stddev << std::endl << std::flush;

  variance = M2 / (voxelCount);
  stddev = std::sqrt(variance);
  std::cout << "Welford Mean: " << mean << ", Variance: " << (variance)
        << " StdDev: " << stddev << std::endl << std::flush;
  std::cout << "Otsu filtered min and max values:" << std::endl;
  std::cout << "  minValue : " <<  minValue << std::endl;
  std::cout << "  maxValue : " <<  maxValue << std::endl;

  std::cout << "Labelled regions have been reported on." << std::endl 
    << std::flush;

  using shiftScaleFilterType = itk::ShiftScaleImageFilter<ImageType, ImageType>;
  auto filter = shiftScaleFilterType::New();
  // outputPixel = (inputPixel + Shift) x Scale
  // Desired range is [0, nbins - 1]
  // mean + 3 std dev -> (nbins - 1)
  // mean - 3 std dev -> 0
  // Outside the range of (mean +/- 3 std dev), map to 0 and nbins - 1.
//  if (minValue < (mean - 3 * stddev)) minValue = (mean - 3 * stddev);
//  if (maxValue > (mean + 3 * stddev)) maxValue = (mean + 3 * stddev);
  minValue = 0;

  std::cout << "Trimmed min and max values:" << std::endl;
  std::cout << "  minValue : " <<  minValue << std::endl;
  std::cout << "  maxValue : " <<  maxValue << std::endl;
  float range = maxValue - minValue;
  std::cout << "range : " <<  range << std::endl;

  float scale = (nbins - 3) / range;
  float shift = (1 / scale) - minValue; // (-1.0) * minValue;
  std::cout << "Shift : " <<  shift << std::endl;
  std::cout << "Scale : " <<  scale << std::endl;
  std::cout << std::flush;

  filter->SetInput( readerRefVol->GetOutput() );
  filter->SetShift( shift );
  filter->SetScale( scale );
  filter->Update();
  std::cout << "Underflow count of voxels below trimmed min: " <<  filter->GetUnderflowCount() << std::endl;
  std::cout << "Overflow count of voxels above trimmed max: " <<  filter->GetOverflowCount() << std::endl;

  // We are also going to set the min and max range to ensure we don't exceed
  // the histogram bin range.
  using ClampFilterType = itk::ClampImageFilter<ImageType, ImageType>;
  auto clampFilter = ClampFilterType::New();
  clampFilter->SetInput( filter->GetOutput() );
  // Clamp the range so that voxels don't map outside the histogram.
  clampFilter->SetBounds(0.0f, nbins - 1.0f);
  // An option: add an additional bin of range in order to collect and count the out of range pixels.
  // clampFilter->SetBounds(0.0f, nbins);
  clampFilter->Update();

  ImageType::Pointer inputRefVolImage = clampFilter->GetOutput();

  // Check the remapped signal intensities range [0, nbins - 1]
//  using StatisticsImageFilterType = itk::StatisticsImageFilter<ImageType>;
//  StatisticsImageFilterType::Pointer statisticsImageFilter =
//      StatisticsImageFilterType::New();
  statisticsImageFilter->SetInput( inputRefVolImage );
  statisticsImageFilter->Update();
 
  /* Not displaying image statistics:
  std::cout << "Statistics for whole rescaled reference volume:" << std::endl;
  printImageStatistics( inputRefVolImage );
  */

  /*
  // Now make a duplicate of the scaled reference image to hold and save.
  auto duplicator2 = DuplicatorType::New();
  duplicator2->SetInputImage( inputRefVolImage );
  duplicator2->Update();
  ImageType::Pointer clonedImage2 = duplicator2->GetOutput();

  // Writing out the rescaledrefvol:
  std::cout << "Writing out a copy of the rescaled ref image" <<
     " as rescaledrefvol.nii"  << std::endl;

  auto writerDup2 = WriterType::New();
  writerDup2->SetFileName( "rescaledrefvol.nii" );
  writerDup2->SetInput( clonedImage2 );

  try
  {
    writerDup2->Update();
  }
  catch (const itk::ExceptionObject & excp)
  {
    std::cerr << "Error: " << excp.GetDescription() << std::endl;
    return EXIT_FAILURE;
  }
  */

  std::vector< ImageType::Pointer > inputSliceImages;
  inputSliceImages.resize( inputSliceFileNames.size() );
  for (unsigned int i = 0; i < inputSliceImages.size(); i++) {
    inputSliceImages[i] = loadImage( inputSliceFileNames[i], shift, scale,
                                nbins );
  }

  // This section of code processes the input transformation.
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
  VersorRigid3DTransformType::Pointer vrTransform =
    VersorRigid3DTransformType::New(); // This is the transform to be read in.
  VersorRigid3DTransformType::Pointer trsf_read =
    VersorRigid3DTransformType::New(); // This is the transform to be read in.
  AffineTransformType::Pointer affine_read =
    AffineTransformType::New(); // This is the transform that may be read in.

  VersorRigid3DTransformType::Pointer finalVRTransform =
    VersorRigid3DTransformType::New(); // This is transform found.
  VersorRigid3DTransformType::ParametersType parmArray( finalVRTransform->GetNumberOfParameters() );

  VersorRigid3DTransformType::Pointer identityVRTransform =
    VersorRigid3DTransformType::New(); // This is to hold an identity transform.


  typedef itk::TransformFileReader::TransformListType * TransformListType;
  TransformListType transforms = trsfreader->GetModifiableTransformList();
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
    // std::cout << "Just loaded a versor rigid 3d transform." << std::endl;
    // trsf_read->Print(std::cout);

    // This is the recommended procedure for initializing when the
    // transforms have different mappings of actual parameter to index
    // in the array of parameters.
    vrTransform->SetCenter( trsf_read->GetCenter() );
    vrTransform->SetTranslation( trsf_read->GetTranslation() );
        // API BUG: itk::VersorRigid3DTransform<double>' has no member
        //   named 'GetRotation'
        // API BUG: ScaleSkewVersor3DTransform has no SetVersor()
    vrTransform->SetRotation( trsf_read->GetVersor() );
    std::cout << "Input parameters : " << vrTransform->GetParameters() << std::endl;
  } else if (!strcmp((*trsfit)->GetNameOfClass(),"AffineTransform")) {
        affine_read = static_cast<AffineTransformType*>((*trsfit).GetPointer());
        // std::cout << "Just loaded an affine transform." << std::endl;
        // affine_read->Print(std::cout);
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

        // Using Eigen3 for the matrix operations
        // ITK-v5.3.0/Modules/ThirdParty/Eigen3

        // Get the 3x3 set of parameters.
        // There is also a set of fixed parameters that represent the origin of the transformation.

        AffineTransformType::ParametersType parms = affine_read->GetParameters();

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
        MatrixType33d PQi = MatrixType33d::Identity();
        MatrixType33d PQit = MatrixType33d::Identity();

        const unsigned int maximumIterations = 100;
        const double eps = 1e-11;

        // Calculate part of the polar decomposition of the matrix M.
        // Compute the matrix NQ as the best orthogonal approximation to the input matrix M.
        for(unsigned int ni = 0; ni < maximumIterations; ni++ )
        {
          // Average current Qi with its inverse transpose
          PQi = PQ.inverse();
          PQit = PQi.transpose();
          NQ = ( PQ + PQit) / 2.0;
          PQNQDiff = NQ - PQ;
          if ( PQNQDiff.norm() < eps )
          {
            // std::cout << "Polar decomposition used " << ni << " iterations " << std::endl;
            break;
          } else {
            PQ = NQ;
          }
        }

        MatrixType33d QMatrix;
        QMatrix = NQ;
        // std::cout << "Initial Matrix = " << std::endl << M << std::endl;
        // std::cout << "Q Matrix = " << std::endl << QMatrix << std::endl;

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

  double mi = 0.0;
  MICostFunction::Pointer  costFunction = MICostFunction::New();

  costFunction->SetInputReferenceVolume(inputRefVolImage);

  for (auto i = inputSliceImages.begin(); i != inputSliceImages.end(); ++i) {
    costFunction->SetInputSliceImage( *i );
  }
  costFunction->SetFixedParameters( vrTransform->GetFixedParameters() );
  costFunction->SetNumberOfBins( nbins );

  finalVRTransform->SetFixedParameters( vrTransform->GetFixedParameters() );
  identityVRTransform->SetFixedParameters( vrTransform->GetFixedParameters() );

  /* The Versor Rigid 3D transform describes a rotation and a translation,
   * described with six parameters, 
   * with a fixed center of rotation, described by three parameters.
   * 
   * The first 3 elements are the components of the versor representation of 
   * 3D rotation. The last 3 parameters are the translation in each dimension.
   * The serialization of the fixed parameters is an array of 3 elements 
   *    defining the center of rotation.
   */


  /* We can exercise the function evaluation at the initial position:
    mi = alignmentEvaluation(
                  inputRefVolImage, 
                  inputSliceImages,
                  identityVRTransform);
    std::cout << "The identity MI value was : " << mi << std::endl << std::flush;

    std::cout << "Cost function at initial parameters: " << 
      costFunction->GetValue( vrTransform->GetParameters() ) << 
      std::endl << std::flush;

  std::cout << "Cost function evaluation was successful." << std::endl
    << std::flush;

  writeHistogramImage( allocImage( sumHistogramPtr ), "initialHistogram.png" );
  std::cout << "Stats of initial histogram: " << std::endl << std::flush;
  printHistogramImageStatistics( allocImage( sumHistogramPtr ) );
  */

//  double startTime = omp_get_wtime();
//  double finishTime = 0.0;

  // Let's try this with NLOPT instead.
  // nlopt::opt opt(nlopt::LN_COBYLA, 6);
//  nlopt::opt opt(nlopt::LN_BOBYQA, 6);
//   nlopt::opt opt(nlopt::LN_SBPLX, 6); // Slowest, but most accurate
  nlopt::algorithm algo;
  try {
    algo = getOptimizerAlgorithm(optimizerName);
  } catch (const std::invalid_argument& e) {
    std::cerr << e.what() << std::endl;
    std::exit(1);
  }
  nlopt::opt opt(algo, 6);
  std::cout << "nlopt optimizer : " << optimizerName << std::endl;

  // Set the initial step size for the optimizer.
  // https://nlopt.readthedocs.io/en/latest/NLopt_Reference/#initial-step-size
  std::vector<double> initial_step_sizes(6);
  // The first 3 elements are the components of the versor representation of
  // 3D rotation. The last 3 parameters are the translation in each dimension.
  // Heuristic for initial step size - set translation step size to 1/10th of the voxel dimension,
  // and rotation step sizes to an edge of brain adjusted similar step size.
  const ImageType::SpacingType & sp = inputRefVolImage->GetSpacing();
  const ImageType::RegionType region = inputRefVolImage->GetLargestPossibleRegion();
  // Set initial step sizes for translation parameters:
  initial_step_sizes[3] = sp[0] * 0.5;
  initial_step_sizes[4] = sp[1] * 0.5;
  initial_step_sizes[5] = sp[2] * 0.5;
  // Set initial step sizes for rotation parameters:
  // Heuristic - set the step size in radians to the number of radians
  // that is approximately a tenth of a voxel step at the edge of the image.
  // Estimate the diameter by a fraction of the distance across the field of view.
  // step size theta = arc length / radius.
  initial_step_sizes[0] = 0.5 * sp[0]/(region.GetSize()[0]*sp[0]*0.5);
  initial_step_sizes[1] = 0.5 * sp[1]/(region.GetSize()[1]*sp[1]*0.5);
  initial_step_sizes[2] = 0.5 * sp[2]/(region.GetSize()[2]*sp[2]*0.5);

  // Set the optimizer initial step sizes
  opt.set_initial_step( initial_step_sizes );
  std::cout << "Initial step size of parameters : " << std::endl;
  std::cout << " initial_step_sizes[0] : " << initial_step_sizes[0] << std::endl;
  std::cout << " initial_step_sizes[1] : " << initial_step_sizes[1] << std::endl;
  std::cout << " initial_step_sizes[2] : " << initial_step_sizes[2] << std::endl;
  std::cout << " initial_step_sizes[3] : " << initial_step_sizes[3] << std::endl;
  std::cout << " initial_step_sizes[4] : " << initial_step_sizes[4] << std::endl;
  std::cout << " initial_step_sizes[5] : " << initial_step_sizes[5] << std::endl;

  // Set the optimizer objective (cost) function
  struct my_func_data_struct myHolder;
  myHolder.costFunctionPointer = costFunction;
  opt.set_max_objective(my_cost_function, &myHolder);

  // JDA: Set the absolute parameter step tolerances
  std::vector<double> xtol_abs_values(6);
  xtol_abs_values[0] = 1.7e-4; // 0.01 deg in rads
  xtol_abs_values[1] = 1.7e-4;
  xtol_abs_values[2] = 1.7e-4;
  xtol_abs_values[3] = 1e-2; // 0.01 mm
  xtol_abs_values[4] = 1e-2;
  xtol_abs_values[5] = 1e-2;
  opt.set_xtol_abs(xtol_abs_values);
  std::cout << "Absolute tolerance on parameters : " << std::endl;
  std::cout << " xtol_abs_values[0] : " << xtol_abs_values[0] << std::endl;
  std::cout << " xtol_abs_values[1] : " << xtol_abs_values[1] << std::endl;
  std::cout << " xtol_abs_values[2] : " << xtol_abs_values[2] << std::endl;
  std::cout << " xtol_abs_values[3] : " << xtol_abs_values[3] << std::endl;
  std::cout << " xtol_abs_values[4] : " << xtol_abs_values[4] << std::endl;
  std::cout << " xtol_abs_values[5] : " << xtol_abs_values[5] << std::endl;

//  // Set relative tolerance on optimization parameters
//  double xtol_value = 1e-4;
//  opt.set_xtol_rel(xtol_value);
//  std::cout << "Relative tolerance on parameters value (xtol_rel) : " << xtol_value << std::endl;

  // Set absolute tolerance on objective function
  double ftol_value = 1e-4;
  opt.set_ftol_abs(ftol_value);
  std::cout << "Absolute tolerance on objective function (ftol_abs) : " << ftol_value << std::endl;

//  // Set relative tolerance on objective function
//  double ftol_value = 1e-4;
//  opt.set_ftol_rel(ftol_value);
//  std::cout << "Relative tolerance on function value (ftol_rel) : " << ftol_value << std::endl;

  // Set maximum number of optimizer evaluations
//  int maxeval_value = 1000;
  opt.set_maxeval(maxeval_value);
  std::cout << "Maximum number of iterations (maxeval) : " << maxeval_value << std::endl;

  std::vector<double> x(6);
  x[0] = vrTransform->GetParameters()[0];
  x[1] = vrTransform->GetParameters()[1];
  x[2] = vrTransform->GetParameters()[2];
  x[3] = vrTransform->GetParameters()[3];
  x[4] = vrTransform->GetParameters()[4];
  x[5] = vrTransform->GetParameters()[5];
  double maxf = 0.0;

  std::vector<double> finalFixedParams(3);
  finalFixedParams[0] = finalVRTransform->GetFixedParameters()[0];
  finalFixedParams[1] = finalVRTransform->GetFixedParameters()[1];
  finalFixedParams[2] = finalVRTransform->GetFixedParameters()[2];

  double startTime = omp_get_wtime();
  double finishTime = 0.0;

  try {
    nlopt::result result = opt.optimize(x, maxf);
    std::cout << "found optimum at f(" << 
        x[0] << "," << x[1] << "," <<
        x[2] << "," << x[3] << "," <<
        x[4] << "," << x[5] << ") = "
        << std::setprecision(10) << maxf << std::endl;
    std::cout << "with fixed parameters : (" <<
        finalFixedParams[0] << "," << finalFixedParams[1] << "," <<
        finalFixedParams[2] << ")"
        << std::setprecision(10) << std::endl;
    std::cout << "Number of function evals is : " << opt.get_numevals() << std::endl;
    std::cout << "Return code is : " << result << std::endl;
  } catch( std::exception &e) {
    std::cout << "nlopt failed." << e.what() << std::endl;
  }

  finishTime = omp_get_wtime();
  std::cout << "NLopt elapsed time (sec) : " << std::setprecision(10) << (finishTime - startTime) << std::endl;

  parmArray[0] = x[0];
  parmArray[1] = x[1];
  parmArray[2] = x[2];
  parmArray[3] = x[3];
  parmArray[4] = x[4];
  parmArray[5] = x[5];

  finalVRTransform->SetParameters( parmArray );

  /*
  writeHistogramImage( allocImage( sumHistogramPtr ), "finalHistogram.nrrd" );
  std::cout << "Stats of final histogram: " << std::endl;
  printHistogramImageStatistics( allocImage( sumHistogramPtr ) );
  */


  // Now write out the final VR Transform
  itk::TransformFileWriter::Pointer trsfWriter;
  trsfWriter = itk::TransformFileWriter::New();
  trsfWriter->SetInput( finalVRTransform );
  std::string outFile = "/data/alignTransform_" + outputTransformLabel + ".tfm";
  trsfWriter->SetFileName(outFile);

//  const char* dataDir = std::getenv("DATA_DIR");
//  std::string baseDir = (dataDir != nullptr) ? dataDir : "/data";
//  std::string outFile = baseDir + "/alignTransform_" + outputTransformLabel + ".tfm";
//  trsfWriter->SetFileName(outFile);

  try {
    trsfWriter->Update();
  }
  catch( itk::ExceptionObject & excp) {
    std::cerr << "Error while saving the transform to file " << std::endl;
    std::cerr << excp << std::endl;
    exit(1);
  }

//  // JDA: Write out the registration trace, RegistrationInfo vector to the parent input directory
//  std::string outputRegInfoPath = "regTrace_" + outputTransformLabel + ".csv";
//  std::ofstream regInfoFile(outputRegInfoPath);
//  if (regInfoFile.is_open()) {
//      // Write CSV header
//      regInfoFile << "EvaluatedTransformParams," // Evaluated transform parameters
//          << "FixedParams,"              // Fixed parameters
//          << "MutualInfo,"               // Mutual information
//          << "InitRuntime,"              // Pre-processing initialization runtime
////          << "time01,"                   // Time 01
////          << "time02,"                   // Time 02 to memset threads
////          << "time03,"                   // Time 03 to memset histograms and initialize counts to 0
//          << "ParallelRuntime,"          // Parallel section runtime
//          << "SerialRuntime"             // Serial section runtime
//          << std::endl;
//
//      for (const auto& info : registrationInfos) {
//          // Write evaluated transform params
//          for (size_t i = 0; i < info.evaluatedTransformParams.size(); ++i) {
//              regInfoFile << info.evaluatedTransformParams[i];
//              if (i < info.evaluatedTransformParams.size() - 1) {
//                  regInfoFile << " ";
//              }
//          }
//          regInfoFile << ",";
//          // Write fixed params
//          for (size_t i = 0; i < info.fixedTransformParams.size(); ++i) {
//              regInfoFile << info.fixedTransformParams[i];
//              if (i < info.fixedTransformParams.size() - 1) {
//                  regInfoFile << " ";
//              }
//          }
//          regInfoFile << ", ";
//          // Write mutual information
//          regInfoFile << info.mutualInformation;
//          regInfoFile << ", ";
//          // Write initialization runtime
//          regInfoFile << info.initTime;
//          regInfoFile << ", ";
////          // Write timestamp 01 runtime
////          regInfoFile << info.time01;
////          regInfoFile << ", ";
////          // Write timestamp 02 runtime
////          regInfoFile << info.time02;
////          regInfoFile << ", ";
////          // Write timestamp 03 runtime
////          regInfoFile << info.time03;
////          regInfoFile << ", ";
//          // Write parallel computation runtime
//          regInfoFile << info.parallelTime;
//          regInfoFile << ", ";
//          // Write serial computation runtime
//          regInfoFile << info.serialTime << std::endl;
//      }
//      regInfoFile.close();
//      std::cout << "Registration trace file saved: " << outputRegInfoPath << std::endl;
//  } else {
//      std::cerr << "Error: Unable to open file for writing registration info." << std::endl;
//      return 1;
//  }

//  // JDA: Write out the joint histogram data, HistogramData to the parent input directory
//    std::string outputHistoPath = "jointHistogram_" + outputTransformLabel + ".csv";
//    std::ofstream histoDataFile(outputHistoPath);
//    if (histoDataFile.is_open()) {
//        // Write CSV header
//        histoDataFile << "j_index,k_index,bincount" << std::endl;
//
//        for (const auto& info : histogramDatas) {
//            // Write J index from joint histogram
//            histoDataFile << info.xindex;
//            histoDataFile << ",";
//            // Write K index from joint histogram
//            histoDataFile << info.yindex;
//            histoDataFile << ",";
//            // Write bin count from J,K in joint histogram
//            histoDataFile << info.bincount << std::endl;
//        }
//        histoDataFile.close();
//        std::cout << "Joint histogram data saved : " << outputHistoPath << std::endl;
//    } else {
//        std::cerr << "Error: Unable to open file for writing histogram data." << std::endl;
//        return 1;
//    }


  exit(0);
}


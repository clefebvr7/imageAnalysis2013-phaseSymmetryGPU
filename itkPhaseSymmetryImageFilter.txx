
#ifndef __itkPhaseSymmetryImageFilter_txx
#define __itkPhaseSymmetryImageFilter_txx
#include "itkPhaseSymmetryImageFilter.h"
#include "itkImageFileWriter.h"
#include <ctime>
#include <iostream>
#include <vector>
#include <iostream>
#include "itkImage.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkImageSliceConstIteratorWithIndex.h"
#include "math.h"
#include "itkVnlFFTRealToComplexConjugateImageFilter.h"
#include "itkVnlFFTComplexConjugateToRealImageFilter.h"
#include "itkImageFileReader.h"
#include "itkConstantPadImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkSliceBySliceImageFilter.h"
#include "itkSimpleFilterWatcher.h"
#include "itkMultiplyImageFilter.h"
#include "itkFFTShiftImageFilter.h"
#include "itkComplexToRealImageFilter.h"
#include "itkFFTWComplexConjugateToRealImageFilter.h"
#include "itkFFTWRealToComplexConjugateImageFilter.h"
#include "itkFFTWComplexToComplexImageFilter.h"
#include "itkMultiplyByConstantImageFilter.h"
#include "itkAbsImageFilter.h"
#include "itkComplexToImaginaryImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkAddConstantToImageFilter.h"
#include "itkDivideImageFilter.h"
#include <algorithm>
#include <complex>
#include "itkFlipImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkImageFileWriter.h"
 
using namespace std;

typedef float  PixelType;
typedef complex<float>  PixelType2;
const unsigned int ImageDimension = 3;
typedef itk::Image< PixelType, ImageDimension > ImageType;

typedef itk::ImageSliceIteratorWithIndex< ImageType > IteratorType;
typedef itk::ImageSliceConstIteratorWithIndex< ImageType > ConstIteratorType;
typedef itk::FFTWRealToComplexConjugateImageFilter < PixelType, ImageDimension> FFTWFilterType;
typedef FFTWFilterType::OutputImageType    SpectralImageType;
typedef itk::ImageSliceIteratorWithIndex< SpectralImageType > IteratorType2;
typedef itk::ImageSliceConstIteratorWithIndex< SpectralImageType > ConstIteratorType2;

#define pi 3.1415926535897932384626433832795
vector<double> azimuths;
vector<double> elevations;


namespace itk
{

	template <class TInputImage, class TOutputImage>
	PhaseSymmetryImageFilter<TInputImage,TOutputImage>::PhaseSymmetryImageFilter()
	{
	//the default values, see itkPhaseSymmetryImageFilter.h for their definitions
		this->m_MinWaveLength = 25.0;
		this->m_SigmaAlpha = 0.45;
		this->m_Nscale = 1;
		this->m_Norient = 6;
		this->m_SigmaFreq = 0.25;
		this->m_Polarity = 1;
		this->m_PercentThreshold = .95;
		this->m_Mult = 3.0;
	}

	template <class TInputImage, class TOutputImage>
	void PhaseSymmetryImageFilter<TInputImage,TOutputImage>::Initialize( void )
	{
		typename TInputImage::SizeType  inputSize;
		typename TInputImage::IndexType inputIndex;
		typename Superclass::OutputImagePointer output = this->GetOutput();
		typename Superclass::InputImagePointer input = const_cast< TInputImage *>( this->GetInput() );

		inputIndex = input->GetLargestPossibleRegion().GetIndex();
		inputSize = input->GetLargestPossibleRegion().GetSize();
		int ndims = int(TInputImage::ImageDimension);
	}

/**************************************************************
********************Helper Functions START*********************
**************************************************************/

int round(double d){
    //round a double to the nearest integer
    return floor(d + 0.5);
}

void bankOrient(int norient){
    //initiate azimuths
    //initiate elevations
    double ArcL = pi / norient;
    for (int p=0; p <= norient/2; p++){
        double elev = p*ArcL;
        double rad = abs(cos(elev));
        double beta = ArcL/rad;
        int Mx = round(2*(pi/beta));
        if (1>Mx){
            Mx = 1;
        }
        for (int z=0; z <= Mx - 1; z++){
            elevations.push_back(elev);
            double az = z*2*pi/Mx;
            azimuths.push_back(az);
            if (p==0 && z==norient-1){
                break;
            }
        }
    }
}

double * sph2cart(double T, double P, int r){
    //sph2cart will convert spherical points to cartien points
    //cart[0] will be x component
    //cart[1] will be y component
    //cart[2] will be z component
    double *cart;
    cart = new double[3];

    cart[0] = sin(T)*cos(P);
    cart[1] = cos(T);
    cart[2] = -sin(T)*sin(P);
    return cart;
}


ImageType::Pointer SquareRootImage(ImageType::Pointer InputImage,ImageType::SizeType size){
    /*take the square root of each pixel, note each pixel is a complex number so sqrt (a + bj) where j
    is an imaginary is sqrt(a^2 + b^2)*/
    ImageType::RegionType region;
    region.SetSize( size );    

    IteratorType outIt( InputImage, region );
    outIt.GoToBegin();
    outIt.SetFirstDirection( 0 );  // 0=x, 1=y, 2=z
    outIt.SetSecondDirection( 1 ); // 0=x, 1=y, 2=z
        
    ImageType::IndexType index;
    
    while( !outIt.IsAtEnd() )
    {
        while( !outIt.IsAtEndOfSlice() )
        {
            while( !outIt.IsAtEndOfLine() )
            {
                float point = outIt.Get();
                float sqRoot = sqrt(point);
                outIt.Set( static_cast<PixelType> (sqRoot) );
                ++outIt;
            }
            outIt.NextLine();
        }
        outIt.NextSlice();
    }
    
    return InputImage;
}

ImageType::Pointer absImage(ImageType::Pointer InputImage,ImageType::SizeType size){
    /*take the absolute value of each pixel, by squaring and taking the square root of each pixel value*/
    ImageType::RegionType region;
    region.SetSize( size );    

    IteratorType outIt( InputImage, region );
    outIt.GoToBegin();
    outIt.SetFirstDirection( 0 );  // 0=x, 1=y, 2=z
    outIt.SetSecondDirection( 1 ); // 0=x, 1=y, 2=z
        
    ImageType::IndexType index;
    
    while( !outIt.IsAtEnd() )
    {
        while( !outIt.IsAtEndOfSlice() )
        {
            while( !outIt.IsAtEndOfLine() )
            {
                float point = outIt.Get();
                float sqRoot = sqrt(pow(point,2));
                outIt.Set( static_cast<PixelType> (sqRoot) );
                ++outIt;
            }
            outIt.NextLine();
        }
        outIt.NextSlice();
    }
    
    return InputImage;
}

ImageType::Pointer addContIm(ImageType::Pointer InputImage,ImageType::SizeType size, float AddConstant){
    //Add constant to each pixels in the image
    //NOTE: I did not use this
    ImageType::RegionType region;
    region.SetSize( size );    

    IteratorType outIt( InputImage, region );
    outIt.GoToBegin();
    outIt.SetFirstDirection( 0 );  // 0=x, 1=y, 2=z
    outIt.SetSecondDirection( 1 ); // 0=x, 1=y, 2=z
        
    ImageType::IndexType index;
    
    while( !outIt.IsAtEnd() )
    {
        while( !outIt.IsAtEndOfSlice() )
        {
            while( !outIt.IsAtEndOfLine() )
            {
                float point = outIt.Get();
                float AddPoint = point + AddConstant;
                outIt.Set( static_cast<PixelType> (AddPoint) );
                ++outIt;
            }
            outIt.NextLine();
        }
        outIt.NextSlice();
    }
    return InputImage;
}


ImageType::Pointer maxPixelZero(ImageType::Pointer InputImage,ImageType::SizeType size, float Tval){
    //Take the maximum value between pixel value minus a value and zero
    ImageType::RegionType region;
    region.SetSize( size );    

    IteratorType outIt( InputImage, region );
    outIt.GoToBegin();
    outIt.SetFirstDirection( 0 );  // 0=x, 1=y, 2=z
    outIt.SetSecondDirection( 1 ); // 0=x, 1=y, 2=z
        
    ImageType::IndexType index;
    
    while( !outIt.IsAtEnd() )
    {
        while( !outIt.IsAtEndOfSlice() )
        {
            while( !outIt.IsAtEndOfLine() )
            {
                float point = outIt.Get();
                float z = 0.0;
                float maxVal = max((point - Tval), z);
                outIt.Set( static_cast<PixelType> (maxVal) );
                ++outIt;
            }
            outIt.NextLine();
        }
        outIt.NextSlice();
    }
    
    return InputImage;
}

ImageType::Pointer maxPixelZeroMod(ImageType::Pointer InputImage,ImageType::SizeType size, float Tval){
    //Take the maximum value between pixel value minus a value and zero. Note: Only needed if using Ilkers' threshold method
    ImageType::RegionType region;
    region.SetSize( size );    

    IteratorType outIt( InputImage, region );
    outIt.GoToBegin();
    outIt.SetFirstDirection( 0 );  // 0=x, 1=y, 2=z
    outIt.SetSecondDirection( 1 ); // 0=x, 1=y, 2=z
        
    ImageType::IndexType index;
    
    while( !outIt.IsAtEnd() )
    {
        while( !outIt.IsAtEndOfSlice() )
        {
            while( !outIt.IsAtEndOfLine() )
            {
                float point = outIt.Get();
                float z = 0.0;
                if (point<Tval){
                    outIt.Set( static_cast<PixelType> (z) );
                }
                    
                ++outIt;
            }
            outIt.NextLine();
        }
        outIt.NextSlice();
    }
    
    return InputImage;
}

float sumPixelImagePow(ImageType::Pointer InputImage, ImageType::SizeType size, int powNum){
    //take each pixel to the power of powNum and then sum all the pixels together

    ImageType::RegionType region;
    region.SetSize( size );
    IteratorType outIt( InputImage, region );
    outIt.GoToBegin();
    outIt.SetFirstDirection( 0 );  // 0=x, 1=y, 2=z
    outIt.SetSecondDirection( 1 ); // 0=x, 1=y, 2=z
                    
    ImageType::IndexType index;
        
    float sum = 0;
    while( !outIt.IsAtEnd() )
    {
        while( !outIt.IsAtEndOfSlice() )
        {
            while( !outIt.IsAtEndOfLine() )
            {
                float point = outIt.Get();
                sum = sum + pow(point, powNum);
                ++outIt;
            }
            outIt.NextLine();
        }
        outIt.NextSlice();
    }
    return sum;
}

// Simple swap function for our in place swapping.
void swap(float &val1, float &val2)
{
    //Quicksort helper function
    float temp = val1;
    val1 = val2;
    val2 = temp;
}

int median3(float *arIntegers,int left,int right)
{
    //Quicksort helper function
     int center = (left+right)/2;

   if(arIntegers[center] < arIntegers[left])
       swap(arIntegers[left],arIntegers[center]);
   if(arIntegers[right] < arIntegers[left])
       swap(arIntegers[left],arIntegers[right]);
   if(arIntegers[right] < arIntegers[center])
       swap(arIntegers[center],arIntegers[right]);

   swap(arIntegers[center],arIntegers[right-1]);

   return center;
}

// This function takes an array (or one half an array) and sorts it.
// It then returns a new pivot index number back to quicksort.

int partition(float *arIntegers, int left, int right, int pivot)
{
    //Quicksort helper function
     float pivotValue = arIntegers[pivot];

     // Swap it out all the way to the end of the array
     // So we know where it always is.
     swap(arIntegers[pivot], arIntegers[right]);
     int storeIndex = left;

     // Move through the array from start to finish comparing each to our
     // pivot value (not index, the value that was located at the pivot index)
     for (int i = left; i < right; i++)
     {
         if (arIntegers[i] <= pivotValue)
         {
             swap(arIntegers[i], arIntegers[storeIndex]);
             storeIndex++;
         }
     }
     swap(arIntegers[storeIndex], arIntegers[right]);
     return storeIndex;
}

// Quicksort controller function, it partitions the different pieces of our array.
void quicksort(float *arIntegers, int left, int right)
{
    //Quick sort
    //note: I obtain quick sort code and all its helper functions from http://www.cplusplus.com/forum/beginner/9388/
    if (right > left)
    {
         int pivotIndex = median3(arIntegers,left,right);
         int pivotNewIndex = partition(arIntegers, left, right, pivotIndex);

         // Recursive call to quicksort to sort each half.
         quicksort(arIntegers, left, pivotNewIndex-1);
         quicksort(arIntegers, pivotNewIndex+1, right);
    }
}




float findMedian(ImageType::Pointer InputImage,ImageType::SizeType size, int NumberOfPixel, double PercentThreshold){
    //Find the median pixel from the whole image

    float *fullPixel = new float[NumberOfPixel];

    //float fullPixel[NumberOfPixel];

    ImageType::RegionType region;
    region.SetSize( size );    

    typedef itk::ImageSliceIteratorWithIndex< ImageType > IteratorType;
    typedef itk::ImageSliceConstIteratorWithIndex< ImageType > ConstIteratorType;
    IteratorType outIt( InputImage, region );
    outIt.GoToBegin();
    outIt.SetFirstDirection( 0 );  // 0=x, 1=y, 2=z
    outIt.SetSecondDirection( 1 ); // 0=x, 1=y, 2=z
        
    ImageType::IndexType index;
    int i =0;
    while( !outIt.IsAtEnd() )
    {
        while( !outIt.IsAtEndOfSlice() )
        {
            while( !outIt.IsAtEndOfLine() )
            {
                float point = outIt.Get();
                fullPixel[i] = point;
                i+=1;
                ++outIt;
            }
            outIt.NextLine();
        }
        outIt.NextSlice();
    }

    quicksort(fullPixel, 0, NumberOfPixel - 1);
    
    int ind = round(NumberOfPixel*PercentThreshold);

    return fullPixel[ind];
}

/**************************************************************
********************Helper Functions END*********************
**************************************************************/

	template <class TInputImage, class TOutputImage>
	void PhaseSymmetryImageFilter<TInputImage,TOutputImage>::GenerateData( void )
	{
	
		typename TInputImage::SizeType  inputSize;
		typename TInputImage::IndexType inputIndex;
		typename Superclass::OutputImagePointer output = this->GetOutput();
		typename Superclass::InputImagePointer input = const_cast< TInputImage *>( this->GetInput() );

		inputIndex = input->GetLargestPossibleRegion().GetIndex();
		inputSize = input->GetLargestPossibleRegion().GetSize();
		int ndims = int(TInputImage::ImageDimension);

	/*******************************************
    Define all functions
    *******************************************/
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typedef itk::VnlFFTRealToComplexConjugateImageFilter<PixelType, ImageDimension>  FFTFilterType;
    typedef itk::MultiplyByConstantImageFilter<ImageType, unsigned char, ImageType> MultiplyConstantFilterType;
    typedef itk::MultiplyImageFilter <ImageType, ImageType > MultiplyImageFilterType;

    typedef itk::ComplexToImaginaryImageFilter<SpectralImageType, ImageType> ImaginaryFilterType3;
    typedef itk::ComplexToRealImageFilter<SpectralImageType, ImageType> RealFilterType3;
    typedef itk::FFTShiftImageFilter<ImageType,ImageType> ShiftFilterType;

    typedef itk::MultiplyImageFilter <SpectralImageType, ImageType> MultiplyImageFilterType2;
    typedef itk::FFTWComplexToComplexImageFilter < PixelType > complexFFTFilterType;
    typedef itk::AddImageFilter <ImageType, ImageType > AddImageFilterType;
    typedef itk::SubtractImageFilter <ImageType, ImageType > SubtractImageFilterType;
    typedef itk::AddConstantToImageFilter <ImageType, float, ImageType> AddConstantToImageFilterType2;
    typedef itk::DivideImageFilter <ImageType, ImageType, ImageType > DivideImageFilterType;

    /*******************************************
    Calculate Filter's Radial and Angular Component
    *******************************************/

	
	double radius; // Radial frequency coordinate
    double z, y, x, rfo;
    double denomRad = (2.0 * pow(log(this->m_SigmaFreq), 2.0));
	int cols = input->GetLargestPossibleRegion().GetSize()[1];
	int rows = input->GetLargestPossibleRegion().GetSize()[0];
	int slices = input->GetLargestPossibleRegion().GetSize()[2];
    double denomSpr = (2 * pow(this->m_SigmaAlpha,2));
    bankOrient(this->m_Norient);
    int sAz = azimuths.size();

    ImageType::Pointer *radialImages;
    radialImages = new ImageType::Pointer[this->m_Nscale];
    
    ImageType::Pointer *angularImages;
    angularImages = new ImageType::Pointer[sAz];

    ImageType::SizeType size;

    size[0] = cols;
    size[1] = rows;
    size[2] = slices;

    ImageType::RegionType region;
    region.SetSize( size );

    int h = 0;
    int maxV = max(sAz, this->m_Nscale);
    for (h; h < maxV; h++){

        ImageType::Pointer angGab = ImageType::New();
        ImageType::Pointer radGab = ImageType::New();
        ImageType::RegionType region;
        region.SetSize( size );
        ImageType::RegionType region2;
        region2.SetSize( size );

        angGab->SetLargestPossibleRegion( region );
        angGab->SetBufferedRegion( region );
        angGab->SetRequestedRegion( region );
        angGab->Allocate();

        radGab->SetLargestPossibleRegion( region2 );
        radGab->SetBufferedRegion( region2 );
        radGab->SetRequestedRegion( region2 );
        radGab->Allocate();


        IteratorType outIt( radGab, region );
        outIt.GoToBegin();
        outIt.SetFirstDirection( 0 );  // 0=x, 1=y, 2=z
        outIt.SetSecondDirection( 1 ); // 0=x, 1=y, 2=z
        ImageType::IndexType index;
        
        IteratorType outIt2( angGab, region );
        outIt2.GoToBegin();
        outIt2.SetFirstDirection( 0 );  // 0=x, 1=y, 2=z
        outIt2.SetSecondDirection( 1 ); // 0=x, 1=y, 2=z
        ImageType::IndexType index2;



        int i =0;
        while( !outIt.IsAtEnd() )
        {
            z = ((-slices/2.0)+i)/(slices/2);
            int j = 0;
            while( !outIt.IsAtEndOfSlice() )
            {
                y = -((-rows/2.0)+j)/(rows/2);
                int k = 0;
                while( !outIt.IsAtEndOfLine() )
                {
                    x = ((-cols/2.0)+k)/(cols/2);
                    
                    radius = sqrt(pow(x,2) + pow(y,2) + pow(z,2));

                    if (radius == 0){
                        radius =0.0001;
                    }

                    if (h < this->m_Nscale){
                        //Filter Radial Component
                        double point;
                        if((i==(slices/2)+1) && (k==(cols/2)+1) && (j==(rows/2)+1)){
                            point=0;
                        }
                        else{
							double m = this->m_Mult;
							double WaveLength = this->m_MinWaveLength * pow(m,h-1);
                            rfo = 2.0/WaveLength;
                            point = exp(-pow(log(radius/rfo), 2.0) / denomRad);
                            outIt.Set( static_cast<PixelType> (point) );
                        }
                        outIt.Set( static_cast<PixelType> (point) );
                    }

                    if (h < sAz){
                        //Filter's angular component
                        double vx = cos(azimuths[h])*cos(elevations[h]);
                        double vy = sin(azimuths[h])*cos(elevations[h]);
                        double vz = sin(elevations[h]);
                        double dAlpha = acos(((vy*y) + (vx*x) + (vz*z))/radius);
                        double point = exp(-pow(dAlpha,2)/denomSpr);
                        outIt2.Set( static_cast<PixelType> (point) );
        
                    }

                    k+=1;
                    ++outIt,++outIt2;
                }
                j+=1;
                outIt.NextLine();
                outIt2.NextLine();
            }
            i+=1;
            outIt.NextSlice();
            outIt2.NextSlice();
        }

        if (h < this->m_Nscale){
            radialImages[h] = radGab;
        }
        if (h < sAz){
            angularImages[h] = angGab;
        }
    }


    /*******************************************
    ********************************************
    ********************************************
    Phase Symmetry
    ********************************************
    *******************************************
    *******************************************/

    /*******************************************
    REAL IMAGE - VNL FFT
    *******************************************/

    //First goal is to take the fft of the given image

    FFTFilterType::Pointer fftFilter = FFTFilterType::New();
    fftFilter->SetInput( input );
  //  fftFilter->SetNumberOfThreads(64);
    try
    {
        fftFilter->Update();
    }
    catch( itk::ExceptionObject & excp )
    {
        std::cerr << "Error: " << std::endl;
        std::cerr << excp << std::endl;
    }                            

    /*******************************************
    Create Empty Matrice for the Sum of amplitude response (sumAn_ThisOrient), so that we can keep count
    *******************************************/

    int orientnew [15] = { 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 };
    int sizeO = 15;
    
    int firstTE = 1;
    ImageType::Pointer totalEnergy = ImageType::New();// Keep track of accumulating weighted phase congruency values (energy)

    totalEnergy->SetLargestPossibleRegion( region );
    totalEnergy->SetBufferedRegion( region );
    totalEnergy->SetRequestedRegion( region );
    totalEnergy->Allocate();

    int firstTSA = 1;
    ImageType::Pointer totalSumAn = ImageType::New();//Keep track of accumulating filter response amplitude values.

    totalSumAn->SetLargestPossibleRegion( region );
    totalSumAn->SetBufferedRegion( region );
    totalSumAn->SetRequestedRegion( region );
    totalSumAn->Allocate();

/* //Ilkers Noise Compensation
	SpectralImageType2::Pointer **EO;
	EO = new SpectralImageType2::Pointer*[sizeO];
	ImageType::Pointer * ifftFilterArray;
	ifftFilterArray = new ImageType::Pointer[nscale];
*/

    ImageType::Pointer * absEOScale1;
    absEOScale1 = new ImageType::Pointer[sizeO]; // save the absolute for EO for scale 1 and for all octaves

    int o =1;
    for (o; o <= sizeO; o++){
        int s =1;

        ImageType::Pointer * imaginaryImageGaborList;
        imaginaryImageGaborList = new ImageType::Pointer[this->m_Nscale];

        ImageType::Pointer * realImageGaborList;
        realImageGaborList = new ImageType::Pointer[this->m_Nscale];

        int firstSTO = 1;
        ImageType::Pointer sumAn_ThisOrient = ImageType::New();

        int firstETO = 1;
        ImageType::Pointer Energy_ThisOrient = ImageType::New();
		
/* //Ilkers Noise Compensation
		EO[o-1] = new SpectralImageType2::Pointer[this->m_Nscale];
		float EM_n=0;
*/

        MultiplyImageFilterType::Pointer RealSquared = MultiplyImageFilterType::New ();

        ImaginaryFilterType3::Pointer ImaginaryImageGabor = ImaginaryFilterType3::New();

        MultiplyImageFilterType::Pointer ImaginarySquared = MultiplyImageFilterType::New ();

        RealFilterType3::Pointer RealImageGabor = RealFilterType3::New();

        for (s; s <= this->m_Nscale; s++){ 
		//for each scale
            /*******************************************
            Create Filter - Multiplying radial component by angular component
            *******************************************/

            MultiplyImageFilterType::Pointer gFilter = MultiplyImageFilterType::New ();
            gFilter->SetInput1(angularImages[orientnew[o-1]]);															//CRASHING HERE!!!!!!!!!!!
            gFilter->SetInput2(radialImages[s-1]);
       //     gFilter->SetNumberOf(64);
            gFilter->Update();
                
            /*******************************************
            Shift filter (filter that was multiplied by radual and angular components)
            *******************************************/
            
            ShiftFilterType::Pointer shifter = ShiftFilterType::New();

            shifter->SetInput(gFilter->GetOutput());
          //  shifter->SetNumberOfThreads(64);

/* //Ilkers Noise Compensation
             typedef itk::VnlFFTRealToComplexConjugateImageFilter<PixelType, ImageDimension>  FFTFilterType;
            FFTFilterType::Pointer inverseFFTFilter = FFTFilterType::New();
            inverseFFTFilter->SetInput( shifter->GetOutput() );
            inverseFFTFilter->Update();

            typedef itk::ComplexToRealImageFilter<SpectralImageType2, ImageType> RealFilterType;
            RealFilterType::Pointer RealInverseFFTFilter = RealFilterType::New();
            RealInverseFFTFilter->SetInput(inverseFFTFilter->GetOutput());

            typedef itk::MultiplyByConstantImageFilter<ImageType, float, ImageType> MultiplyConstantFilterType;
            MultiplyConstantFilterType::Pointer ifftFil = MultiplyConstantFilterType::New();
            ifftFil->SetInput(reader2->GetOutput());
            float val = rows*cols*slices;
            float sqVal = sqrt(val);
            ifftFil->SetConstant(sqVal);
            ifftFil->Update();

            ifftFilterArray[s-1] = ifftFil->GetOutput();
*/

            /*******************************************
            Multiply the log-gabor filter and the fft of the image
            Then take the products ifft
            *******************************************/
                                                                            
            MultiplyImageFilterType2::Pointer imNfi = MultiplyImageFilterType2::New ();
            imNfi->SetInput1(fftFilter->GetOutput());
            imNfi->SetInput2(shifter->GetOutput());
         //   imNfi->SetNumberOfThreads(64);
            imNfi->Update();

            /*******************************************
            Take the inverse of (Log-Gabor)*(FFT Image)
            *******************************************/
            complexFFTFilterType::Pointer complexFFT = complexFFTFilterType::New();
            complexFFT->SetInput( imNfi->GetOutput() ); // compute forward FFT
            //complexFFT->SetInput( imNfi ); // compute forward FFT
            //complexFFT->SetNumberOfThreads(64);
            complexFFT->Update();

/* //Ilkers Noise Compensation
			EO[o-1][s-1] = complexFFT->GetOutput();
*/

            /*******************************************
            Update sumAn_ThisOrient
                -sqrt((realPart)^2 + (imaginaryPart)^2)
            *******************************************/

            RealImageGabor->SetInput(complexFFT->GetOutput());
           // RealImageGabor->SetNumberOfThreads(64);
            ImaginaryImageGabor->SetInput(complexFFT->GetOutput());
         //  ImaginaryImageGabor->SetNumberOf(64);

            //Save these values, they will be used again in polarity and save time not computing them
            imaginaryImageGaborList[s-1] = ImaginaryImageGabor->GetOutput();
            realImageGaborList[s-1] = RealImageGabor->GetOutput();

			//Square the real values
            RealSquared->SetInput1(RealImageGabor->GetOutput());
            RealSquared->SetInput2(RealImageGabor->GetOutput());
            //RealSquared->SetNumberOfThreads(64);

			//Square the imaginary values
            ImaginarySquared->SetInput1(ImaginaryImageGabor->GetOutput());
            ImaginarySquared->SetInput2(ImaginaryImageGabor->GetOutput());
           // ImaginarySquared->SetNumberOf(64);


		   //Add the squared real values, and the squared imaginary values
            AddImageFilterType::Pointer addImagReal = AddImageFilterType::New ();
            addImagReal->SetInput1(RealSquared->GetOutput());
            addImagReal->SetInput2(ImaginarySquared->GetOutput());
            //addImagReal->SetNumberOfThreads(64);
            addImagReal->Update();


            ImageType::Pointer An = SquareRootImage(addImagReal->GetOutput(), size);

            //add the calculated An to the sumAn_ThisOrient

            if (s == 1){
                absEOScale1[o-1] = An;    //need for noise compensation
            }

            if (firstSTO == 1){
                sumAn_ThisOrient = An;
                firstSTO = 0;
            }
            else{
                AddImageFilterType::Pointer addAn = AddImageFilterType::New ();
                addAn->SetInput1(sumAn_ThisOrient);
                addAn->SetInput2(An);
               // addAn->SetNumberOfThreads(64);
                addAn->Update();
                sumAn_ThisOrient = addAn->GetOutput();
            }

            /*******************************************
            Squares each pixles and then sum all the pixels in the image (only if s==1)
            *******************************************/
/* //Ilkers Noise Compensation
            if (s==1){

                EM_n = sumPixelImagePow(shifter->GetOutput(), size, 2);
            }
*/
        }

        /*******************************************
        Calculate the phase symmetry measure
        note: we already calculated
            real(EO[o][s]) with the variable name RealImageGabor->GetOutput()
            imaginary(EO[o][s]) with the variable name ImaginaryImageGabor->GetOutput()
            (real(EO[o][s]))^2 with variable name ImaginarySquared->GetOutput()
            (imaginary(EO[o][s]))^2 with variable name ImaginarySquared->GetOutput()
            to get the absolute value for each, we only need to call the SquareRootImage
        *******************************************/
        if (this->m_Polarity==0){
		//look for 'white' and 'black' spots
            s=1;
            for (s; s <= this->m_Nscale; s++){

                ImageType::Pointer absReal = absImage(realImageGaborList[s-1], size);
                ImageType::Pointer absImaginary = absImage(imaginaryImageGaborList[s-1], size);

                SubtractImageFilterType::Pointer subtractFilter = SubtractImageFilterType::New ();
                subtractFilter->SetInput1(absReal);
                subtractFilter->SetInput2(absImaginary);
//                subtractFilter->SetNumberOfThreads(64);
                subtractFilter->Update();

                if (firstETO==1){
                    Energy_ThisOrient = subtractFilter->GetOutput();
                    firstETO = 0;
                }
                else{
                    AddImageFilterType::Pointer addEnergy = AddImageFilterType::New ();
                    addEnergy->SetInput1(Energy_ThisOrient);
                    addEnergy->SetInput2(subtractFilter->GetOutput());
                 //   addEnergy->SetNumberOfThreads(64);
                    addEnergy->Update();
                    Energy_ThisOrient = addEnergy->GetOutput();
                }
            }
        }

        else if (this->m_Polarity==1){
		//Just look for 'white' spots
            s=1;
            for (s; s <= this->m_Nscale; s++){

                ImageType::Pointer absImaginary = absImage(imaginaryImageGaborList[s-1], size);

                SubtractImageFilterType::Pointer subtractFilter = SubtractImageFilterType::New ();
                subtractFilter->SetInput1(realImageGaborList[s-1]);
                subtractFilter->SetInput2(absImaginary);
                //subtractFilter->SetNumberOfThreads(64);
                subtractFilter->Update();

                if (firstETO==1){
                    Energy_ThisOrient = subtractFilter->GetOutput();
                    firstETO = 0;
                }
                else{
                    AddImageFilterType::Pointer addEnergy = AddImageFilterType::New ();
                    addEnergy->SetInput1(Energy_ThisOrient);
                    addEnergy->SetInput2(subtractFilter->GetOutput());
                  //  addEnergy->SetNumberOfThreads(64);
                    addEnergy->Update();
                    Energy_ThisOrient = addEnergy->GetOutput();
                }
            }
        }

        else if (this->m_Polarity==-1){
		//Just look for 'black' spots
            s=1;
            for (s; s <= this->m_Nscale; s++){

                ImageType::Pointer absImaginary = absImage(imaginaryImageGaborList[s-1], size);

                SubtractImageFilterType::Pointer subtractFilter = SubtractImageFilterType::New ();
                subtractFilter->SetInput1(realImageGaborList[s-1]);
                subtractFilter->SetInput2(absImaginary);
          //      subtractFilter->SetNumberOfThreads(64);
                subtractFilter->Update();

                if (firstETO==1){
                    Energy_ThisOrient = subtractFilter->GetOutput();
                    firstETO = 0;
                }
                else{
                    SubtractImageFilterType::Pointer subtractFilter2 = SubtractImageFilterType::New ();
                    subtractFilter2->SetInput1(Energy_ThisOrient);
                    subtractFilter2->SetInput2(subtractFilter->GetOutput());
                //    subtractFilter2->SetNumberOfThreads(64);
                    subtractFilter2->Update();
                    Energy_ThisOrient = subtractFilter2->GetOutput();
                }
            }
        }

		/* //Ilkers Noise Compensation        
        typedef itk::MultiplyImageFilter <ImageType, ImageType> MultiplyImageFilterType5;
        MultiplyImageFilterType5::Pointer absEOSquared = MultiplyImageFilterType5::New ();
        absEOSquared->SetInput1(absEOScale1[o-1]);
        absEOSquared->SetInput2(absEOScale1[o-1]);
        absEOSquared->Update();

        float medianE2n = findMedian(absEOSquared->GetOutput(), size, rows*cols*slices);
        float meanE2n = (-medianE2n)/log(0.5);
        float noisePower = meanE2n/EM_n;

        ImageType::Pointer EstSumAn2 = blankIm->GetOutput();
        s=1;    
        for (s; s <= nscale; s++){
        
            typedef itk::MultiplyImageFilter <ImageType, ImageType> MultiplyImageFilterType5;
            MultiplyImageFilterType5::Pointer ifftFilterSquared = MultiplyImageFilterType5::New ();
            ifftFilterSquared->SetInput1(ifftFilterArray[s-1]);
            ifftFilterSquared->SetInput2(ifftFilterArray[s-1]);

            typedef itk::AddImageFilter <ImageType, ImageType > AddImageFilterType2;
            AddImageFilterType2::Pointer addifftFilt = AddImageFilterType2::New ();
            addifftFilt->SetInput1(EstSumAn2);
            addifftFilt->SetInput2(ifftFilterSquared->GetOutput());
            addifftFilt->Update();

            EstSumAn2 = addifftFilt->GetOutput();
        }

        ImageType::Pointer EstSumAiAj = blankIm->GetOutput();
        int si=1;    
        for (si; si <= (nscale - 1); si++){
            int sj= si + 1;
            for (sj; sj <= nscale; sj++){

                typedef itk::MultiplyImageFilter <ImageType, ImageType> MultiplyImageFilterType5;
                MultiplyImageFilterType5::Pointer multiSiSj = MultiplyImageFilterType5::New ();
                multiSiSj->SetInput1(ifftFilterArray[si-1]);
                multiSiSj->SetInput2(ifftFilterArray[sj-1]);

                typedef itk::AddImageFilter <ImageType, ImageType > AddImageFilterType2;
                AddImageFilterType2::Pointer addSiSjEst = AddImageFilterType2::New ();
                addSiSjEst->SetInput1(EstSumAiAj);
                addSiSjEst->SetInput2(multiSiSj->GetOutput());
                addSiSjEst->Update();

                EstSumAiAj = addSiSjEst->GetOutput();
            }
        }

        float sumEstSumAn2 = sumPixelImagePow(EstSumAn2, size, 1);
        float sumEstSumAiAj = sumPixelImagePow(EstSumAiAj, size, 1);

        //DO calculate EstNoiseEnergy2   --------  absEOScale1[o-1]

        float EstNoiseEnergy2 = (2*noisePower*sumEstSumAn2) + (4*noisePower*sumEstSumAiAj);

        float tau = sqrt(EstNoiseEnergy2/2);
        float EstNoiseEnergy = tau*(sqrt(pi/2));
        float EstNoiseEnergySigma = sqrt((2-(pi/2))*(pow(tau,2)));

        float T = EstNoiseEnergy + (k*EstNoiseEnergySigma);

        T = (T/1.7);
		
		Energy_ThisOrient = maxPixelZero(Energy_ThisOrient, size, T);
*/

		
		float thres = findMedian(Energy_ThisOrient, size, rows*cols*slices, this->m_PercentThreshold);
		Energy_ThisOrient = maxPixelZeroMod(Energy_ThisOrient, size, thres);
        

        /*******************************************
        Update accumulator for sumAn and totalEnergy
        *******************************************/
        if (firstTSA==1){
            totalSumAn = sumAn_ThisOrient;
            firstTSA=0;
        }
        else{
            AddImageFilterType::Pointer NewSumAn = AddImageFilterType::New ();
            NewSumAn->SetInput1(totalSumAn);
            NewSumAn->SetInput2(sumAn_ThisOrient);
          //  NewSumAn->SetNumberOfThreads(64);
            NewSumAn->Update();        
            totalSumAn = NewSumAn->GetOutput();
        }

        if (firstTE==1){
            totalEnergy = Energy_ThisOrient;
            firstTE=0;
        }
        else{
            AddImageFilterType::Pointer NewEnergy = AddImageFilterType::New ();
            NewEnergy->SetInput1(totalEnergy);
            NewEnergy->SetInput2(Energy_ThisOrient);
        //    NewEnergy->SetNumberOfThreads(64);
            NewEnergy->Update();        
            totalEnergy = NewEnergy->GetOutput();
        }
    }        

    float epsilon = 0.0001;

    AddConstantToImageFilterType2::Pointer finalTotalSumAn = AddConstantToImageFilterType2::New();
    finalTotalSumAn->SetInput(totalSumAn);
    finalTotalSumAn->SetConstant(epsilon);
 //   finalTotalSumAn->SetNumberOfThreads(64);
    finalTotalSumAn->Update();

    DivideImageFilterType::Pointer phaseSymmIm = DivideImageFilterType::New ();
    phaseSymmIm->SetInput1(totalEnergy);
    phaseSymmIm->SetInput2(finalTotalSumAn->GetOutput());
    //phaseSymmIm->SetNumberOfThreads(64);
    phaseSymmIm->Update();

	this->GraftOutput(phaseSymmIm->GetOutput());

		itkDebugMacro("GenerateOutputInformation End");
	}


	template <class TInputImage, class TOutputImage>
	void PhaseSymmetryImageFilter<TInputImage,TOutputImage>::GenerateInputRequestedRegion()
	{
		itkDebugMacro("GenerateInputRequestedRegion Start");
		Superclass::GenerateInputRequestedRegion();

		if ( this->GetInput() )
		{
			typename TInputImage::RegionType RequestedRegion;
			typename TInputImage::SizeType  inputSize;
			typename TInputImage::IndexType inputIndex;
			typename TInputImage::SizeType  inputLargSize;
			typename TInputImage::IndexType inputLargIndex;
			typename TOutputImage::SizeType  outputSize;
			typename TOutputImage::IndexType outputIndex;

			outputIndex = this->GetOutput()->GetRequestedRegion().GetIndex();
			outputSize = this->GetOutput()->GetRequestedRegion().GetSize();
			inputLargSize = this->GetInput()->GetLargestPossibleRegion().GetSize();
			inputLargIndex = this->GetInput()->GetLargestPossibleRegion().GetIndex();

			for(unsigned int i=0; i<TInputImage::ImageDimension; i++)
			{
				inputSize[i] = outputSize[i];
				inputIndex[i] = outputIndex[i];
			}

			RequestedRegion.SetSize(inputSize);
			RequestedRegion.SetIndex(inputIndex);
			InputImagePointer input = const_cast< TInputImage *> ( this->GetInput() );
			input->SetRequestedRegion (RequestedRegion);
		}


		itkDebugMacro("GenerateInputRequestedRegion End");
	}


	/**
	* GenerateData Performs the accumulation
	*/
	template <class TInputImage, class TOutputImage>
	void PhaseSymmetryImageFilter<TInputImage,TOutputImage>::GenerateOutputInformation( void )
	{
		typename TOutputImage::RegionType outputRegion;
		typename TInputImage::IndexType inputIndex;
		typename TInputImage::SizeType  inputSize;
		typename TOutputImage::SizeType  outputSize;
		typename TOutputImage::IndexType outputIndex;
		typename TInputImage::SpacingType inSpacing;
		typename TInputImage::PointType inOrigin;
		typename TOutputImage::SpacingType outSpacing;
		typename TOutputImage::PointType outOrigin;

		// Get pointers to the input and output
		typename Superclass::OutputImagePointer output = this->GetOutput();
		typename Superclass::InputImagePointer input = const_cast< TInputImage *>( this->GetInput() );

		//Return if input and output are both null
		if( !input || !output )
		{
			return;
		}

		inputIndex = input->GetLargestPossibleRegion().GetIndex();
		inputSize = input->GetLargestPossibleRegion().GetSize();
		inSpacing = input->GetSpacing();
		inOrigin = input->GetOrigin();

		// Set the LargestPossibleRegion of the output.

		for(unsigned int i = 0; i<InputImageDimension; i++)
		{
			outputSize[i]  = inputSize[i];
			outputIndex[i] = inputIndex[i];
			outSpacing[i] = inSpacing[i];
			outOrigin[i]  = inOrigin[i];
		}

		//Set the size of the output region
		outputRegion.SetSize(outputSize);
		//Set the index of the output region
		outputRegion.SetIndex(outputIndex);
		//Set the origin and spacing
		output->SetOrigin(outOrigin);
		output->SetSpacing(outSpacing);
		//Set the largest po
		output->SetLargestPossibleRegion(outputRegion);
	}


	template <class TInputImage, class TOutputImage>
	void PhaseSymmetryImageFilter<TInputImage,TOutputImage>::PrintSelf(std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf(os,indent);

		//  os << indent << " Integral Filter Normalize By: " << m_Cutoff << std::endl;

	}



} // end namespace itk




#endif

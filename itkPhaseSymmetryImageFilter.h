#ifndef __itkPhaseSymmetryImageFilter_h
#define __itkPhaseSymmetryImageFilter_h

#include "itkRealAndImaginaryToComplexImageFilter.h"

#include <vector>
#include <complex>

namespace itk
{


	template <class TInputImage, class TOutputImage>
	class ITK_EXPORT PhaseSymmetryImageFilter : public ImageToImageFilter<TInputImage,TOutputImage>
	{
	public:
		/** Standard class typedefs. */
		typedef PhaseSymmetryImageFilter                           Self;
		typedef ImageToImageFilter<TInputImage,TOutputImage>  Superclass;
		typedef SmartPointer<Self>                            Pointer;
		typedef SmartPointer<const Self>                      ConstPointer;

		/** Method for creation through the object factory. */
		itkNewMacro(Self);

		/** Run-time type information (and related methods). */
		itkTypeMacro(PhaseSymmetryImageFilter, ImageToImageFilter);

		/** Some convenient typedefs. */
		typedef TInputImage                            InputImageType;
		typedef typename    InputImageType::Pointer    InputImagePointer;
		typedef typename    InputImageType::RegionType InputImageRegionType;
		typedef typename    InputImageType::PixelType  InputImagePixelType;

		typedef TOutputImage                              OutputImageType;
		typedef typename     OutputImageType::Pointer     OutputImagePointer;
		typedef typename     OutputImageType::RegionType  OutputImageRegionType;
		typedef typename     OutputImageType::PixelType   OutputImagePixelType;



		typedef OutputImagePixelType ComplexPixelComponentType;
		typedef OutputImagePixelType ImagePixelType;
		typedef std::complex< ComplexPixelComponentType >      ComplexPixelType;
		
		typedef itk::Image<ImagePixelType,TInputImage::ImageDimension> FloatImageType;




 itkSetMacro( MinWaveLength, double );
  itkGetConstMacro( MinWaveLength, double );

  itkSetMacro( Nscale, int );
  itkGetConstMacro( Nscale, int );

  itkSetMacro( SigmaAlpha, double );
  itkGetConstMacro( SigmaAlpha, double );
  
  itkSetMacro( SigmaFreq, double );
  itkGetConstMacro( SigmaFreq, double );

  itkSetMacro( Norient, int );
  itkGetConstMacro( Norient, int );

  itkSetMacro( Polarity, int );
  itkGetConstMacro( Polarity, int );

	itkSetMacro( Mult, double );
	itkGetConstMacro( Mult, double );

		itkSetMacro( PercentThreshold, double );
		itkGetConstMacro( PercentThreshold, double );

		void Initialize();
		/** Input and output images must be the same dimension, or the output's
		dimension must be one less than that of the input. */
#ifdef ITK_USE_CONCEPT_CHECKING
		/** Begin concept checking */
		itkConceptMacro(ImageDimensionCheck,
			(Concept::SameDimension<itkGetStaticConstMacro(InputImageDimension),
			itkGetStaticConstMacro(OutputImageDimension)>));
		/** End concept checking */
#endif


	protected:
		PhaseSymmetryImageFilter();
		virtual ~PhaseSymmetryImageFilter() {};
		void PrintSelf(std::ostream& os, Indent indent) const;

		/** Apply changes to the output image information. */
		virtual void GenerateOutputInformation();

		/** Apply changes to the input image requested region. */
		virtual void GenerateInputRequestedRegion();

		void GenerateData(void);

	private:
		PhaseSymmetryImageFilter(const Self&); //purposely not implemented
		void operator=(const Self&); //purposely not implemented

		//Comments below are from Ilkers Phase symmetry code
		double m_MinWaveLength; //Wavelength of smallest scale filter.
		double m_SigmaAlpha; // Ratio of angular interval between filter orientations and the standard deviation of the angular Gaussian function used to construct filters in the frequency plane.
		int m_Nscale;  //Number of wavelet scales. 
		int m_Norient; //Number of filter orientations.
		double m_SigmaFreq; //Ratio of the standard deviation of the Gaussian describing the log Gabor filter's transfer function in the frequency domain to the filter center frequency. 
		int m_Polarity; //Look for both black and white spots of symmetrry
		double m_Mult; //Scaling factor between successive filters.
		
		double m_PercentThreshold; //The percent threshold of the number of pixel values that are kept in the phase symmetry. Needs to between 0 and 1; if chose 0.5 then the top 50% of the values are kept and the rest of made zero.
		

		typename FloatImageType::Pointer m_PhaseSymmetry;




	};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPhaseSymmetryImageFilter.txx"
#endif

#endif

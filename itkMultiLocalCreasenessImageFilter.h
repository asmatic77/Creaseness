
#ifndef __itkMultiLocalCreasenessImageFilter_h
#define __itkMultiLocalCreasenessImageFilter_h

#include "itkRecursiveGaussianImageFilter.h"
#include "itkNthElementImageAdaptor.h"
#include "itkImage.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkPixelTraits.h"
#include "itkProgressAccumulator.h"

namespace itk
{

/** \class itkMultiLocalCreasenessImageFilter
 *
 * \brief Computes the multilocal creaseness measure [López99] 
 * 
 * 
 * \ingroup GradientFilters   
 * \ingroup Singlethreaded
 */
// NOTE that the ITK_TYPENAME macro has to be used here in lieu 
// of "typename" because VC++ doesn't like the typename keyword 
// on the defaults of template parameters
	template <typename TInputImage, typename TOutputImage = Image < 
		ITK_TYPENAME NumericTraits< ITK_TYPENAME TInputImage::PixelType>::RealType,
		::itk::GetImageDimension<TInputImage>::ImageDimension > > 
class ITK_EXPORT MultiLocalCreasenessImageFilter:
    public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef MultiLocalCreasenessImageFilter  Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;
    
  /** Pixel Type of the input image */
  typedef TInputImage                                    InputImageType;
  typedef typename TInputImage::PixelType                PixelType;
  typedef typename NumericTraits<PixelType>::RealType    RealType;

  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Define the image type for internal computations 
      RealType is usually 'double' in NumericTraits. 
      Here we prefer float in order to save memory.  */

  typedef double                                            InternalRealType;
  typedef Image<InternalRealType, 
                itkGetStaticConstMacro(ImageDimension) >   RealImageType;

  /** Type of the output Image */
  typedef TOutputImage                                      OutputImageType;
  typedef typename          OutputImageType::PixelType      OutputPixelType;
  typedef typename PixelTraits<OutputPixelType>::ValueType  OutputComponentType;
	typedef typename OutputImageType::Pointer OutputImagePointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** MultiLocalCreasenessImageFilter needs all of the input to produce an
   * output. Therefore, MultiLocalCreasenessImageFilter needs to provide
   * an implementation for GenerateInputRequestedRegion in order to inform
   * the pipeline execution model.
   * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion() throw(InvalidRequestedRegionError);

protected:
  MultiLocalCreasenessImageFilter();
  virtual ~MultiLocalCreasenessImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  /** Generate Data */
  void GenerateData( void );

  // Override since the filter produces the entire dataset
  void EnlargeOutputRequestedRegion(DataObject *output);

private:
  MultiLocalCreasenessImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
    
  OutputImagePointer                  m_outputImage;

  /** Normalize the image across scale space */
  bool m_NormalizeAcrossScale; 

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultiLocalCreasenessImageFilter.txx"
#endif

#endif

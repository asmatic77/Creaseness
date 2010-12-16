#include "itkMultiLocalCreasenessImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkImage.h"

int main(int argc, char* argv []  ) 
{
	
	typedef itk::Image<float,2> float2dImg;
	typedef itk::Image<float,3> float3dImg;
	typedef itk::Image<unsigned char, 2> uchar2dImg;

	typedef itk::ImageFileReader<float2dImg> ReaderFilter;
	ReaderFilter::Pointer reader = ReaderFilter::New();
	//reader->SetFileName("C:/Users/sergio.vera/Documents/MATLAB/PhD/distmap.png");
  reader->SetFileName("C:/Users/svera/Documents/MATLAB/PhD/distmap.png");
	try {
    /*
		typedef itk::RescaleIntensityImageFilter<float2dImg,float2dImg> RescaleFilter;
		RescaleFilter::Pointer rescale = RescaleFilter::New();
		rescale->SetInput(reader->GetOutput());
		rescale->SetOutputMinimum(0);
		rescale->SetOutputMaximum(1);
		typedef itk::MultiLocalCreasenessImageFilter<float2dImg,float2dImg> RidgeFilter;
		RidgeFilter::Pointer creases = RidgeFilter::New();
		creases->SetInput(rescale->GetOutput());
		creases->Update();
		typedef itk::ImageFileWriter<float2dImg> WriterFilter;
		WriterFilter::Pointer writer = WriterFilter::New();
		writer->SetFileName("ridges2d.mhd");
		writer->SetInput(creases->GetOutput());
		writer->Update();
*/
		
		typedef itk::ImageFileReader<float3dImg> ReaderFilter3;
		ReaderFilter3::Pointer reader3 = ReaderFilter3::New();
		//reader3->SetFileName("Z:/Investigación/Datos/Shared/Synthetic/sinteticvessel.mhd");
		reader3->SetFileName("arbol3d.mhd");

		typedef itk::InvertIntensityImageFilter<float3dImg,float3dImg> InvertFilter;
		InvertFilter::Pointer inverter = InvertFilter::New();
		inverter->SetMaximum(1.0);
		inverter->SetInput(reader3->GetOutput());

		typedef itk::DanielssonDistanceMapImageFilter<float3dImg,float3dImg> DistMap3;
		DistMap3::Pointer distmap3 = DistMap3::New();
		distmap3->SetInput(inverter->GetOutput());
		distmap3->SetUseImageSpacing(true);		
		distmap3->ReleaseDataFlagOn();		
		distmap3->Update();
		
		typedef itk::ImageFileWriter<float3dImg> WriterFilter3;
		WriterFilter3::Pointer writerDist = WriterFilter3::New();
		writerDist->SetFileName("Distmap.mhd");
		writerDist->SetInput(distmap3->GetOutput());
		writerDist->Update();
/*
		typedef itk::RescaleIntensityImageFilter<float3dImg,float3dImg> RescaleFilter3;
		RescaleFilter3::Pointer rescale3 = RescaleFilter3::New();
		rescale3->SetInput(distmap3->GetOutput());
		rescale3->SetOutputMinimum(0);
		rescale3->SetOutputMaximum(1);
 */
		typedef itk::MultiLocalCreasenessImageFilter<float3dImg,float3dImg> RidgeFilter3;
		RidgeFilter3::Pointer creases3 = RidgeFilter3::New();
		creases3->SetInput(distmap3->GetOutput());
		creases3->Update();
		
		WriterFilter3::Pointer writer3 = WriterFilter3::New();
		writer3->SetFileName("ridges3d.mhd");
		writer3->SetInput(creases3->GetOutput());
		writer3->Update();
		

	}catch(itk::ExceptionObject& e) {
		std::cout << e.what() << std::endl;
		return -1;
	}catch(...) {
		std::cout << "Problema\n";
		return -2;
	}
	return 0;
}
#include "tclap/CmdLine.h"
#include "itkMesh.h"
#include "itkQuadEdgeMesh.h"
#include "itkSTLMeshIOFactory.h"
#include "itkSTLMeshIO.h"
#include "itkMeshFileReader.h"
#include "itkMeshFileWriter.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
 
#include "itkCastImageFilter.h"
#include "itkTriangleMeshToBinaryImageFilter.h"
 

typedef class CmdLineType
{
public:
  std::string InIm, OutputIm;
  float buffer;
} CmdLineType;

void ParseCmdLine(int argc, char* argv[],
                  CmdLineType &CmdLineObj
                  )
{
  using namespace TCLAP;
  try
    {
    // Define the command line object.
    CmdLine cmd("loadMesh ", ' ', "0.9");

    ValueArg<std::string> inArg("i","input","input image",true,"result","string");
    cmd.add( inArg );
    ValueArg<std::string> outArg("o","output","output image", true,"","string");
    cmd.add( outArg );

    ValueArg<float> bArg("b","buffer","padding around bounding box", false, 5,"float");
    cmd.add(bArg);

    // Parse the args.
    cmd.parse( argc, argv );

    CmdLineObj.InIm = inArg.getValue();
    CmdLineObj.OutputIm = outArg.getValue();
    CmdLineObj.buffer = bArg.getValue();
    }
  catch (ArgException &e)  // catch any exc eptions
    {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}


////////////////////////////////////////////////////////
int main(int argc, char * argv[])
{

  const int dim = 3;
  CmdLineType CmdLineObj;
  ParseCmdLine(argc, argv, CmdLineObj);

  constexpr unsigned int Dimension = 3;
  using PixelType = float;

  using MeshType = itk::Mesh<PixelType, Dimension>;
 
  using QEMeshType = itk::QuadEdgeMesh<PixelType, Dimension>;

  itk::STLMeshIOFactory::RegisterOneFactory();

  using ReaderType = itk::MeshFileReader<MeshType>;
  using WriterType = itk::MeshFileWriter<MeshType>;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName(CmdLineObj.InIm);
  reader->Update();
  MeshType * mesh = reader->GetOutput();


  const MeshType::BoundingBoxType * BB = mesh->GetBoundingBox();
 
  std::cout << *BB << std::endl;
  std::cout << "bounds: " << BB->GetBounds() << std::endl;
  std::cout << "center: " << BB->GetCenter() << std::endl;
  using OutputPixelType = unsigned char;
  using OutputImageType = itk::Image<OutputPixelType, Dimension>;
 
  // bounding box characteristics to inform
  // image parameters
  float offset[3] = {0};
  float bbsize[3] = {0};
  // bit of trial and error to get the stl directions
  // sorted out. Bounding box is xpos1, xpos2, ypos1, ypos2, zpos1, zpos2.
  const float buffer = CmdLineObj.buffer;
  offset[0] = BB->GetBounds()[0] - buffer;
  offset[1] = BB->GetBounds()[2] - buffer;
  offset[2] = BB->GetBounds()[4] - buffer;

  bbsize[0] =  BB->GetBounds()[1] - BB->GetBounds()[0] + 2 * buffer;
  bbsize[1] =  BB->GetBounds()[3] - BB->GetBounds()[2] + 2 * buffer;
  bbsize[2] =  BB->GetBounds()[5] - BB->GetBounds()[4] + 2 * buffer;

  OutputImageType::SpacingType sp;
  sp.Fill(0.5);

  OutputImageType::SizeType sz;
  sz[0] = round(bbsize[0]/sp[0]);
  sz[1] = round(bbsize[1]/sp[1]);
  sz[2] = round(bbsize[2]/sp[2]);

 

  using FilterType = itk::TriangleMeshToBinaryImageFilter<MeshType, OutputImageType>;
  auto filter = FilterType::New();
  filter->SetInput(reader->GetOutput());
  filter->SetSize(sz);
  filter->SetSpacing(sp);
  filter->SetOrigin(offset);

  filter->SetInsideValue(itk::NumericTraits<OutputPixelType>::max());
  try
  {
    filter->Update();
  }
  catch (const itk::ExceptionObject & error)
  {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
  }
 
  try
  {
    itk::WriteImage(filter->GetOutput(), CmdLineObj.OutputIm);
  }
  catch (const itk::ExceptionObject & error)
  {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
  }
 
  return EXIT_SUCCESS;
}

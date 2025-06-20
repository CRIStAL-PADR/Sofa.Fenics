#include <sofa/fenics/forcefield/Material.h>
#include <sofa/helper/system/DynamicLibrary.h>
#include <SofaPython3/PythonEnvironment.h>
#include <sofa/helper/system/FileSystem.h>

namespace sofa::fenics
{

using sofa::helper::system::DynamicLibrary;

UfcxMaterial::UfcxMaterial() :
    d_filename(initData(&d_filename, std::string(""),  "filename",  "Fenics material description for the hyperelastic force field."))
   ,d_element(initData(&d_element, std::string(""), "element", "Type of element among [tetrahedron, hexahedron]."))
   ,d_element_order(initData(&d_element_order, 1, "element_order",  "Order of the element [1 or 2]."))
{
    d_filename.setGroup("YOLO");
}

typedef const char* (*f1)();

void UfcxMaterial::init()
{
    if(d_componentState.getValue() == sofa::core::objectmodel::ComponentState::Valid)
        return;

    auto filename = d_filename.getFullPath();
    auto extension = d_filename.getExtension();
    std::cout << "FILE is " << filename << " extensions " << extension << std::endl;

    std::string binaryMaterial="";
    if(extension == ".c" || sofa::helper::system::FileSystem::exists(filename))
    {
        std::cout << "Compiling material from " << filename << std::endl;

        std::cout << "Compiling file to binary " << filename << std::endl;

        sofapython3::PythonEnvironment::executePython([&filename,&binaryMaterial](){
            std::cout << "YO LO" << filename << std::endl;
            auto res = pybind11::eval("Sofa.livecoding.compile('"+filename+"')");
            binaryMaterial = pybind11::cast<std::string>(res);
        });
    }else{
        binaryMaterial = filename;
    }

    std::cout << "Loading file " << filename << std::endl;
    std::cout << "Loading binary file " << binaryMaterial << std::endl;

    auto handle = DynamicLibrary::load(binaryMaterial);

    ufcx_form* ufcxFormF = *reinterpret_cast<ufcx_form**>(DynamicLibrary::getSymbolAddress(handle, "form_SaintVenantKirchhoff_Tetra_F"));
    if( ufcxFormF->num_integrals(ufcx_integral_type::cell) != 1 )
    {
        msg_error() << "Only single integrales are supported right now in SOFA got " << ufcxFormF->num_integrals(ufcx_integral_type::cell) << msgendl
                    << "please contact author to request the feature";
        d_componentState = core::objectmodel::ComponentState::Invalid;
    }

    // Retrieve the constants from ufcx and map them to a data field
    auto data = new Data<std::string>("",true,false);
    data->setName("signature");
    data->setGroup("Ufcx");
    data->setDisplayed(true);
    data->setValue(ufcxFormF->signature);
    addData(data, data->getName());

    // Do the same for the field
    for(int i=0;i<ufcxFormF->num_coefficients;i++)
    {
        auto data = new Data<float>("",true,false);
        data->setName(ufcxFormF->coefficient_name_map()[i]);
        data->setGroup("Ufcx");
        data->setDisplayed(true);
        addData(data, data->getName());
    }

    for(int i=0;i<ufcxFormF->num_constants;i++)
    {
        auto data = new Data<float>("",true,false);
        data->setName(ufcxFormF->constant_name_map()[i]);
        data->setGroup("Ufcx Constants");
        data->setDisplayed(true);
        addData(data, data->getName());
    }

    // Retrieve
    ufcxComputeF = ufcxFormF->integrals(ufcx_integral_type::cell)[0];

    d_componentState = core::objectmodel::ComponentState::Valid;
}

const ufcx_integral* UfcxMaterial::getFunctionF(){ return ufcxComputeF; }
const ufcx_integral* UfcxMaterial::getFunctionK(){ return ufcxComputeJ; }

std::vector<double> UfcxMaterial::getParameters(){
    std::vector<double> tmp;
    return tmp;
}

}

#include <sofa/core/ObjectFactory.h>
namespace sofa::fenics
{

void registerUfcxMaterial(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects( sofa::core::ObjectRegistrationData("A Fenics based hyperelastic forcefield.")
        .add<UfcxMaterial>());
}

}

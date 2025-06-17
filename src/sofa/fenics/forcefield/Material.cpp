#include <sofa/fenics/forcefield/Material.h>

namespace sofa::fenics
{

UfcxMaterial::UfcxMaterial() :
    d_material_file(initData(&d_material_file, std::string(""),  "material_file",  "Fenics material description for the hyperelastic force field."))
   ,d_element(initData(&d_element, std::string(""), "element", "Type of element among [tetrahedron, hexahedron]."))
   ,d_element_order(initData(&d_element_order, 1, "element_order",  "Order of the element [1 or 2]."))
{

}

const ufcx_integral* UfcxMaterial::getFunctionF(){ return nullptr; }
const ufcx_integral* UfcxMaterial::getFunctionK(){ return nullptr; }

std::vector<double> UfcxMaterial::getParameters(){
    std::vector<double> tmp;
    return tmp;
}

}

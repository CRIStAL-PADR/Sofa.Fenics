#pragma once

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/DataFileName.h>

#include <sofa/fenics/config.h>
#include "sofa/fenics/ufcx.h"

namespace sofa::fenics
{

using sofa::core::objectmodel::Data;
using sofa::core::objectmodel::DataFileName;

class UfcxMaterial : public sofa::core::objectmodel::BaseObject
{
public:
    DataFileName d_material_file;
    Data<std::string> d_element;
    Data<int> d_element_order;

    const ufcx_integral* getFunctionF();
    const ufcx_integral* getFunctionK();

    std::vector<double> getParameters();

private:
    using Self      = sofa::fenics::UfcxMaterial;

    UfcxMaterial();
};

}

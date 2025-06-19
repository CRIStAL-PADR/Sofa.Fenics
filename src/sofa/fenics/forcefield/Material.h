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
    SOFA_CLASS(UfcxMaterial, BaseObject);

    DataFileName d_filename;
    Data<std::string> d_element;
    Data<int> d_element_order;

    UfcxMaterial();

    void init() override;

    // interface to retrieve the ufcx informations
    const ufcx_integral* getFunctionF();
    const ufcx_integral* getFunctionK();
    const ufcx_integral* getFunctionPi();
    std::vector<double> getParameters();

private:
    using Self      = sofa::fenics::UfcxMaterial;

    ufcx_integral* ufcxComputeF {nullptr};
    ufcx_integral* ufcxComputeJ {nullptr};
    ufcx_integral* ufcxComputePi {nullptr};
};

}

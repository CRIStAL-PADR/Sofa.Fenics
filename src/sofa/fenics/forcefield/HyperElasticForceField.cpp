#include <sofa/helper/AdvancedTimer.h>
#include "HyperElasticForceField.h"
#include "../ufcx.h"
namespace sofa::fenics
{

using Eigen::Matrix;

template<class DataTypes> HyperElasticForceField<DataTypes>::HyperElasticForceField() :
    l_topology(initLink("topology", "link to the container holding the topology"))
   ,l_material(initLink("material", "link to the UfcxMaterial to use as mechanical law"))
{
}

template<class DataTypes>
void HyperElasticForceField<DataTypes>::init(){
    if(!l_topology.get())
    {
        msg_error() << "Missing topology container";
        Inherit1::d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
    }
    if(!l_material.get())
    {
        msg_error() << "Missing material";
        Inherit1::d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
    }

    if(!l_material->isComponentStateValid())
    {
        msg_error() << "Invalid material";
        Inherit1::d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
    }

    std::cout << " SET COMPOENT TO VALID" << std::endl;
    Inherit1::d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
}

template<class DataTypes>
void HyperElasticForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* mparams, sofa::core::MultiVecDerivId fId)
{
    std::cout << "ADD FORCE... " << std::endl;
    if(!Inherit1::isComponentStateValid())
         return;

    std::cout << "ADD FORCE WITH COMPONENT STATE... " << std::endl;

    if (!Inherit1::mstate)
        return;

    std::cout << "ADDING A VALID " << std::endl;
}

template<class DataTypes>
SReal HyperElasticForceField<DataTypes>::getPotentialEnergy(
    const sofa::core::MechanicalParams* /* mparams */,
    const sofa::core::objectmodel::Data<VecCoord>& /* d_x */) const
{}

template<class DataTypes>
void HyperElasticForceField<DataTypes>::addForce(
        const sofa::core::MechanicalParams* mparams,
        sofa::core::objectmodel::Data<VecDeriv>& d_f,
        const sofa::core::objectmodel::Data<VecCoord>& d_x,
        const sofa::core::objectmodel::Data<VecDeriv>& d_v)
{

        using namespace sofa::core::objectmodel;
       using namespace sofa::helper;

       SOFA_UNUSED(mparams);
       SOFA_UNUSED(d_v);


       //if(d_material->MaterialIsAvailable()) {
       //    msg_error() << "No material named " << d_material->getMaterialName() << " available";
       //    return;
       //}

       constexpr int Dimension = 3;
       constexpr int NumberOfNodesPerElement = 4;

       ReadAccessor<Data<VecCoord>> sofa_x = d_x;
       ReadAccessor<Data<VecCoord>> sofa_x0 = this->mstate->readRestPositions();
       WriteAccessor<Data<VecDeriv>> sofa_f = d_f;
       const auto elements = l_topology->getTetrahedra();

       if (sofa_x.size() != sofa_f.size())
           return;

       if (sofa_x.empty() || elements.empty())
           return;

       const auto nb_nodes = sofa_x.size();

        Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>>    X       (sofa_x.ref().data()->data(),  nb_nodes, Dimension);
        Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>>    X0      (sofa_x0.ref().data()->data(), nb_nodes, Dimension);

       // Compute the displacement with respect to the rest position
       const auto u =  X - X0;

       // Get FEniCS F
       const ufcx_integral *computeF = l_material->getFunctionF();

       // Get constants from the l_material
       const std::vector<Real> constants = l_material->getParameters();

       sofa::helper::AdvancedTimer::stepBegin("HyperelasticForcefield_FEniCS::addForce");

       for (auto element : elements)
       {
           // Fetch the initial and current positions of the element's nodes
           Matrix<Real, NumberOfNodesPerElement, Dimension> current_nodes_position;
           Matrix<Real, NumberOfNodesPerElement, Dimension> coefficients;

           for (std::size_t i = 0; i < NumberOfNodesPerElement; ++i) {
               const auto t = element[i];
               current_nodes_position.row(i).noalias() = X0.row(t);
               coefficients.row(i).noalias() = u.row(t);
           }

           // Compute the nodal forces
           Matrix<Real, Eigen::Dynamic, NumberOfNodesPerElement, Dimension> nodal_forces;
           nodal_forces.fill(0);

           computeF->tabulate_tensor_float64(nodal_forces.data(), coefficients.data(), constants.data(), current_nodes_position.data(), nullptr, nullptr);

           for (size_t i = 0; i < NumberOfNodesPerElement; ++i) {
               const auto t = element[i];
               for (size_t j = 0; j < Dimension; ++j) {
                   sofa_f[t][j] -= nodal_forces(i,j);
               }
           }
       }

       sofa::helper::AdvancedTimer::stepEnd("HyperelasticForcefield_FEniCS::addForce");
}

template<class DataTypes>
void HyperElasticForceField<DataTypes>::addDForce(
        const sofa::core::MechanicalParams* /*mparams*/,
        sofa::core::objectmodel::Data<VecDeriv>& /*d_df*/,
        const sofa::core::objectmodel::Data<VecDeriv>& /*d_dx*/)
{}

template<class DataTypes>
void HyperElasticForceField<DataTypes>::addKToMatrix(BaseMatrix * /*matrix*/, SReal /*kFact*/, unsigned int & /*offset*/)
{}


template<class DataTypes>
void HyperElasticForceField<DataTypes>::assemble_stiffness(){

}

}


#include <sofa/core/ObjectFactory.h>
namespace sofa::fenics
{

void registerHyperElasticForceField(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects( sofa::core::ObjectRegistrationData("A Fenics based hyperelastic forcefield.")
        .add< sofa::fenics::HyperElasticForceField<defaulttype::Vec3Types> >(true));
}

}

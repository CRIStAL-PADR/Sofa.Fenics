#include "HyperElasticForceField.h"
#include "sofa/core/ComponentLibrary.h"
#include "sofa/helper/AdvancedTimer.h"
namespace sofa::fenics::forcefield
{

template<class DataTypes> HyperElasticForceField<DataTypes>::HyperElasticForceField() :
    d_material_file(initData(&d_material_file, std::string(""),  "material_file",  "Fenics material description for the hyperelastic force field."))
   ,d_element(initData(&d_element, std::string(""), "element", "Type of element among [tetrahedron, hexahedron]."))
   ,d_element_order(initData(&d_element_order, 1, "element_order",  "Order of the element [1 or 2]."))
{
}

template<class DataTypes>
void HyperElasticForceField<DataTypes>::init(){}

template<class DataTypes>
void HyperElasticForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* mparams, sofa::core::MultiVecDerivId fId)
{
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

       if (!this->mstate)
           return;

       //if(d_material->MaterialIsAvailable()) {
       //    msg_error() << "No material named " << d_material->getMaterialName() << " available";
       //    return;
       //}

       ReadAccessor<Data<VecCoord>> sofa_x = d_x;
       ReadAccessor<Data<VecCoord>> sofa_x0 = this->mstate->readRestPositions();
       WriteAccessor<Data<VecDeriv>> sofa_f = d_f;

       if (sofa_x.size() != sofa_f.size())
           return;
       const auto nb_nodes = sofa_x.size();
       const auto nb_elements = this->number_of_elements();

       if (nb_nodes == 0 || nb_elements == 0)
           return;


       Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>>    X       (sofa_x.ref().data()->data(),  nb_nodes, Dimension);
       Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Dimension, Eigen::RowMajor>>    X0      (sofa_x0.ref().data()->data(), nb_nodes, Dimension);

       // Compute the displacement with respect to the rest position
       const auto u =  X - X0;

       // Get FEniCS F
       const ufcx_integral *integral = material->FEniCS_F();


       // Get constants from the material
       const double constants_mooney[3] = {
                                       material->getMooneyRivlinConstants()(0, 0), // C01
                                       material->getMooneyRivlinConstants()(0, 1),  // C10
                                       material->getMooneyRivlinConstants()(0, 2)   // K
                                   };

       const double constants_ogden[9] = {
                                       material->getOgdenConstants()(0, 0), // bulk modulus
                                       material->getOgdenConstants()(0, 1), // a
                                       material->getOgdenConstants()(0, 2), // b
                                       material->getOgdenConstants()(0, 3), // a_f
                                       material->getOgdenConstants()(0, 4), // b_f
                                       material->getOgdenConstants()(0, 5), // a_s
                                       material->getOgdenConstants()(0, 6), // b_s
                                       material->getOgdenConstants()(0, 7), // a_fs
                                       material->getOgdenConstants()(0, 8), // b_fs
                                   };

       const double constants_else[2] = {
                                       material->getYoungModulusAndPoissonRatio()(0, 0), // Young Modulus
                                       material->getYoungModulusAndPoissonRatio()(0, 1)  // Poisson ratio
                                   };
       const double* constants;
       if(material->getMaterialName() == "MooneyRivlin") {
           constants = constants_mooney;
       } else if(material->getMaterialName() == "Ogden") {
           constants = constants_ogden;
       } else {
           constants = constants_else;
       }

       sofa::helper::AdvancedTimer::stepBegin("HyperelasticForcefield_FEniCS::addForce");

       for (std::size_t element_id = 0; element_id < nb_elements; ++element_id) {

           // Fetch the node indices of the element
           auto node_indices = this->topology()->domain()->element_indices(element_id);

           // Fetch the initial and current positions of the element's nodes
           Matrix<NumberOfNodesPerElement, Dimension> current_nodes_position;
           Matrix<NumberOfNodesPerElement, Dimension> coefficients;


           for (std::size_t i = 0; i < NumberOfNodesPerElement; ++i) {
               current_nodes_position.row(i).noalias() = X0.row(node_indices[i]);
               coefficients.row(i).noalias() = u.row(node_indices[i]);

           }

           // Compute the nodal forces
           Matrix<NumberOfNodesPerElement, Dimension> nodal_forces;
           nodal_forces.fill(0);

           integral->tabulate_tensor_float64(nodal_forces.data(), coefficients.data(), constants, current_nodes_position.data(), nullptr, nullptr);

           for (size_t i = 0; i < NumberOfNodesPerElement; ++i) {
               for (size_t j = 0; j < Dimension; ++j) {
                   sofa_f[node_indices[i]][j] -= nodal_forces(i,j);
               }
           }
       }

       sofa::helper::AdvancedTimer::stepEnd("HyperelasticForcefield_FEniCS::addForce");

       // This is the only I found to detect when a stiffness matrix reassembly is needed for calls to addDForce
       K_is_up_to_date = false;
       eigenvalues_are_up_to_date = false;


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
        .add< sofa::fenics::forcefield::HyperElasticForceField<defaulttype::Vec3Types> >(true));
}

}

#pragma once


#include <sofa/fenics/config.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/objectmodel/DataFileName.h>
#include <sofa/core/topology/BaseTopology.h>
#include <sofa/linearalgebra/BaseMatrix.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <functional>
#include <array>

#include <sofa/fenics/forcefield/Material.h>

namespace sofa::fenics
{

using sofa::linearalgebra::BaseMatrix;
using sofa::core::behavior::ForceField;
using sofa::core::objectmodel::Data;
using sofa::core::objectmodel::DataFileName;

template<class DataTypes>
class HyperElasticForceField : public ForceField<DataTypes>
{
private:
    using Self      = sofa::fenics::HyperElasticForceField<DataTypes>;
    using VecCoord  = typename DataTypes::VecCoord;
    using VecDeriv  = typename DataTypes::VecDeriv;
    using Coord     = typename DataTypes::Coord;
    using Deriv     = typename DataTypes::Deriv;
    using Real      = typename DataTypes::Real;

    template <typename ObjectType>
    using Link = sofa::core::objectmodel::SingleLink<Self, ObjectType, sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK>;

public:
    SOFA_CLASS(SOFA_TEMPLATE(HyperElasticForceField, DataTypes), SOFA_TEMPLATE(ForceField, DataTypes));
    
    HyperElasticForceField();
    
    void init() override;    
    void addForce(const sofa::core::MechanicalParams* mparams, sofa::core::MultiVecDerivId fId ) override;

    void addForce(
        const sofa::core::MechanicalParams* mparams,
        sofa::core::objectmodel::Data<VecDeriv>& d_f,
        const sofa::core::objectmodel::Data<VecCoord>& d_x,
        const sofa::core::objectmodel::Data<VecDeriv>& d_v) override;

    void addDForce(
        const sofa::core::MechanicalParams* /*mparams*/,
        sofa::core::objectmodel::Data<VecDeriv>& /*d_df*/,
        const sofa::core::objectmodel::Data<VecDeriv>& /*d_dx*/) override;

    SReal getPotentialEnergy(
        const sofa::core::MechanicalParams* /* mparams */,
        const sofa::core::objectmodel::Data<VecCoord>& /* d_x */) const override;

    void addKToMatrix(BaseMatrix * /*matrix*/, SReal /*kFact*/, unsigned int & /*offset*/) override;

    /**
     *  Assemble the stiffness matrix K.
     *
     *  Since the stiffness matrix is function of the position vector x, i.e. K(x), this method will
     *  use the mechanical state vector used in the last call to addForce (usually Position or FreePosition).
     *  If another state vector should be used as the x, use instead the update_stiffness(x) method.
     *
     *  A reference to the assembled stiffness matrix K as a column major sparse matrix can be later
     *  obtained using the method K().
     *
     */
    void assemble_stiffness();

    /**
     *  Assemble the stiffness matrix K.
     *
     *  Since the stiffness matrix is function of the position vector x, i.e. K(x), this method will
     *  use the data vector x passed as parameter. If the
     *
     *  A reference to the assembled stiffness matrix K as a column major sparse matrix can be later
     *  obtained using the method K().
     *
     */
    //virtual void assemble_stiffness(const sofa::core::objectmodel::Data<VecCoord> & x);

    /**
     *  Assemble the stiffness matrix K.
     *
     *  Since the stiffness matrix is function of the position vector x, i.e. K(x), this method will
     *  use the position vector x passed as a Eigen matrix nx3 parameter with n the number of nodes.
     *
     *  A reference to the assembled stiffness matrix K as a column major sparse matrix can be later
     *  obtained using the method K().
     *
     */
     //template <typename Derived>
    //void assemble_stiffness(const Eigen::MatrixBase<Derived> & x);

    //template <typename Derived>
    //void assemble_stiffness(const Eigen::MatrixBase<Derived> & x, const Eigen::MatrixBase<Derived> & x0);

    Link<UfcxMaterial> l_material;
    Link<core::topology::TopologyContainer> l_topology;

private:

    // These private methods are implemented but can be overridden

};

}

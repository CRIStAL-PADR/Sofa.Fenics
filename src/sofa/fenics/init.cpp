#include <sofa/core/ObjectFactory.h>
#include <sofa/fenics/config.h>

namespace sofa::fenics
{

    void registerHyperElasticForceField(sofa::core::ObjectFactory* factory);
}

extern "C" {
    SOFA_FENICS_API void initExternalModule();
    SOFA_FENICS_API const char* getModuleName();
    SOFA_FENICS_API const char* getModuleVersion();
    SOFA_FENICS_API const char* getModuleLicense();
    SOFA_FENICS_API const char* getModuleDescription();
    SOFA_FENICS_API void registerObjects(sofa::core::ObjectFactory* factory);
}


void initExternalModule()
{
    static bool first = true;
    if (first)
    {
        first = false;
    }
}

const char* getModuleName()
{
    return sofa_tostring(SOFA_TARGET);
}

const char* getModuleLicense()
{
    return "";
}

const char* getModuleVersion()
{
    return sofa_tostring(SOFA_FENICS_VERSION);
}

const char* getModuleDescription()
{
    return "Use Fenics to define mechanical laws";
}

void registerObjects(sofa::core::ObjectFactory* factory)
{
    sofa::fenics::registerHyperElasticForceField(factory);
}

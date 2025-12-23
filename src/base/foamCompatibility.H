#ifndef foamCompatibility_H
#define foamCompatibility_H

// --- Field Access Macros ---
// solids4foam uses macros like these to handle the .boundaryFieldRef() change
#if defined(OPENFOAM_COM) || defined(OPENFOAM_ORG)
    #define getBoundaryField(field) field.boundaryFieldRef()
    #define getInternalField(field) field.primitiveFieldRef()
#elif defined(FOAMEXTEND)
    #define getBoundaryField(field) field.boundaryField()
    #define getInternalField(field) field.internalField()
#endif

// --- Dictionary Access ---
// Handling the shift from .lookup() to .get<T>()
#if defined(FOAMEXTEND)
    #define dictGet(dict, type, name) dict.lookupType<type>(name)
#else
    #define dictGet(dict, type, name) dict.get<type>(name)
#endif

// --- Topology Changes or Meshing ---
// foam-extend and OpenFOAM.com have different class names for mesh subsetting
#if defined(FOAMEXTEND)
    #define fvMeshSubsetType subMesh
#else
    #define fvMeshSubsetType fvMeshSubset
#endif

// --- Version Specific Code Blocks ---
// Template for when to execute entirely different logic
//#if defined(FOAMEXTEND)
//    #define IF_EXTEND(code) code
//    #define IF_NOT_EXTEND(code)
//#else
//    #define IF_EXTEND(code)
//    #define IF_NOT_EXTEND(code) code
//#endif

#endif

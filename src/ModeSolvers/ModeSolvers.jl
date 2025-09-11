module ModeSolvers

using Interpolations
using Trapz
using LinearAlgebra
using SparseArrays
using Bessels
using Gridap
using GridapGmsh
using Gridap.Fields
using Gridap.Geometry
using Gridap.TensorValues
using Gridap.CellData
using Gridap.Fields
using Gridap.ReferenceFEs
using HCubature
using StaticArrays
using NearestNeighbors
using Roots
using Unitful

import OpticalFibers.PhysicalData: c,mu0,eps0,h,Z0,c_Unitful,mu0_Unitful,eps0_Unitful,h_Unitful,Z0_Unitful
import OpticalFibers: eigs_MUMPS,eigs_LU,eigs_CUDA,get_companion,add_cylindrical_PML,add_twist_PML,add_rectangular_PML,nb_args,UnitfulField,FunctionField,FEMField,ArrayField,UnitfulModel
import Base: *, /, -, +, sin, cos, sqrt, abs, abs2, real, imag, conj

export FiberEMField
export ScalarFiberEMField
export ScalarFiberEMField1D
export ScalarFiberEMField2D
export VectorFiberEMField
export Mode
export getEMField
export nonLinearCoefficient
export multi_step_fiber_modes
export MFD
export Aeff
export normalize!
export normalize
export overlap
export FEM1D
export FEM2D
export FEM2D_anisotropic
export FEM2D_periodic
export convertTo2D
export convertToVector

include("Mode_Field.jl")
include("multi_step_fiber.jl");
include("FEM.jl");

end # module

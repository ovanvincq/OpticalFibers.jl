module ModeSolvers

using Interpolations
using Trapz
using LinearAlgebra
using SparseArrays
using Bessels
using Gridap
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

import OpticalFibers.PhysicalData: c,mu0,eps0,h,Z0
import OpticalFibers: eigs_MUMPS,eigs_LU,tensor3,inverse,get_companion,add_cylindrical_PML,add_rectangular_PML,nb_args, integrate1D,integrate2D,tensorComponent,isaNumber,isaFunction

export Field
export ScalarField1D
export ScalarField2D
export VectorField2D
export ScalarFieldFEM1D
export ScalarFieldFEM2D
export VectorFieldFEM2D
export ScalarFieldMatrix1D
export ScalarFieldMatrix2D
export VectorFieldMatrix2D
export ScalarFieldFunction1D
export ScalarFieldFunction2D
export VectorFieldFunction2D
export Mode
export PoyntingVector
export normalize!
export overlap
export nonLinearCoefficient
export multi_step_fiber_modes
export MFD
export Aeff
export FEM1D
export FEM2D
export FEM2D_anisotropic
export FEM2D_periodic
export getField
export computeField
export losses
export isvalidField
export isvalidMode
export getValue
export getField
export convertToMatrix
export convertTo2D
export convertToVector

include("Mode_Field.jl")
include("multi_step_fiber.jl");
include("FEM.jl");

end # module

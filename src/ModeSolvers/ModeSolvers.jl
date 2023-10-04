module ModeSolvers

using Interpolations;
using Trapz;
using LinearAlgebra;
using SparseArrays;
using Arpack;
using SpecialFunctions;
using Gridap;
using Gridap.Fields
using Gridap.Geometry
using Gridap.TensorValues
using Gridap.CellData

import OpticalFibers.PhysicalData: c,mu0,eps0,h,Z0
import OpticalFibers: indexGaussSampling1D,Laplacian,indexGaussSampling2D,eigs_MUMPS,eigs_LU,tensor3,inverse,get_companion

export ScalarMode1D
export ScalarMode2D
export ScalarModeFEM
export VectorModeFEM
export VectorMode
export ScalarField
export VectorField
export ScalarFieldFEM
export VectorFieldFEM
export PoyntingVector
export normalize!
export overlap
export FD
export nonLinearCoefficient
export multi_step_fiber_modes
export MFD
export Aeff
export Interpolation
export FEM

include("Mode_Field.jl")
include("FD.jl")
include("multi_step_fiber.jl");
include("FEM.jl");

end # module

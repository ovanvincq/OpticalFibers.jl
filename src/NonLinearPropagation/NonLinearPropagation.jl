module NonLinearPropagation

export GNLSE_param
export GNLSE_mode
export GNLSE_Raman
export propagation

import OpticalFibers.PhysicalData: c,mu0,eps0,h,Z0

include("structs_GNLSE.jl")
include("Propagation_CPU.jl")

end
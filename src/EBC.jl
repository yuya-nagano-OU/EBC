using WignerD
using SparseArrays
using LinearAlgebra
using StaticArrays
using Healpix


include("./function/convolver.jl")
include("./function/mapmaker.jl")
include("./function/rotaters.jl")
include("./function/structures.jl")
include("./function/tools.jl")
include("./function/sph_functions.jl")



export simple_rotater, swignerd_calc
export gen_ConvolutionParams
export alm_idx, for_healpy_order, unique_theta_num, unique_theta_val, store_wignerd


# rotator
include("./function2/rotator.jl")
export WignerD_calculator!, WignerD_calculator_fast, global_wignerD, local_effective_wignerD 
include("./function2/wignerd.jl")
export wignerd_mdown_from_l

#sphcoeff
include("./function2/sphcoeff.jl")
export ConvolutionSky, ConvolutionBeam
export lmr_idx, alm_idx

#convolution
include("./function2/convolution.jl")
export ConvolutionCalculate



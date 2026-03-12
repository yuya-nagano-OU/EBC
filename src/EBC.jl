using WignerD
using SparseArrays
using LinearAlgebra
using StaticArrays
using Healpix


include("./function/convolution.jl")
include("./function/sphcoeff.jl")
include("./function/rotator.jl")
include("./function/tools.jl")

#sphcoeff.jl
export ConvolutionSky, ConvolutionBeam
export alm_idx
export slice_spin_alm_by_l
export slice_spin_blm_by_l

#convolution.jl
export ConvolutionCalculate
export convolver_1pixel, compute_pixel_convolution

#rotator.jl
export local_effective_wignerD_conj_reduced_formapmake
export global_wignerD_conj
export calc_local_euiler_angles
export global_wignerd, global_d2D_conj

#tools.jl
export ring_pixel_range, pixels_in_ring, first_pixel_in_ring, unique_theta_val
export build_index_dict

#=
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
export local_effective_wignerD_conj, local_effective_wignerD_conj_reduced
export calc_local_euiler_angles
include("./function2/wignerd.jl")
export wignerd_mdown_from_l

#sphcoeff
include("./function2/sphcoeff.jl")
export ConvolutionSky, ConvolutionBeam
export lmr_idx, alm_idx

#convolution
include("./function2/convolution.jl")
export ConvolutionCalculate
=#


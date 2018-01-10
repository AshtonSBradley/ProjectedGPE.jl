"""

`kernel = scatteringkernel(k,σ)`

Construct the scattering kernel in momentum space and returns on the same `k`-grid.
This is used to construct both the energy-damping term, and noise in the energy-damped
stochastic projected Gross-Pitaevskii equation.

`k` the momentum-space grid in where scattering is to be evaluated.

`σ` width parameter for `k`-space degrees of freedom in the case of N=1,2 spatial dimensions.

Low dimensional SPGPE scattering kernel functions:
[Low-D stochastic projected GPE](https://arxiv.org/abs/1507.02023)

"""
function scatteringkernel(k,σ=1)
if ndims(k) == 1
  return erfcx.(abs.(k)*σ/sqrt(2))/sqrt(8*π*σ^2)
elseif ndims(k) == 2
  return exp.(abs.(k*σ/2)^2)*besselk(0,abs.(k*σ/2)^2)/(2*π)
elseif ndims(k) == 3
 return 1./abs.(k)
end
end

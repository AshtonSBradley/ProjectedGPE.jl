"""

`S(k,σ)=scatteringkernel(k,σ)`
Evaluates the scattering kernel in momentum space.
This is used to construct both the energy-damping term, and noise in the energy-damped
stochastic projected Gross-Pitaevskii equation.

## Arguments
* `k` - the momentum-space grid in where scattering is to be evaluated
* `σ` - width parameter for `k`-space degrees of freedom in the case of N=1,2 spatial dimensions.
"""
scatteringkernel(k,σ=1)
if ndims(k) == 1
  return erfcx(abs(k)*σ/sqrt(2))/sqrt(8*π*σ^2)
elseif ndims(k) == 2
  return exp((k*σ/2)^2)*besselk(0,(k*σ/2)^2)/(2*π)
else
 return 1./abs(k)
end

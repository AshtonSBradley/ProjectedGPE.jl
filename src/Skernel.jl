"""
    kernel = Skernel(k,σ=1.0)

Construct the scattering kernel functions in momentum space.
The kernel function is used to construct both the energy-damping term, and noise in the energy-damped
stochastic projected Gross-Pitaevskii equation.

`k`: momentum-space grid in where scattering is to be evaluated.

`σ`: width parameter for `k`-space degrees of freedom in the case of N=1,2 spatial dimensions.

External links

[Low-Dimensional Stochastic Projected Gross-Pitaevskii Equation, Bradley, Rooney, MacDonald, Physical Review A 92, 033631 (2015)](https://arxiv.org/abs/1507.02023)

"""
function Skernel1(k,σ=1.0)
    erfcx(abs(k)*σ/sqrt(2))/sqrt(8*π*σ^2)
end
function Skernel2(k,σ=1.0)
    k==0.0 ? 0.0 : exp(abs(k*σ/2)^2)*besselk(0,abs(k*σ/2)^2)/(2*π)
end
function Skernel3(k)
    k==0.0 ? 0.0 : 1/abs(k)
end

Skernel(k::Array{Float64,1},σ=1.0) = Skernel1.(k,σ)
Skernel(k::Array{Float64,2},σ=1.0) = Skernel2.(k,σ)
Skernel(k::Array{Float64,3},σ=1.0) = Skernel3.(k)

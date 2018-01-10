"""
    T = eigmat(basis,M,x,ω=1.0,α=1.0)

Create a matrix of orthonormal mode functions for the chosen basis.
The matrix is in a form that allows the transformation
from the mode coefficients to the spatial grid `x`:

```math
\\psi(x) = \\sum_j c_j\\phi_j(x) \\equiv Tc
```

where `c` is a column vector of coefficients in the basis.

At present `basis = "Hermite"` is implemented.

`basis`: set of eigenfunctions representing the c-field.

`M`: number of modes in the basis.

`x`: spatial grid to which the coefficients are mapped.

`ω`: mode frequency, in units of the chosen reference frequency.

`α`: extra input for the `laguerre` basis.

Defaults of the last two arguments are 1.0 and 0.0 respectively.
"""

function eigmat(basis,M,x,ω=1.0,α=0.0)
  if basis=="Hermite"
    M > 371 && error("Quadrature does not converge for M > 371.")
    x = convert(Vector{BigFloat},x)
    ψ0 = exp.(-(√ω*x).^2/2)*BigFloat((ω/π)^(1/4))
    ψ1 = BigFloat(sqrt(2))*exp.(-(√ω*x).^2/2).*(√ω*x)*BigFloat((ω/π)^(1/4))
    T=zeros(x*ones(1,M))
    n=convert(Vector{BigFloat},collect(0:M-1))
    T[:,1] = ψ0
    T[:,2] = ψ1
    for m=1:M-2
      T[:,m+2]=sqrt(2/(n[m+2]))*(√ω*x).*T[:,m+1]-sqrt(n[m+1]/n[m+2])T[:,m]
    end
    x = convert(Vector{Float64},x)
    T = convert(Matrix{Float64},T)

  elseif basis=="Laguerre"
    error(basis," basis not implemented")
    #M>26 && error("recursion unstable for M >=27")
    #alternate approach using ApproxFun.jl -same instability as
    #previous approach in logspace
    #=
    T=zeros(x*ones(1,M))
    n=0:BigFloat(M)-1;
    Lnorm = exp.(lgamma.(n+α+1)-lgamma.(n+1))/ω

    for j = 0:M-1
      c=zeros(M)
      c[j+1]=1.
      Lna=Fun(Laguerre(α),c./sqrt(Lnorm));
      T[:,j+1]=Lna.(ω*x).*exp(-ω*x/2).*(ω*x).^(α/2)
    end
    =#
  else
    error(basis," basis not implemented")
  end
  return T
end

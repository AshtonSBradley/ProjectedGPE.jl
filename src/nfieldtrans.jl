"""

x,w,T = nfieldtrans(basis,M,n,ω=1)

Constructs transforms and associated arrays for precise numerical quadrature evaluation of n-field integrals,
starting from a representation of the quantum state with respect to a particular basis of eigenstates.

## Arguments
 - `n` order of the field product.
 - `M` number of modes in the c-field.
 - `basis` string argument denoting the basis of eigenstates representing c-field state.
 - `T` is the linear transformation matrix that affects the mapping
 - `ω` frequency of the oscillator states in given direction, relative to chosen reference frequency.
`T*c` = `` ψ(x)≡ ∑_{j=1}^{M}c_jϕ_j(x)``

for a state represented by `M` coefficients, the number of modes in the c-field.
 - `x` is the quadrature grid onto which ``ψ(x)`` is mapped
 - `w` are weights such that the an exact integral may be carreid out.
The integral must be a product of order `n` in the field ``ψ``.

The integrals is performed by
1. Transforming to the quadrature grid using `T`
2. Constructing the product and then evaluating the sum: `∑ⱼwⱼ*ψ(xⱼ)ⁿ=sum(w*ψ^n)`

# Examples

## C-field population
Compute the number of particles in the C-field, for a state of `M` modes.

```
julia> M = 30;
c = randn(M)+im*randn(M);
x,w,T=nfieldtrans("Hermite",M,2);
ψ = T*c;
N = sum(w.*abs(ψ).^2)
73.24196674113007
```

Computes the integral ``N = ∫dx |ψ(x)|^2`` as may be checked by direct summation:

```
julia> sum(abs(c).^2)
73.24196675017353
```

 Relative error of numerical integrals will usually be smaller than ~1e-10.

## Interaction energy
The most common integral of this type involves a four-field product, `n=4`.
The PGPE interaction energy is of this form, as is the nonlinear term in the PGPE;
the exact propagation of the PGPE requires repeated use of this 4-field transformation:

```
julia> x,w,T=nfieldtrans("Hermite",M,4);
ψ = T*c;
Uint = sum(w.*abs(ψ).^4)
552.9762736751692
```
computes the integral ``U_{\int}≡∫ dx|ψ|^4`` to accuracy very close to working precision.

**Warning:** Using the transform `T` for physical analysis in position space should be avoided
as the `x` grid is non-uniformly spaced.
Instead, use `eigmat.jl` to create a transform to a specific position grid.
"""

function nfieldtrans(basis,M,K,ω=1.,α=0.)
    iseven(K*M) ? n=Int(K*M/2) : n=Int((K*M+1)/2)
    if basis=="Hermite"
    x, w = gausshermite(n)
    w    = w.*exp.(x.^2)/sqrt(K*ω/2)
    T    = eigmat("Hermite",M,x/sqrt(K*ω/2),ω)
    elseif basis=="Laguerre"
      error(basis,"basis not implemented yet.")
      #=
      #sort out stable recursion for Laguerre
      M>26 && error("recursion unstable for M >=27")
      #needs proper testing
      x, w = gausslaguerre(n)
      w    = w.*exp(x)/(K*ω/2)./x.^α #careful wieght check needed for order K
      T    = eigmat("Laguerre",M,x/(K*ω/2),ω,α)
      =#
    else
        error(basis," basis not implemeted yet.")
    end
    return x,w,T
end

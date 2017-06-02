"""

x,w,T = nfieldtrans(basis,M,n,ω)

Constructs transforms and associated arrays for precise numerical quadrature evaluation of n-field integrals,
starting from a representation of the quantum state with respect to a particular basis of eigenstates.

## Arguments
 - `n` order of the field product.
 - `M` number of modes in the c-field.
 - `basis` is an optional string argument denoting the basis of eigenstates representing c-field state. Default is "hermite"
 - `T` is the linear transformation matrix that affects the mapping

`T*c` = `` ψ(x)≡ ∑_{j=1}^{M}c_jϕ_j(x)``

for a state represented by `M` coefficients, the number of modes in the c-field.
 - `x` is the quadrature grid onto which ``ψ(x)`` is mapped
 - `w` are weights such that the an exact integral may be carreid out.
The integral must be a product of order `n` in the field ``ψ``, and it is assumed that `n` is even.

The integrals is performed by
1. Transforming to the quadrature grid using `T`
2. Constructing the product and then evaluating the sum: `∑ⱼwⱼ*ψ(xⱼ)ⁿ=sum(w*ψ^n)`

# Examples

## C-field population
Compute the number of particles in the C-field, for a state of `M` modes.

```
julia> M = 30;
c = randn(M)+im*randn(M);
x,w,T=nfieldtrans(M,2);
ψ = T*c;
N = sum(w.*abs(ψ).^2)
73.24196674113007
```

Compuates the integral ``N = ∫dx |ψ(x)|^2`` as may be checked by direct summation:

```
julia> sum(abs(c).^2)
73.24196675017353
```

The accuracy seen in this example (10 digits) is a worst-case scenario - a random superposition of all modes is a high temperature limit.
Accuracy will normally approach 15 digits.

## Interaction energy
The most common integral of this type involves a four-field product, `n=4`.
The PGPE interaction energy is of this form, as is the nonlinear term in the PGPE;
the exact propagation of the PGPE requires repeated use of this 4-field transformation:

```
julia> x,w,T=nfieldtrans(4,M);
ψ = T*c;
Uint = sum(w.*abs(ψ).^4)
552.9762736751692
```
computes the integral ``U_{\int}≡∫ dx|ψ|^4`` to accuracy very close to working precision.

**Warning:** Using the transform `T` for physical analysis in position space should be avoided
as the `x` grid is non-uniformly spaced.
Instead, use `eigmat.jl` to create a transform to a specific position grid.
"""

function nfieldtrans(basis,M,n,ω=1)
    iseven(M)  ? nothing : error("M must be an even integer ")
    if basis=="hermite"
    x, w = gausshermite(n*M/2)
    w    = w.*exp(x.^2)/sqrt(n*ω/2)
    T    = eigmat("hermite",M,x/sqrt(n*ω/2),n*ω/2)
    return x,w,T
    else
        error(basis," basis not implemeted yet.")
    end
end

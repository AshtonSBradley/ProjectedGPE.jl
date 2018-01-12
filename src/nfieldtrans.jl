"""
    x,w,T = nfieldtrans(M,ω=1.0;n=4,basis="Hermite",α=0.0)

Construct transforms and associated arrays for numerical quadrature evaluation of n-field integrals,
starting from a representation of the quantum state with respect to a particular basis of eigenstates.

### Arguments
`basis`: string argument denoting the basis of eigenstates representing c-field state.

`M`: number of modes in the c-field.

`ω`: frequency of the oscillator states in given direction, relative to chosen reference frequency.

`n`: order of the field product.

`basis`: name of orthonormal eigenstate basis. Currently `"Hermite"` is implemented.

`α`: extra parameter for Laguerre basis.

### Outputs
`T`: linear transformation matrix that takes spectral state coefficienets to a quadrature grid. In 1D:

```math
Tc = \\psi(x) \\equiv \\sum_{j=1}^{M}c_j\\phi_j(x)
```

for a state represented by `M` coefficients, the number of modes in the c-field.

 * `x` is the quadrature grid onto which ``\\psi(x)`` is mapped
 * `w` are weights such that the an exact integral may be carreid out.

The integral must be a product of order `n` in the field ``\\psi``.

The integral is performed by

1. Transforming to the quadrature grid using `T`

2. Constructing the product and then evaluating the sum:
```math
\\sum_jw_j*\\psi(x_j)^n = sum(w.*\\psi.^n)
```

## Examples

### C-field population
Compute the number of particles in the C-field, for a state of `M` modes:

```julia
M = 30
c = randn(M)+im*randn(M)
x,w,T = nfieldtrans(M,2)
ψ = T*c
N = sum(w.*abs2.(ψ))
73.24196674113007
```

Computes the integral ``N = \\int dx |\\psi(x)|^2`` as may be checked by direct summation:

```julia
sum(abs2.(c))
73.24196675017353
```

Relative error of numerical integrals will usually be smaller than ~1e-10.

### Interaction energy
The most common integral of this type involves a four-field product, `n=4`.
The PGPE interaction energy is of this form, as is the nonlinear term in the PGPE;
the exact propagation of the PGPE requires repeated use of this 4-field transformation:

```julia
x,w,T = nfieldtrans(M,4)
ψ = T*c
Uint = sum(w.*abs.(ψ).^4)
552.9762736751692
```
computes the integral ``U_{\\rm int}\\equiv \\int dx|\\psi(x)|^4`` to accuracy very close to working precision.

"""

function nfieldtrans(M,ω=1.0;n=4,basis="Hermite",α=0.0)
    iseven(n*M) ? K=Int(n*M/2) : K=Int((n*M+1)/2)
    if basis=="Hermite"
    x, w = gausshermite(K)
    w    = @. exp(log(w)+x^2)/sqrt(n*ω/2)
    T    = eigmat(M,x/sqrt(n*ω/2),ω=ω,basis=basis)
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

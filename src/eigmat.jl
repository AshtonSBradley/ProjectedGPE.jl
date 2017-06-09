"""

T    = eigmat(basis,M,x,ω,α)

Returns a matrix of mode functions of the chosen basis.
The matrix is in a form that allows the transformation
from the mode coefficients to the spatial grid `x`:

`ψ(x)=∑ⱼ cⱼϕⱼ(x) ≡ T*c`

where `c` is a column vector of coefficients in the basis.

At present `basis = "Hermite"` is implemented.

## Arguments
-  `basis` is the set of eigenfunctions representing the c-field
- `M` is the number of modes in the spatial direction denoted by
- `x`, the spatial grid to which the coefficients are mapped.
- `ω` is the *relative* frequency, in units of the chosen reference frequency.
- `α` is an extra input for the `laguerre` basis.
Defaults of the last two arguments are 1.

"""

function eigmat(basis,M,x,ω=1.0,α=1.0)
  if basis=="Hermite"
    x,w=gausshermite(M)
    T=zeros(x*ones(1,M))
    n=0:BigFloat(M)-1;
    Hnorm = sqrt(BigFloat(π)/ω)*BigFloat(2).^n.*factorial.(n)

    for j = 0:M-1
    c=zeros(M)
    c[j+1]=1.
    Hn=Fun(Hermite(),c./sqrt(Hnorm));
    T[:,j+1]=Hn.(sqrt(ω)*x).*exp(-ω*x.^2/2)
    end

  elseif basis=="Laguerre"
    x,w=gausslaguerre(M)
    T=zeros(x*ones(1,M))
    n=0:BigFloat(M)-1;
    Hnorm = sqrt(BigFloat(π)/ω)*BigFloat(2).^n.*factorial.(n)

    for j = 0:M-1
    c=zeros(M)
    c[j+1]=1.
    Hn=Fun(Hermite(),c./sqrt(Hnorm));
    T[:,j+1]=Hn.(sqrt(ω)*x).*exp(-ω*x.^2/2)
    end
  else
    error("basis not implemented")
  end
  return T
end

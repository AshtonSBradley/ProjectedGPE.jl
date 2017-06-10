"""

T    = eigmat(basis,M,x,ω,α)

Returns a matrix of mode functions of the chosen basis.
The matrix is in a form that allows the transformation
from the mode coefficients to the spatial grid `x`:

`ψ(x)=∑ⱼ cⱼϕⱼ(x) ≡ T*c`

where `c` is a column vector of coefficients in the basis.

At present `basis = "Hermite"` is implemented.

## Arguments
- `basis` is the set of eigenfunctions representing the c-field
- `M` is the number of modes in the spatial direction denoted by
- `x`, the spatial grid to which the coefficients are mapped.
- `ω` is the *relative* frequency, in units of the chosen reference frequency.
- `α` is an extra input for the `laguerre` basis.
Defaults of the last two arguments are 1. and 0. respectively.

"""

function eigmat(basis,M,x,ω=1.0,α=0.0)
  if basis=="Hermite"
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
    M>26 && error("recursion unstable for M >=27")
    #alternate approach using ApproxFun.jl -same instability as previous approach 
    T=zeros(x*ones(1,M))
    n=0:BigFloat(M)-1;
    Lnorm = exp.(lgamma.(n+α+1)-lgamma.(n+1))/ω

    for j = 0:M-1
      c=zeros(M)
      c[j+1]=1.
      Lna=Fun(Laguerre(α),c./sqrt(Lnorm));
      T[:,j+1]=Lna.(sqrt(ω)*x).*exp(-sqrt(ω)*x/2).*(sqrt(ω)*x).^(α/2)
    end
    #==
    #old approach uses log-space but is also unstable for M>=27
    T=zeros(x*ones(1,M))
    x=complex(x,0)
  	laguerre_nm1 = real(exp((α*log(x)-lgamma(α+1))/2))
	  norm12 = exp((+lgamma(1+α)-lgamma(α+2))/2)
  	laguerre_n = norm12*(1+α-sqrt(ω)*x).*laguerre_nm1.*exp(-sqrt(ω)*x/(2*M))

  	T[:,1] = laguerre_nm1.*exp(0*sqrt(ω)*x/(2*M)-sqrt(ω)*x/2)
  	T[:,2] = laguerre_n.*exp(1*sqrt(ω)*x/(2*M)-sqrt(ω)*x/2)

	if M > 2
  		for nn = 1:M-2
			norm12 = exp((log(nn+1)+lgamma(nn+1+α)-lgamma(nn+α+2))/2)
			norm13 = exp((log((nn+1)*nn)+lgamma(nn+α)-lgamma(nn+α+2))/2)
	 		laguerre_np1 = (norm12.*exp(-sqrt(ω)*x/(2*M)).*(2*nn+1+α-sqrt(ω)*x)./(nn+1).*laguerre_n - norm13*exp(-2*sqrt(ω)*x/(2*M)).*((nn+α)/(nn+1)).*laguerre_nm1)
			laguerre_nm1 = laguerre_n
			laguerre_n   = laguerre_np1
			T[:,nn+2]=laguerre_np1.*exp((nn+1)*sqrt(ω)*x/(2*M)-sqrt(ω)*x/2)
  		end
  	end
	T = T*ω
==#
  else
    error("basis not implemented")
  end
  return T
end

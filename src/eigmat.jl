"""

T    = eigmat(basis,M,x,f,α)

Returns a matrix of mode functions of the chosen basis.
The matrix is in a form that allows the transformation
from the mode coefficients to the spatial grid `x`:

`ψ(x)=∑ⱼ cⱼϕⱼ(x) ≡ T*c`

where `c` is a column vector of coefficients in the basis.

At present `basis = "hermite"` is implemented.

## Arguments
-  `basis` is the set of eigenfunctions representing the c-field
- `M` is the number of modes in the spatial direction denoted by
- `x`, the spatial grid to which the coefficients are mapped.
- `f` is the *relative* frequence, in units of the chosen reference frequency.
- `α` is an extra input for the `laguerre` basis.
Defaults of the last two arguments are 1.

"""

function eigmat(basis,M,x,ω=1.0,α=1.0)
N = length(x)
x = x[:]

wfn_mat = zeros(N,M)

if basis=="hermite"
  hermite_nm1 = ones(size(x))
  hermite_n = complex(sqrt(2*ω)*x.*exp(-ω*x.^2/(2*M)),0)
  indx0 = find(abs(hermite_nm1).<1e-6) #avoid log of zero !
  indx1 = find(abs(hermite_nm1).>=1e-6)

  wfn_mat[indx1,1]=real(exp(log(hermite_nm1[indx1])-ω*(1-0/M)*x[indx1].^2/2))
  wfn_mat[indx0,1]=real(hermite_nm1[indx0].*exp(-ω*x[indx0].^2/2))

  if M > 1
    indx0 = find(abs(hermite_n).<1e-6)
    indx1 = find(abs(hermite_n).>=1e-6)
    wfn_mat[indx1,2]=real(exp(log(hermite_n[indx1])-ω*(1-1/M)*x[indx1].^2/2))
    wfn_mat[indx0,2]=real(hermite_n[indx0].*exp(-ω*(1-1/M)*x[indx0].^2/2))
  end

  if M > 2
    for nn = 1:M-2
        hermite_np1 = sqrt(2*ω/(nn+1))*x.*exp(-ω*x.^2/(2*M)).*hermite_n - sqrt(nn/(nn+1))*exp(-2*ω*x.^2/(2*M)).*hermite_nm1
        hermite_nm1 = hermite_n
        hermite_n   = hermite_np1   #H_n(n=nn+1)
        indx0 = find(abs(hermite_np1).<1e-6)
        indx1 = find(abs(hermite_np1).>=1e-6)
        wfn_mat[indx1,nn+2]=real(exp(log(hermite_np1[indx1])-ω*(1-(nn+1)/M)*x[indx1].^2/2))
        wfn_mat[indx0,nn+2]=real(hermite_np1[indx0].*exp(-ω*(1-(nn+1)/M)*x[indx0].^2/2))
    end
  end
wfn_mat = wfn_mat*(ω/π)^(1/4)

elseif basis=="laguerre"
    x=complex(x,0)
  	laguerre_nm1 = real(exp((α*log(x)-lgamma(α+1))/2))
	  norm12 = exp((+lgamma(1+α)-lgamma(α+2))/2)
  	laguerre_n = norm12*(1+α-ω*x).*laguerre_nm1.*exp(-ω*x/(2*M))

  	wfn_mat[:,1] = laguerre_nm1.*exp(0*ω*x/(2*M)-ω*x/2)
  	wfn_mat[:,2] = laguerre_n.*exp(1*ω*x/(2*M)-ω*x/2)

	if M > 2
  		for nn = 1:M-2
			norm12 = exp((log(nn+1)+lgamma(nn+1+α)-lgamma(nn+α+2))/2)
			norm13 = exp((log((nn+1)*nn)+lgamma(nn+α)-lgamma(nn+α+2))/2)
	 		laguerre_np1 = (norm12.*exp(-ω*x/(2*M)).*(2*nn+1+α-ω*x)./(nn+1).*laguerre_n - norm13*exp(-2*ω*x/(2*M)).*((nn+α)/(nn+1)).*laguerre_nm1)
			laguerre_nm1 = laguerre_n
			laguerre_n   = laguerre_np1
			wfn_mat[:,nn+2]=laguerre_np1.*exp((nn+1)*ω*x/(2*M)-ω*x/2)
  		end
  	end
	wfn_mat = wfn_mat*ω
else
  error("basis not implemented")
end
#wfn_mat=transpose(wfn_mat)
end

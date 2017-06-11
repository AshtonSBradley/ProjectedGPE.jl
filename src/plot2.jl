#plot solution for 2d
## Transform to cartesian grid
using Interact
R = sqrt(2μ)
xMax=1.5*R
Nx = 500
x = collect(linspace(-xMax,xMax,Nx))
Tx = eigmat("hermite",M,x,1.0,1.0)

#Plot
f=figure(figsize=(9,2))
@manipulate for i=1:length(t) withfig(f,clear=true) do
        ψ = Tx*sol[i];
        θ = unwrap(angle(ψ));
        subplot(1,2,1);plot(x/R,g*abs(ψ).^2,x/R,ones(x)*μ,"g:")
        ylim(-1,1.2*μ);xlim(-xMax/R,xMax/R)
        xlabel(L"x/R");ylabel(L"g|ψ|^2")
        subplot(1,2,2);plot(x/R,θ)
        xlim(-xMax/R,xMax/R);#ylim(-pi*1.1,pi*1.1);
        xlabel(L"x/R");ylabel(L"angle(ψ)")
    end
end

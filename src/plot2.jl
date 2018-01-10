#plot solution for 2d
## Transform to cartesian grid
@unpack ωx,ωy,ωz,Γ̄,M̄,g,t0,μ,ti,tf,Nt,t,dt = siminfo
R = sqrt(2μ/ωx^2)
xMax=1.5R
Nx = 500
x = collect(linspace(-xMax,xMax,Nx))
Tx = eigmat(M,x)

#Plot
f=figure(figsize=(9,2))
@manipulate for i=1:length(t) withfig(f,clear=true) do
        ψ = Tx*sol[i];
        θ = unwrap(angle(ψ));
        subplot(1,2,1);plot(x/R,g*abs.(ψ).^2,x/R,ones(x)*μ,"g:")
        ylim(-1,1.2*μ);xlim(-xMax/R,xMax/R)
        xlabel(L"x/R");ylabel(L"g|ψ|^2")
        subplot(1,2,2);plot(x/R,θ)
        xlim(-xMax/R,xMax/R);#ylim(-pi*1.1,pi*1.1);
        xlabel(L"x/R");ylabel(L"angle(ψ)")
    end
end

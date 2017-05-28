using ProjectedGPE
using Base.Test

# write your own tests here
function testvortexposition(Nv)
Lx = 150.
Ly = 50.
Nx = 600
Ny = 200
x = linspace(-Lx/2,Lx/2,Nx)';x=x[:]
y = linspace(-Ly/2,Ly/2,Ny)'
dx = x[2]-x[1]
dy = y[2]-y[1]

#randomly distributed vortices and charges
Nv = 10
testvort=zeros(Nv,3)

#makes sure vortices are away from edges
k=1
while k<=Nv
a = -Lx/2+Lx*rand()
b = -Ly/2+Ly*rand()
σ = rand([-1,1],1)
    if (-Lx/2+dx<a<Lx/2-dx && -Ly/2+dy<b<Ly/2-dy)
        testvort[k,:] = [a b σ]
        k+=1
    end
end

testvort = sortrows(testvort)

#construct phase
phase = zeros(Nx,Ny)
for j=1:Nv
    phase += testvort[j,3]*atan2(ones(x)*(y-testvort[j,2]),(x-testvort[j,1])*ones(y))
end


#construct wavefunction
psi = ones(x*y).*exp(im*phase)
vpos = findvortices(x,y,psi)

#check detection to grid resolution
vortfound = 0
for j=1:Nv
    if abs(testvort[j,1]-vpos[j,1])<dx && abs(testvort[j,2]-vpos[j,2])<dy
        vortfound+=1
    end
end
test = vortfound==Nv
return test
#imshow(phase)
#colorbar()
end

@test testvortexposition(30)

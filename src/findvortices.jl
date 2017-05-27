"""
 `xp,yp,xn,yn = findvortices(x,y,ψ,normal)`
 Locates vortices as 2π phase windings around plaquettes on a cartesian spatial field.
 Coorindates are returned of positive (xp,yp) and negative (xn,yn) vortices.

 ## Arguments
  * Spatial vectors x and y correspond to the column and row directsions of ψ, respectively.
  * ψ is the wavefunction on a cartesian grid.
  * `normal` refers the direction of the normal vector to the slice plane. 
  The normal vector to slice plane can point along `slice = x,y,z`.

If the field is 2D then returns vortex coordinates.
If the field is 3D then returns 2D slices of coordinates normal to direction
`slice`.
"""

function findvortices(x,y,ψ,normal="z")

phase = angle(ψ)
Nx,Ny = size(ψ)

# x corresponds to column of ψ
# y is a row vector

X = x*ones(y)
Y = ones(x)*y
vortexgrid = zeros(Nx,Ny)

for i = 1:Nx-1
    for j = 1:Nx-1

            m = 0;
            Δ = phase[i+1,j]-phase[i,j]
            abs(Δ) > π ? m += sign(Δ) : nothing

            Δ = phase[i+1,j+1]-phase[i+1,j]
            abs(Δ) > π ? m += sign(Δ) : nothing

            Δ = phase[i,j+1]-phase[i+1,j+1]
            abs(Δ) > π ? m += sign(Δ) : nothing

            Δ = phase[i,j]-phase[i,j+1]
            abs(Δ) > π ? m += sign(Δ) : nothing

            vortexgrid[i,j] = m
        end
    end
    X = x*ones(y)
    Y = ones(x)*y
    pos = vortexgrid.>0
    neg = vortexgrid.<0
    indp = find(pos)
    xp = X[indp];xp = xp[:]
    yp = Y[indp];yp = yp[:]
    indn = find(neg)
    xn = X[indn];xn = xn[:]
    yn = Y[indn];yn = yn[:]
    return xp, yp, xn, yn
end

"""
 ## findvortices(ψ,slice)
 Locates vortices as 2π phase windings around plaquettes on a cartesian spatial field.

If the field is 2D then returns vortex coordinates.
If the field is 3D then returns 2D slices of coordinates normal to direction
`slice`. The normal vector to slice plane can point along `slice = x,y,z`
"""

function findvortices(ψ,normal="z")

Phase = angle(ψ);

grid_length = size(ψ);
X = grid_length[1] - 1
Y = grid_length[2] - 1

Vortex_Grid = zeros(grid_length[1],grid_length[2]);


 for ii = 1:X-1
   for jj = 1:Y-1

            Alpha1 = Phase[ii,jj];
            β_1 = Phase[ii,jj+1];
       m = 0;

                if β_1 - Alpha1 > pi;
          m = m -1;
                    elseif β_1 - Alpha1 < -pi;
          m = m + 1;
       end

            Alpha2 = Phase[ii,jj+1];
            Beta2 = Phase[ii+1,jj+1];

         if Beta2 - Alpha2 > pi;
          m = m -1;
       elseif Beta2 - Alpha2 < -pi;
          m = m + 1;
       end


            Alpha3 = Phase[ii+1,jj+1];
            Beta3 = Phase[ii+1,jj];

         if Beta3 - Alpha3 > pi;
          m = m -1;
       elseif Beta3 - Alpha3 < -pi;
          m = m + 1;
       end


            Alpha4 = Phase[ii+1,jj];
            Beta4 = Phase[ii,jj];

         if Beta4 - Alpha4 > pi;
          m = m -1;
       elseif Beta4 - Alpha4 < -pi;
          m = m + 1;
         end

            Vortex_Grid[ii,jj] = m;

   end
 end
        return Vortex_Grid
end

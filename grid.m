function grid

global nc ncm k kout ka iord itest stab % era USE NUM_PAR
global b c diverg c1 c2 x0 dx           % era USE GEOM_PAR
global x                                % era USE GEOM_VAR
      
      b=0.0;       % left boundary position
      c=1.0;       % right boundary position
      dx=(c-b)/nc; % calcolo del passo per mesh equispaziata

      for n=1:nc
        x(n) = b + 0.5d0*dx + dx*(n-1);  % cell center position
      end
          
end
function U=vortex(u,x,y,sign)

% input: u: phi field
%        x
%        y
%        sign: (+1) for counterclockwise vortex
% ouput: U: phi field with a vortex at (x,y) whose L1 radius is infinity

    U=u;
    [m,n]=size(u);
    for i=1:m
        for j=1:n
            alpha=atan2(i-x,j-y)+pi/2*sign;
            U(i,j)=U(i,j)+alpha;
        end
    end

end
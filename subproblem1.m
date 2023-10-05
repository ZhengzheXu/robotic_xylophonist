function solution = subproblem1(w,r,p,q)
    u = p - r;
    v = q - r;
    u_proj = u - w*w'*u;
    v_proj = v - w*w'*v;
    if u_proj == 0
        disp('Infinite solutions');
    end
    a = w'*(cross(u_proj,v_proj));
    b = u_proj'*v_proj;
    solution = atan2d(a,b);
end
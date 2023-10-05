function solution = subproblem3(w, r, p, q, d)
    u = p - r;
    v = q - r;
    
    up = u - w*w'*u;
    vp = v - w*w'*v;
    dp_sq = d^2 - abs(dot(w,p-q))^2;

    value1 = dot(w,cross(up,vp));
    value2 = dot(up,vp);
    value3 = (dot(up,up) + dot(vp,vp) - dp_sq)/(2*norm(up)*norm(vp));
    theta0 = atan2d(value1,value2);
    if value3 == 1
        solution = theta0;
        return;
    end
    theta1 = theta0 + acosd(value3);
    theta2 = theta0 - acosd(value3);
    solution = [theta1 theta2];
end
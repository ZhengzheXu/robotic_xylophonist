% AUTHOR:       Xu Zhengzhe
% 
% ABSRACT:      
% 
% INPUT: xi: 6x1 twist vector
%        theta: rotation angle
% OUTPUT: g: 4x4 transformation matrix
% 

function g = Transformation(xi, theta)
    theta = theta * pi / 180;
    v = xi(1:3);
    w = xi(4:6);
    w_hat = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
    R = eye(3) + sin(theta)*w_hat + (1-cos(theta))*w_hat^2;
    p = (eye(3) - R)*cross(w,v) + w * w' * v * theta;
    g = [R p; 0 0 0 1];
end
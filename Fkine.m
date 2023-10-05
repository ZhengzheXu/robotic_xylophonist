% AUTHOR :  Xu Zhengzhe
%
% ABSTRACT： 这是计算机器人正解函数，通用函数，需调用Transformation函数
% 
% INPUT： Xi      机器人各关节运动旋量, 6XN矩阵，单位m和rad
%         theta   机器人关节位移，1xN向量，单位m和rad
%         g0      机器人基准参考位姿， 4X4矩阵
% OUTPUT: g_st    机器人末端位姿位姿， 4X4矩阵
% 
function g_st = Fkine(Xi,theta,g0)
    n = size(Xi,2);
    g_st = g0;
    for i = n:-1:1
        g_st = Transformation(Xi(:,i),theta(i))*g_st;
    end
end
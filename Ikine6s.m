% AUTHOR : Xu Zhengzhe
%
% ABSTRACT：这是六轴机器人逆运动学函数，适用于六轴机器人
% 
% INPUT：g0               机器人初始末端姿态矩阵，4x4矩阵
%        gd               机器人目标末端姿态矩阵，4x4矩阵
%        config           机器人杆长信息
%        
% OUTPUT: theta           关节转角矩阵，1X6，单位rad      
%

function theta = Ikine6s(g0,gd,config)
    L1 = config{1}; L2 = config{2}; L3 = config{3}; L4 = config{4};
    q1 = config{5}; q2 = config{6}; q3 = config{7}; q4 = config{8};
    q5 = config{9}; q6 = config{10};
    w1 = config{11}; w2 = config{12}; w3 = config{13}; w4 = config{14};
    w5 = config{15}; w6 = config{16};
    theta1_range = config{17}; theta2_range = config{18}; theta3_range = config{19};
    theta4_range = config{20}; theta5_range = config{21}; theta6_range = config{22};

    xi1 = [-cross(w1,q1); w1];
    xi2 = [-cross(w2,q2); w2];
    xi3 = [-cross(w3,q3); w3];
    xi4 = [-cross(w4,q4); w4];
    xi5 = [-cross(w5,q5); w5];
    xi6 = [-cross(w6,q6); w6];

    g1 = gd/g0;
    pw = [0 0 L1+L2+L3]'; % wrist center
    pb = [0 0 L1]'; % base center

    theta = [];

    % disp("---------Start to solve IK problem---------");
    % disp("inital position: ");
    % disp(g0);
    % disp("target position: ");
    % disp(gd);

    % STEP1: solve for theta3
    % disp("[STEP1] Solve for theta3");
    delta = g1*[pw; 1] - [pb; 1];
    d = norm(delta);
    solution1 = subproblem3(w3, q3, pw, pb, d);
    solution1 = [solution1, solution1+360, solution1-360];

    % STEP2: solve for theta1, theta2
    q = g1*[pw; 1];
    q = q(1:3);
    for i = 1:length(solution1)
        this_theta3 = solution1(i);
        if this_theta3 < theta3_range(1) || this_theta3 > theta3_range(2)
            % disp("theta3 out of range");
            continue
        end
        % disp("theta3: " + this_theta3);

        this_g3 = Transformation(xi3, this_theta3);
        this_g3_inv = inv(this_g3);

        p = Transformation(xi3, this_theta3)*[pw; 1];
        p = p(1:3);
        % disp("[STEP2] Solve for theta1, theta2");
        solution2 = subproblem2(w1,w2,q2,p,q);
        if isempty(solution2)
            continue
        end
        solution2 = [solution2, solution2+360, solution2-360];

        for j = 1:size(solution2,2)
            this_theta1 = solution2(1,j);
            if this_theta1 < theta1_range(1) || this_theta1 > theta1_range(2)
                % disp("theta1 out of range");
                continue
            end

            this_g1 = Transformation(xi1, this_theta1);
            this_g1_inv = inv(this_g1);

            this_theta2 = solution2(2,j);
            if this_theta2 < theta2_range(1) || this_theta2 > theta2_range(2)
                % disp("theta2 out of range");
                continue
            end

            this_g2 = Transformation(xi2, this_theta2);
            this_g2_inv = inv(this_g2);

            % STEP3: solve for theta4
            % g1g2g3g4g5g6 gst0 = gd
            % g1g2g3g4g5g6 = gd * gst0^-1 = (g1)
            % 已知g1g2g3，求g4g5g6
            % g4g5g6 = (g1g2g3)^-1 * gd * gst0^-1, 右边都是已知量
            % 对于只有g4g5g6，它的y轴坐标永远是0
            % 因此求(g1g2g3)^-1 * gd * gst0^-1的第2个坐标
            % disp("[STEP3] Solve for theta4");
            g4g5g6 = this_g3_inv * this_g2_inv * this_g1_inv * g1;
            qq = g4g5g6 * [0.01,0.01,0.01,1]';
            x = qq(1);
            y = qq(2);
            if x==0
                theta40 =-90;
            else
                theta40 = atand(y/x);
            end

            for m=-2:1:2
                theta4 = theta40 + m*180;
                if theta4 >= theta4_range(1) && theta4 <= theta4_range(2)
                    this_g4 = Transformation(xi4, theta4);
                    this_g4_inv = inv(this_g4);
                    % STEP4: solve for theta5
                    % disp("[STEP4] Solve for theta5");
                    % using subproblem 1
                    rhs = this_g4_inv * this_g3_inv * this_g2_inv * this_g1_inv * g1 * [q6; 1];
                    qqq = rhs(1:3);
                    solution5 = subproblem1(w5,q5,q6,qqq);
                    solution5 = [solution5, solution5+360, solution5-360];
                    % STEP5: solve for theta6
                    % disp("[STEP5] Solve for theta6");
                    for n=1:length(solution5)
                        this_theta5 = solution5(n);
                        if this_theta5 < theta5_range(1) || this_theta5 > theta5_range(2)
                            % disp("theta5 out of range");
                            continue
                        end
                        this_g5 = Transformation(xi5, this_theta5);
                        this_g5_inv = inv(this_g5);
                        p22 = [100 0 0]';
                        rhs2 = this_g5_inv * this_g4_inv * this_g3_inv * this_g2_inv * this_g1_inv * g1 * [p22;1];
                        q22 = rhs2(1:3);
                        solution6 = subproblem1(w6,q6,p22,q22);
                        solution6 = [solution6, solution6+360, solution6-360];
                        for nn=1:length(solution6)
                            this_theta_6 = solution6(nn);
                            if this_theta_6 < theta6_range(1) || this_theta_6 > theta6_range(2)
                                % disp("theta6 out of range");
                                continue
                            end
                            theta = [theta; this_theta1 this_theta2 this_theta3 theta4 this_theta5 this_theta_6];
                            % disp("solution found")
                            % disp(theta)
                        end
                    end
                end
            end
        end
    end

    % disp("---------End of IK problem---------");
    % if isempty(theta)
    %     disp("No solution found");
    % else
    %     disp("Solution (deg): ");
    %     disp(theta);
    %     disp("Solution (rad): ");
    %     disp(theta*pi/180);
%     end

end





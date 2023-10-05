clc;clear;close all

% Define the mapping from notes to stack numbers
note_map = containers.Map({'G3', 'A3' 'B3' 'C4' 'D4' 'E4' 'F4' 'G4' 'A4' 'B4' 'C5' 'D5' 'E5' 'F5' 'G5'}, ...
                          {15, 14, 13, 12, 11, 10, 9, ...
                           8, 7, 6, 5, 4, 3, 2, 1});

% Canon
bpm = 35;
beat = 4; % 4/4
duration = [
            2,2,...
            2,2,...
            2,2,...
            2,2,...
            2,2,...
            2,2,...
            2,2,...
            2,2,...
            2,4,4,...
            2,4,4,...
            2,4,4,...
            2,4,4,...
            8,8,8,8,4,4,...
            4,4,8,8,8,8,...
            8,8,8,8,8,8,8,8,...
            8,8,4,8,8,8,8,...
            8,8,8,8,4,4,...
            4,4,8,8,8,8,...
            8,8,8,8,8,8,8,8,...
            8,8,4,16/3,16,8,16,16,...
            16/3,16,16,16,16,4,16,16,16,16,16,...
            16/3,16/3,16,16,8,16,16/3,8,...
            4,8,8,4,8,8,...
            4,16,16,8,16/3,8,16,16,16,...
            16/3,16,16,16,16,4,16,16,16,16,16,...
            16/3,16/3,16,16,8,16,8,8,...
            4,8,8,4,8,8,...
            4,16,16,8,16/3,8,16,16,16,...
            8,16,16,8,16,16,16,16,16,16,16,16,16,16,...
            8,16,16,8,16,16,16,16,16,16,16,16,16,16,...
            8,16,16,8,16,16,16,16,16,16,16,16,16,16,...
            8,16,16,8,16,16,16,16,16,16,16,16,16,16,...
            8,16,16,8,16,16,16,16,16,16,16,16,16,16,...
            8,16,16,8,16,16,16,16,16,16,16,16,16,16,...
            8,16,16,8,16,16,16,16,16,16,16,16,16,16,...
            8,16,16,8,16,16,16,16,16,16,16,16,16,16,...
            8,16,16,8,16,16,16,16,16,16,16,16,16,16,...
            8,16,16,8,16,16,16,16,16,16,16,16,16,16,...
            8,16,16,8,16,16,16,16,16,16,16,16,16,16,...
            8,16,16,8,16,16,16,16,16,16,16,16,16,16,...
            1
            ];
melody = {
        'E4','D4',...
        'C4','B3',...
        'A3','G3',...
        'A3','B3',...
        'E5','D5',...
        'C5','B4',...
        'A4','G4',...
        'A4','B4',...
        'E5','D5','F4',...
        'C5','B4','G4',...
        'A4','G4','E4',...
        'A4','G4','F4',...
        'C5','B4','C5','E4','G4','B4',...
        'C5','E5','G4','E4','G4','A4',...
        'F5','E5','D5','F5','E5','D5','C5','B4',...
        'A4','F4','C5','B4','G4','C5','B4',...
        'C5','B4','C5','E4','G4','B4',...
        'C5','E5','G4','E4','G4','A4',...
        'F5','E5','D5','F5','E5','D5','C5','B4',...
        'A4','G4','F4','C5','G4','B4','D5','G5',...
        'E5','G4','E5','D5','C5','D5','E5','F5','E5','D5','E5',...
        'C5','C5','B4','C5','B4','G4','E4','G4',...
        'A4','B4','C5','G4','E4','G4',...
        'A4','F4','A4','C5','C5','B4','C5','D5','G4',...
        'E5','G4','E5','D5','C5','D5','E5','F5','E5','D5','E5',...
        'C5','C5','B4','C5','B4','G4','E4','G4',...
        'A4','B4','C5','G4','E4','G4',...
        'A4','F4','A4','C5','C5','B4','C5','D5','G4',...
        'G5','E5','F5','G5','E5','F5','G5','B4','A4','B4','C5','D5','E5','F5',...
        'E5','C5','D5','E5','E4','F4','G4','A4','G4','F4','G4','C5','B4','C5',...
        'A4','C5','B4','A4','G4','F4','G4','F4','E4','F4','G4','A4','B4','C5',...
        'A4','C5','B4','C5','B4','C5','B4','A4','B4','C5','D5','E5','F5','G5',...
        'G5','E5','F5','G5','E5','F5','G5','B4','A4','B4','C5','D5','E5','F5',...
        'E5','C5','D5','E5','E4','F4','G4','A4','G4','F4','G4','C5','B4','C5',...
        'A4','C5','B4','A4','G4','F4','G4','F4','E4','F4','G4','A4','B4','C5',...
        'A4','C5','B4','C5','B4','C5','B4','A4','B4','C5','D5','E5','F5','G5',...
        'E5','C5','D5','E5','D5','C5','D5','B4','C5','D5','E5','D5','C5','B4',...
        'C5','A4','B4','C5','C5','D5','E5','F5','E5','D5','E5','C5','B4','C5',...
        'A4','C5','B4','A4','G4','F4','G4','F4','E4','F4','G4','A4','B4','C5',...
        'A4','C5','B4','C5','B4','A4','B4','C5','D5','C5','B4','C5','A4','B4',...
        'C5'
        };

% Stack position
left_upper = [458 162.15];
right_lower = [458 -192];
stack_size = [60 20];
row = 15;
y_step = (left_upper(2)-right_lower(2))/(row-1);
knock_height = 554;
raise_height = knock_height + 20;
knock_pitch = 154;
raise_pitch = knock_pitch + 12;
eul = [175 156 50];

% Calculate the position of each block
block_position = zeros(row,3);
for i = 1:row
    block_position(i,:) = [left_upper(1), left_upper(2)-(i-1)*y_step, knock_height];
end

% Store the position of each note
waypoints = zeros(length(melody),3);
pitch_waypoints = zeros(length(melody),1);
for i = 1:length(melody)
    if ~strcmp(melody{i},' ')
        disp(melody{i});
        waypoints(i,:) = block_position(note_map(melody{i}),:);
        
    end
    pitch_waypoints(i) = knock_pitch;
end

% Deal with the empty notes
for i = 1:length(melody)
    if strcmp(melody{i},' ')
        if i == 1
            % Find the first non-empty note
            for j = i+1:length(melody)
                if ~strcmp(melody{j},' ')
                    waypoints(i,:) = block_position(note_map(melody{j}),:);
                    waypoints(i,3) = raise_height;
                    break;
                end
            end
        else
            % Middle position
            if i == length(melody)
                waypoints(i,:) = block_position(note_map(melody{i-1}),:);
                waypoints(i,3) = raise_height;
            else
                waypoints(i,:) = (waypoints(i-1,:)+waypoints(i+1,:))/2;
                waypoints(i,3) = raise_height;
            end
        end
        pitch_waypoints(i) = raise_pitch;
    end
end

% Insert raise_height waypoints between each two waypoints
% (x,y) is the average of the two waypoints, z is raise_height
insert_waypoints = [];
for i = 1:length(melody)
    insert_waypoints = [insert_waypoints; waypoints(i,1:2)];
end
insert_waypoints = [insert_waypoints,ones(size(insert_waypoints,1),1)*raise_height];
pitch_insert_waypoints = ones(size(insert_waypoints,1),1)*raise_pitch;

% Merge new_waypoints and Waypoints
new_waypoints = [];
for i = 1:length(melody)
    new_waypoints = [new_waypoints;insert_waypoints(i,:);waypoints(i,:)];
end
% new_waypoints = [new_waypoints;waypoints(end,:)];
waypoints = new_waypoints;

new_pitch_waypoints = [];
for i = 1:length(melody)
    new_pitch_waypoints = [new_pitch_waypoints;pitch_insert_waypoints(i);pitch_waypoints(i)];
end

pitch_waypoints = new_pitch_waypoints;

% Recover position
recover_waypoints = [left_upper(1),0,raise_height+20];
waypoints = [recover_waypoints;waypoints;recover_waypoints];
pitch_waypoints = [raise_pitch;pitch_waypoints;raise_pitch];

% Deal with duration
dur_sec = 60/bpm*beat./duration;
min_duration = min(dur_sec);
knock_time = min_duration/1.5;
dur_wo_knock = dur_sec - knock_time;

% Merge dur_wo_knock and knock_time
new_dur_sec = [];
for i = 1:length(dur_wo_knock)
    new_dur_sec = [new_dur_sec;knock_time;dur_wo_knock(i)];
end
new_dur_sec = new_dur_sec(1:end-1);
dur_sec = new_dur_sec;

% Add recover time
recover_time = 3;
recover_time2 = 5;
dur_sec = [recover_time;dur_sec];
dur_sec = [dur_sec;recover_time2];

% Interpolation
time_step = 0.005;
dt = time_step;
interp_waypoints = [];
interp_pitch_waypoints = [];

% Calculate the time of each node
time_node = [0];
for i = 1:length(dur_sec)
    time_node = [time_node,time_node(end)+round(dur_sec(i)/time_step)*time_step];
end

% Use pchip interpolation to interpolate the waypoints of each axis
for i = 1:3
    this_interp = pchip(time_node,waypoints(:,i),0:time_step:time_node(end)-time_step);
    interp_waypoints = [interp_waypoints,this_interp'];
end

this_pitch_interp = pchip(time_node,pitch_waypoints,0:time_step:time_node(end)-time_step);
interp_pitch_waypoints = [interp_pitch_waypoints,this_pitch_interp'];

points_num = size(interp_waypoints,1);
time = linspace(0,points_num*time_step,points_num);
num_kp = points_num;

figure(2)
plot3(interp_waypoints(:,1),interp_waypoints(:,2),interp_waypoints(:,3),'r*')
hold on
plot3(interp_waypoints(:,1),interp_waypoints(:,2),interp_waypoints(:,3),'r')

% % % Draw the stacks
% % count = 1;
% % for i = 1:row
% %         pos = [left_upper(1), left_upper(2)-(i-1)*y_step-stack_size(2)/2, stack_size(1), stack_size(2)];
% %         % 画出矩形
% %         rectangle('Position',pos,'Curvature',[0.2 0.2])
% %         % 画出矩形的中心点
% %         hold on
% %         plot(pos(1)+stack_size(1)/2,pos(2)+stack_size(2)/2,'r*')
% %         hold on
% %         % 画出矩形的编号
% %         text(pos(1)+stack_size(1)/2,pos(2)+stack_size(2)/2+stack_size(2)/6,num2str(count))
% %         count = count + 1;
% % end
% % axis equal

% % [p1,p2,p3,p4,p5,p6] = calKeyPoints(begin_2d_pos, ending_2d_pos,knock_height,raise_height);
% % points = [p1',p2',p3',p4',p5',p6'];

% Define the manipulator parameters
L1 = 491; % mm
L2 = 450; % mm
L3 = 450; % mm
L4 = 84; % mm

q1 = [0 0 0]';
q2 = [0 0 L1]';
q3 = [0 0 L1+L2]';
q4 = [0 0 L1+L2+L3]';
q5 = [0 0 L1+L2+L3]';
q6 = [0 0 L1+L2+L3+L4]';

w1 = [0 0 1]';
w2 = [0 1 0]';
w3 = [0 1 0]';
w4 = [0 0 1]';
w5 = [0 1 0]';
w6 = [0 0 1]';

theta1_range = [-170 170]; % deg
theta2_range = [-120 120]; % deg
theta3_range = [-228 148]; % deg
theta4_range = [-170 170]; % deg
theta5_range = [-120 120]; % deg
theta6_range = [-360 360]; % deg

config = {L1 L2 L3 L4...
    q1 q2 q3 q4 q5 q6...
    w1 w2 w3 w4 w5 w6...
    theta1_range theta2_range theta3_range...
    theta4_range theta5_range theta6_range};

% Define the initial position
p_st0 = [0 0 L1+L2+L3+L4]'; % mm
rot_st0 = [0 0 -180]; % deg
R_st0 = eul2rotm(rot_st0*pi/180,"ZYZ");
g0 = [R_st0, p_st0; 0 0 0 1];

solutions = cell(1,num_kp);
init_joint = [0.324 7.962 109 -2.29 81.5 -38.8];

for i=1:num_kp
    disp("processing percent:")
    disp(i/num_kp*100)
    % Define the target position
    p_d = interp_waypoints(i,:)';
    this_eul = [eul(1),interp_pitch_waypoints(i),eul(3)];
    R_d = eul2rotm(this_eul*pi/180,"ZYZ");
    gd = [R_d, p_d; 0 0 0 1];

    % Calculate the inverse kinematics
    solution = Ikine6s(g0,gd,config);
    if isempty(solution)
        disp('No solution for the point:')
        disp(i)
        disp(interp_waypoints(i,:))
        return
    end

    % Find the closest solution 
    if i ==1
        this_joint = init_joint;
    else
        this_joint = solutions{i-1};
    end
    % disp("size(solution,1)")
    % disp(size(solution,1))
    for j=1:size(solution,1)
        this_solution = solution(j,:);
        this_norm = norm(this_solution-this_joint);
        if j == 1
            min_norm = this_norm;
            min_solution = this_solution;
        else
            if this_norm < min_norm
                min_norm = this_norm;
                min_solution = this_solution;
            end
        end
    end
    solutions{i} = min_solution;
end

disp('=====================================')
disp('The path planning is done.')

sol_mat = cell2mat(solutions);
sol_mat = reshape(sol_mat,6,num_kp)';
theta1_key = sol_mat(:,1);
theta2_key = sol_mat(:,2);
theta3_key = sol_mat(:,3);
theta4_key = sol_mat(:,4);
theta5_key = sol_mat(:,5);
theta6_key = sol_mat(:,6);
theta_key = [theta1_key,theta2_key,theta3_key,theta4_key,theta5_key,theta6_key];

% Trajectory Generation
t = 0:dt:(num_kp-1)*dt; 
q = theta_key; % Joint angle vector
q_smooth = q;
for i=1:6
    q_smooth(:,i) = smooth(q(:,i),25);
end
q = q_smooth;

qd = zeros(length(t),6); % Joint angular velocity
for i=1:6
    qd(:,i) = [diff(q(:,i))/dt; qd(end,i)];
end

qdd = zeros(length(t),6); % Joint angular acceleration
for i=1:6
    qdd(:,i) = [diff(qd(:,i))/dt; qdd(end,i)];
end

q = deg2rad(q);
qd = deg2rad(qd);
qdd = deg2rad(qdd);
t = t';

figure(4)
col = 3; row = 2;
for i=1:col
    for j=1:row
        subplot(row,col,(j-1)*col+i)
        plot(t,q_smooth(:,i+(j-1)*col),'b','LineWidth',2)
        title(['theta',num2str(i+(j-1)*col)])
    end
end

figure(5)
col = 3; row = 2;
for i=1:col
    for j=1:row
        subplot(row,col,(j-1)*col+i)
        plot(t,qd(:,i+(j-1)*col),'b','LineWidth',2)
        title(['Angular velocity of theta',num2str(i+(j-1)*col)])
    end
end

figure(6)
col = 3; row = 2;
for i=1:col
    for j=1:row
        subplot(row,col,(j-1)*col+i)
        plot(t,qdd(:,i+(j-1)*col),'b','LineWidth',2)
        title(['Angular acceleration of theta',num2str(i+(j-1)*col)])
    end
end

disp("total time: ")
disp(t(end))

% Write to point file
q_deg = rad2deg(q);
points_data = q_deg;
dt_new = 0.001;
points_data = interp1(1:size(points_data,1),points_data,linspace(1,size(points_data,1),t(end)/dt_new),'spline');
time_data = 0:dt_new:t(end)-dt_new;
% figure(7)
% plot(time_data,points_data)

% Save to file
fileID = fopen('xzz_slow.txt','w');
fprintf(fileID,'%f %f %f %f %f %f\r\n',points_data');
fclose(fileID);
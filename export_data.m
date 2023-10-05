% 创建点位文件

points_data = q*180/pi;
% 对points_data进行插值，每两个点之间插入999个点
points_data = interp1(1:size(points_data,1),points_data,1:0.001:size(points_data,1));
figure(1)
% 画出6个关节的角度变化曲线
subplot(3,2,1)
plot(points_data(:,1))
title('joint1')
subplot(3,2,2)
plot(points_data(:,2))
title('joint2')
subplot(3,2,3)
plot(points_data(:,3))
title('joint3')
subplot(3,2,4)
plot(points_data(:,4))
title('joint4')
subplot(3,2,5)
plot(points_data(:,5))
title('joint5')
subplot(3,2,6)
plot(points_data(:,6))
title('joint6')

% 保存到文件
fileID = fopen('points1.txt','w');
fprintf(fileID,'%f %f %f %f %f %f\r\n',points_data');
fclose(fileID);

fileID = fopen('points2.txt','w');
fprintf(fileID,'%f %f %f %f %f %f\r',points_data');
fclose(fileID);

fileID = fopen('points3.txt','w');
fprintf(fileID,'%f %f %f %f %f %f\n',points_data');
fclose(fileID);

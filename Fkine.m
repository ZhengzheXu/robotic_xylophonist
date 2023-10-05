% AUTHOR :  Xu Zhengzhe
%
% ABSTRACT�� ���Ǽ�����������⺯����ͨ�ú����������Transformation����
% 
% INPUT�� Xi      �����˸��ؽ��˶�����, 6XN���󣬵�λm��rad
%         theta   �����˹ؽ�λ�ƣ�1xN��������λm��rad
%         g0      �����˻�׼�ο�λ�ˣ� 4X4����
% OUTPUT: g_st    ������ĩ��λ��λ�ˣ� 4X4����
% 
function g_st = Fkine(Xi,theta,g0)
    n = size(Xi,2);
    g_st = g0;
    for i = n:-1:1
        g_st = Transformation(Xi(:,i),theta(i))*g_st;
    end
end
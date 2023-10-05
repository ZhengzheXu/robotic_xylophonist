function [p1,p2,p3,p4,p5,p6] = calKeyPoints(start_2d_pos, end_2d_pos,pick_height,raise_height)
    p1 = [start_2d_pos(1), start_2d_pos(2), raise_height];
    p2 = [start_2d_pos(1), start_2d_pos(2), pick_height];
    p3 = p1;
    p4 = [end_2d_pos(1), end_2d_pos(2), raise_height];
    p5 = [end_2d_pos(1), end_2d_pos(2), pick_height];
    p6 = p4;
end

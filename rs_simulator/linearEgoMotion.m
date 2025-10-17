% INPUT ARGUMENTS
%   -R0: rotation matrix for 1st row - R_world_to_camera 
%   -t0: translation vector for 1st tow - t_world_to_camera  
%   -w: roatation velocity - w_world_to_camera 
%   -d: translation velocity - d_world_to_camera  
%   -cameraParams: MATLAB define cameraParameters
% OUTPUT ARGUMENTS
%   -Rot_Rows: n cells of (3x3) rotation matrix - R_world_to_camera 
%   -trans_Rows: n cells of (3x1) transaltion vecotr - t_world_to_camera  

function [Rot_Rows, trans_Rows] = linearEgoMotion(R0,t0,w,d, cameraParams)

%% Data Prepare 
% image size
width = cameraParams.IntrinsicMatrix(3,1) * 2;
height = cameraParams.IntrinsicMatrix(3,2) * 2;

%% Ego-motion matrices generation 
% Rot_Rows = {};
% trans_Rows = {};
% for row_index = 1 : height
%     Rot_Rows{row_index} = (eye(3) + X_(w) * row_index) * R0; 
%     trans_Rows{row_index} = t0 + d * row_index;
% end




Rot_Rows = {};
trans_Rows = {};
last_r = eye(3);
last_t = t0;
if norm(w) == 0
    one_r = eye(3);
else
    one_r = axang2rotm([(w / norm(w))' norm(w)]);
end


for row_index = 1 : height
    last_t = (d + last_t);
    last_r = one_r * last_r;
    cur_t = last_t;
    cur_r = last_r * R0;
    Rot_Rows{row_index} = cur_r; 
    trans_Rows{row_index} = cur_t;
end

% Rot_Rows = {};
% trans_Rows = {};
% last_r = eye(3);
% last_t = zeros(3, 1);
% if norm(w) == 0
%     one_r = eye(3);
% else
%     one_r = axang2rotm([(w / norm(w))' norm(w)]);
% end
% 
% 
% for row_index = 1 : height
%     last_t = last_t + last_r * d;
%     last_r = one_r * last_r;
%     cur_t = last_t + t0;
%     cur_r = last_r * R0;
%     Rot_Rows{row_index} = cur_r; 
%     trans_Rows{row_index} = cur_t;
% end
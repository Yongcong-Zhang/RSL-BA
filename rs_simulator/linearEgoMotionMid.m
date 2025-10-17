% INPUT ARGUMENTS
%   -R0: rotation matrix for 1st row - R_world_to_camera 
%   -t0: translation vector for 1st tow - t_world_to_camera  
%   -w: roatation velocity - w_world_to_camera 
%   -d: translation velocity - d_world_to_camera  
%   -cameraParams: MATLAB define cameraParameters
% OUTPUT ARGUMENTS
%   -Rot_Rows: n cells of (3x3) rotation matrix - R_world_to_camera 
%   -trans_Rows: n cells of (3x1) transaltion vecotr - t_world_to_camera  

function [Rot_Rows, trans_Rows] = linearEgoMotionMid(R0,t0,w,d, cameraParams)

%% Data Prepare 
% image size
width = cameraParams.IntrinsicMatrix(3,1) * 2;
height = cameraParams.IntrinsicMatrix(3,2) * 2;

%% Ego-motion matrices generation 
Rot_Rows = {};
trans_Rows = {};
for row_index = 1: height
    Rot_Rows{row_index} = (eye(3) + X_(w) * (row_index - (height / 2))) * R0; 
    trans_Rows{row_index} = t0 + d * (row_index - (height / 2));
end

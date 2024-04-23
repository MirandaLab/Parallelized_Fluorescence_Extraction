function [Spherical_vol_cyt] = OAM_220820_Get_Sphere_Vol_cyt(mask_cyt_Red)
%% uses a binary mask of a cell or cellular structures as a matrix to calculate the volumen of a sphere with the same equivalent diameter
mask_cyt_RedA=(bwlabel(mask_cyt_Red));
ed=regionprops(mask_cyt_RedA,'EquivDiameter'); % Computed as sqrt(4*Area/pi)
Spherical_vol_cyt=0.523*(ed.EquivDiameter)^3;

%figure;imagesc(mask_cyt_RedA)
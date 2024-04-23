function [Spherical_vol_nuc] = OAM_220820_Get_Sphere_Vol_nuc(mask_nuc_Red)
%% uses a binary mask of a cell or cellular structures as a matrix to calculate the volumen of a sphere with the same equivalent diameter

mask_nuc_RedA=(bwlabel(mask_nuc_Red));
ed=regionprops(mask_nuc_RedA,'EquivDiameter'); % Computed as sqrt(4*Area/pi)
Spherical_vol_nuc=0.523*(ed.EquivDiameter)^3;


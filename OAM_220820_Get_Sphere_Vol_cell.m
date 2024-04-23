function [Spherical_vol_cell] = OAM_220820_Get_Sphere_Vol_cell(ccell)
%% uses a binary mask of a cell or cellular structures as a matrix to calculate the volumen of a sphere with the same equivalent diameter
ccellA=(bwlabel(ccell));
ed=regionprops(ccellA,'EquivDiameter'); % Computed as sqrt(4*Area/pi)
Spherical_vol_cell=0.523*(ed.EquivDiameter)^3;


clear all
parpool('local',4)

type='.tif';
path_h='C:\Users\vulpine\Desktop\A1\';
foldrs=dir(path_h);
% remove all files (isdir property is 0)
foldrs = foldrs([foldrs(:).isdir]);
% remove '.' and '..' 
foldrs = foldrs(~ismember({foldrs(:).name},{'.','..','Tracks','MAT_Tracks','TET_Tracks'}));

fold_name = { foldrs.name };

cha='.png';

%============load segmentation=============================================
for it0=1:size(fold_name,2)
    path_im=[path_h fold_name{it0}];
    path_seg=[path_h 'Tracks\'];

    foldrs1=dir(fullfile(path_im,['*' cha]));
    fold_name1 = { foldrs1.name };
    final_time=size(fold_name1,2);

    load([path_seg  fold_name{it0} '_ART_Track_' num2str(final_time)])
    
    x_size =size(Mask2{1,1},1);   %1040;
    y_size =size(Mask2{1,1},2);   %1388;
    
    cell_Vol                      =zeros(no_obj,final_time);
    %CELL morphology=======================================================
    cell_area                    =zeros(no_obj,final_time); % area in pixels
    cell_perimeter=zeros(no_obj,final_time);
    cell_solidity=zeros(no_obj, final_time);
    cell_convex_area=zeros(no_obj, final_time);
    cell_equiv_diameter=zeros(no_obj, final_time);
    cell_extent=zeros(no_obj, final_time);
    cell_eulernumber=zeros(no_obj, final_time);
    cell_convex_hull=zeros(no_obj,final_time);
    cell_eccentricity=zeros(no_obj,final_time);
    %texture
    cell_Contrast=zeros(no_obj,final_time);
    cell_Homogeneity=zeros(no_obj,final_time);
    cell_Correlation=zeros(no_obj,final_time);
    cell_Energy=zeros(no_obj,final_time);
    cell_Variance=zeros(no_obj,final_time);
    %There are more Haralick measures if you want to add later.

    for c_time=final_time:-1:1
        disp(c_time)
        Lcells=Mask2{1,c_time}; %  imagesc(Lcells) a broadcast variable
        
        image_number=sprintf('%09d',c_time-1);

        parfor cell_no=1:no_obj
            if c_time>=cell_exists(1,cell_no) &&  sum(sum(Lcells==cell_no))~=0 % >= must be to avoid firsst column being zeros
                disp(cell_no)
                ccell=(Lcells==cell_no);
                ccell=bwmorph(ccell,"clean");% get rid of small peripheral pixels imagesc(ccell)
%                 ccell=imopen((Lcells==cell_no),ones(3)); %%  imagesc(ccell-ccell1)
                cell_Vol(cell_no, c_time) = OAM_220820_Get_Sphere_Vol_cell(ccell);
                cell_margin=1;
                [x_cn,y_cn]=get_wind_coord1(ccell,cell_margin);
                
           
                    p_Ccell=double(ccell);
%                     imagesc(p_Ccell)
                    ccell2=double(p_Ccell(y_cn,x_cn));
                   % figure(4)
                   % imagesc(ccell2)
                    

%                         imagesc(ccell2)
                    if isempty(ccell2)==0 && sum(ccell2(:))~=0
                        props1=regionprops(ccell2,'EulerNumber');
                        cell_eulernumber(cell_no,c_time)=props1.EulerNumber;
%                         figure(5)
%                         imagesc(ccell2)
                        %=============  data for vacuoles  ============================================
                        cell_area(cell_no,c_time)       = sum(ccell2(:));
                        props=regionprops(ccell2,'Perimeter','Solidity','Area','ConvexArea','Solidity','EquivDiameter','Extent','Eccentricity','ConvexHull');
                        cell_perimeter(cell_no,c_time)=props.Perimeter; % Distance around the boundary of the region returned as a scalar
                        cell_solidity(cell_no, c_time)=props.Solidity; %  Area/ConvexArea.
                        cell_convex_area(cell_no, c_time)=props.ConvexArea; % 	Number of pixels in 'ConvexImage'
                        cell_equiv_diameter(cell_no, c_time)=props.EquivDiameter; % Diameter of a circle with the same area as the region,sqrt(4*Area/pi).
                        cell_extent(cell_no, c_time)=props.Extent; % Ratio of pixels in the region to pixels in the total bounding box
                        cell_convex_hull(cell_no, c_time)=props.ConvexArea; % perimeter sqrt / sphere surface area, all in pixels
                        cell_eccentricity(cell_no,c_time)=props.Eccentricity; % 0= circle, 1= extremely elogated elipse (a line)
                                             
                        glcms = graycomatrix(ccell2,'NumLevels',length(unique(ccell2)),'GrayLimits',[],'Offset',[0 1; -1 1;-1 0;-1 -1]); %gray-level coocuurrence matrices
                        stats = graycoprops(glcms,{'Contrast','Homogeneity','Correlation','Energy'});
%                         figure(6)
%                         imagesc(ccell2)
                        
                        %Texture Features
                        cell_Contrast(cell_no,c_time)=mean(stats.Contrast); %  intensity contrast between a pixel and its neighbor()
                        cell_Homogeneity(cell_no,c_time)=mean(stats.Homogeneity); % correlated a pixel is to its neighbor
                        cell_Correlation(cell_no,c_time)=mean(stats.Correlation); % 
                        cell_Energy(cell_no,c_time)=mean(stats.Energy);% sum of squared elements in the closeness of the distribution of elements in the GLCM to the GLCM diagonal.
                        
                       else
                        disp('No cell')
                        disp(['cell_no' num2str(cell_no)] )
                        disp(['time' num2str(c_time)] )
                    end
                    %==============================================================================
  
            end %end if cell exists
        end %parfor

    end %time-loop
    
    %add all the fields to all_obj structure
    all_obj.appr_vol =cell_Vol;
    all_obj.cell_area                       =cell_area;
    %additional
    all_obj.cell_perimeter=cell_perimeter;
    all_obj.cell_solidity=cell_solidity;
    all_obj.cell_convex_area=cell_convex_area;
    all_obj.cell_equiv_diameter=cell_equiv_diameter;
    all_obj.cell_extent=cell_extent;
    all_obj.cell_eulernumber=cell_eulernumber;
    all_obj.cell_roundness_convex_hull=cell_convex_hull;
    all_obj.cell_eccentricity=cell_eccentricity;
    all_obj.cell_Contrast=cell_Contrast;
    all_obj.cell_Homogeneity=cell_Homogeneity;
    all_obj.cell_Correlation=cell_Correlation;
    all_obj.cell_Energy=cell_Energy;
       
    name1=[fold_name{it0} '_morpho_Seg_' num2str(final_time) ];
    save(fullfile(path_seg,name1), 'all_obj');
end


delete(gcp('nocreate'))




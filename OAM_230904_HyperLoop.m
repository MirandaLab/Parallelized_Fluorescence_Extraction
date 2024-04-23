




clear 
I_mean_modifier=1.5;
peak_cutoff=0.75;
path_h='D:\Master_London\London_Exps\OAM211015_SK1_HTB1_RIM11';
foldrs=dir(path_h);
% remove all files (isdir property is 0)
foldrs = foldrs([foldrs(:).isdir]);
% remove '.' and '..' 
foldrs = foldrs(~ismember({foldrs(:).name},{'.','..','Tracks','New','X'}));
fold_name = { foldrs.name };


fold3=['D:\Master_London\London_Exps\OAM211015_SK1_HTB1_RIM11\' fold_name{3} ];
path_seg='D:\Master_London\Extraction\OAM211015_SK1_HTB1_RIM11\';

load([path_seg 'OAM211015_SK1_HTB1_RIM11_' fold_name{3} '_110'])

    x_size = size(all_obj.cells(:,:,1),1);
    y_size = size(all_obj.cells(:,:,1),2);



%if its the first one
ALLDATA{1,1}='Channels';
ALLDATA{2,1}='cell_vol';

file_n=dir(fullfile(fold3, '*.tif'));
file_n2={file_n.name};

Name=cell(1,1);
    for it0=1:7
        Name{it0,1}=char(file_n2{it0}(1,14:end)); % 14 is the position in the file_name after img_000000000
    end

channels=unique(Name);

for it1=1:size(channels,1) % channels

file1=dir(fullfile(fold3, ['*' channels{it1}]));
file2={file1.name};
% allocation 
cell_Vol1=zeros(no_obj,size(all_obj.twoD_area,2)); 
max_nuc_int1=zeros(no_obj,size(all_obj.twoD_area,2)); 
mean_cell_Fl1=zeros(no_obj,size(all_obj.twoD_area,2)); 
Conc_T_cell_Fl1=zeros(no_obj,size(all_obj.twoD_area,2)); 
mem_area1=zeros(no_obj,size(all_obj.twoD_area,2)); 
nuc_area1=zeros(no_obj,size(all_obj.twoD_area,2)); 
cyt_area1=zeros(no_obj,size(all_obj.twoD_area,2)); 
mean_fl_mem1=zeros(no_obj,size(all_obj.twoD_area,2)); 
std_fl_mem1=zeros(no_obj,size(all_obj.twoD_area,2)); 
tot_Fl_mem1=zeros(no_obj,size(all_obj.twoD_area,2)); 
tot_Fl_cyt1=zeros(no_obj,size(all_obj.twoD_area,2)); 
tot_Fl_nuc1=zeros(no_obj,size(all_obj.twoD_area,2)); 
mean_int_per_area_C1=zeros(no_obj,size(all_obj.twoD_area,2)); 
mean_int_per_area_N1=zeros(no_obj,size(all_obj.twoD_area,2)); 
nuc_Vol1=zeros(no_obj,size(all_obj.twoD_area,2)); 
cyt_Vol1=zeros(no_obj,size(all_obj.twoD_area,2)); 
cyt_Vol_sub1=zeros(no_obj,size(all_obj.twoD_area,2)); 
FL_Conc_T1=zeros(no_obj,size(all_obj.twoD_area,2)); 
FL_Conc_C1=zeros(no_obj,size(all_obj.twoD_area,2)); 
FL_Conc_N1=zeros(no_obj,size(all_obj.twoD_area,2)); 
FL_mean_int_N_thr1=zeros(no_obj,size(all_obj.twoD_area,2)); 
FL_mean_int_C_thr1=zeros(no_obj,size(all_obj.twoD_area,2)); 


for c_time=1:10%;size(all_obj.twoD_area,2) % images

cell_Vol=zeros(no_obj,1); 
max_nuc_int=zeros(no_obj,1); 

    Lcells=all_obj.cells(:,:,c_time); % figure;imagesc(Lcells) % or Mask2{1,it2}
     I = imread([fold3 '\' file2{c_time}]);
     I = double(I); % figure;imagesc(I)
     I = medfilt2(I,'symmetric'); % filtering 
     bck = (I.*(~Lcells)); %figure;imagesc(bck)
     backgr =   median(bck(bck~=0));% 
     I = (I-backgr); % background correction

%      parvar=find(cell_exists(:,2)<=c_time);

  parfor cell_no = 1:10%no_obj % cell_no = 4
        
%        if ctime==1 || ctime>cell_exists(cell_no,2)
                ccell = (Lcells == cell_no); %  figure;imagesc(ccell); figure;imagesc(Lcells);
      if sum(ccell(:))~=0
                cell_margin = 5;
                [x_cn, y_cn] = get_wind_coord1(ccell, cell_margin);
                ccell=ccell(y_cn, x_cn);%  figure;imagesc(ccell);
                cell_Vol(cell_no,1) = OAM_220820_Get_Sphere_Vol_cell(ccell);

                    I_cell = I(y_cn,x_cn); % figure;imagesc(I_cell)
                    put_I = ccell.*I_cell; % figure;imagesc(put_I)
                    max_nuc_int(cell_no,1) = max(put_I(:)); % figure;imagesc(put)

                    mean_cell_Fl(cell_no,1) = sum((put_I(:)))./sum(ccell(:)); % total intensity/number of cell pixels
                    Conc_T_cell_Fl(cell_no,1) = sum((put_I(:)))./cell_Vol(cell_no,1); % total intensity / approx. volume
                    % Obtain nucleus by Gaussian fit  imagesc(mask_nuc)
                    % this gaussian corrects for multiple equally intense brightest pixels 
                    %[mask_nuc] = OAM_231006_Gaussian_nuclear_fit(I_cell,peak_cutoff); 
                    [mask_nuc] = OAM_221006_Gaussian_nuclear_fit(I_cell,peak_cutoff,x_size,y_size,Lcells,ccell);  
                 
                    % visualize:  figure;imagesc(mask_nuc)
                    mask_mem = bwmorph(ccell, 'remove'); % figure;imagesc(mask_mem)
                    mem_area(cell_no,1) = sum(mask_mem(:));  % area of the membrane

if mask_nuc~=0
%                     mask_cyt = double(bwmorph(ccell.*bwmorph(~mask_nuc,'thicken',1),'erode',1));  
                    mask_cyt = double(ccell-mask_nuc);%.*bwmorph(~mask_nuc,'thicken',1),'erode',1));  
else
 mask_cyt =nan;
end
                    % visualize:  figure;imagesc(p_nuc) figure;imagesc(mask_cyt)
                    nuc_area(cell_no,1) = sum(mask_nuc(:));
                    cyt_area(cell_no,1) = sum(mask_cyt(:));
                    mem_fl = mask_mem.*I_cell; % figure;imagesc(mem_fl);
                    mean_fl_mem(cell_no,1) = median(mem_fl(mem_fl~=0));
                    std_fl_mem(cell_no,1) = std(mem_fl(mem_fl~=0)); % figure;imagesc(mem+ccell)
                    tot_Fl_mem(cell_no,1) = sum(mem_fl(:));
                    tot_Fl_cyt(cell_no,1) = sum(sum(mask_cyt.*I_cell));
                    tot_Fl_nuc(cell_no,1) = sum(sum(mask_nuc.*I_cell));
                    mean_int_per_area_C(cell_no,1) = sum(sum(mask_cyt.*I_cell))./sum(mask_cyt(:));
                    mean_int_per_area_N(cell_no,1) = sum(sum(mask_nuc.*I_cell))./nuc_area(cell_no,1);
                    nuc_Vol(cell_no,1) = OAM_220820_Get_Sphere_Vol_nuc(mask_nuc);
                    cyt_Vol(cell_no,1) = OAM_220820_Get_Sphere_Vol_cyt(mask_cyt);% figure;imagesc(mask_cyt)
                    cyt_Vol_sub(cell_no,1)  = cell_Vol(cell_no,1) - nuc_Vol(cell_no,1);
                    FL_Conc_T(cell_no,1)     = sum(put_I(:))./ cell_Vol(cell_no,1);
                    FL_Conc_C(cell_no,1)     = tot_Fl_cyt(cell_no,1)./cyt_Vol(cell_no,1);
                    FL_Conc_N(cell_no,1)     = tot_Fl_nuc(cell_no,1)./nuc_Vol(cell_no,1);
                    % Obtain nucleus by threshold constructed with "I_mean_modifier"
                    ccell = put_I;  % visualize: figure;imagesc(put)
                    put_mod = (ccell>(I_mean_modifier.*mean(ccell(ccell>0))));
                    put_mod = bwareaopen(put_mod,5,4);
                    put_mod = imfill(put_mod,'holes');  % figure;imagesc(put_mod)
                    FL_mean_int_N_thr(cell_no,1) = sum(sum(put_mod.*ccell))./sum(put_mod(:));
                    no = ccell.*(~put_mod); % visualize: imagesc(no)
                    FL_mean_int_C_thr(cell_no,1) = sum(no(no>0))./sum(no(no>0)>0);

       else 
       end
% 
   end %parfor

cell_Vol1(:,c_time)=cell_Vol;
max_nuc_int1(:,c_time)=max_nuc_int;
mean_cell_Fl1(:,c_time)=mean_cell_Fl; 
Conc_T_cell_Fl1(:,c_time)=Conc_T_cell_Fl; 
mem_area1(:,c_time)=mem_area; 
nuc_area1(:,c_time)=nuc_area; 
cyt_area1(:,c_time)=cyt_area; 
mean_fl_mem1(:,c_time)=mean_fl_mem; 
std_fl_mem1(:,c_time)=std_fl_mem; 
tot_Fl_mem1(:,c_time)=tot_Fl_mem; 
tot_Fl_cyt1(:,c_time)=tot_Fl_cyt; 
tot_Fl_nuc1(:,c_time)=tot_Fl_nuc; 
mean_int_per_area_C1(:,c_time)=mean_int_per_area_C; 
mean_int_per_area_N1(:,c_time)=mean_int_per_area_N; 
nuc_Vol1(:,c_time)=nuc_Vol; 
cyt_Vol1(:,c_time)=cyt_Vol; 
cyt_Vol_sub1(:,c_time)=cyt_Vol_sub; 
FL_Conc_T1(:,c_time)=FL_Conc_T; 
FL_Conc_C1(:,c_time)=FL_Conc_C; 
FL_Conc_N1(:,c_time)=FL_Conc_N; 
FL_mean_int_N_thr1(:,c_time)=FL_mean_int_N_thr; 
FL_mean_int_C_thr1(:,c_time)=FL_mean_int_C_thr; 


end     %cell - for

%if its the first one
ALLDATA{1,it1+1}=channels{it1};
ALLDATA{2,it1+1}=cell_Vol1;

ALLDATA{3,it1+1}=max_nuc_int1; 
ALLDATA{4,it1+1}=mean_cell_Fl1;
ALLDATA{5,it1+1}=Conc_T_cell_Fl1;
ALLDATA{6,it1+1}=mem_area1;
ALLDATA{7,it1+1}=nuc_area1;
ALLDATA{8,it1+1}=cyt_area1;
ALLDATA{9,it1+1}=mean_fl_mem1;
ALLDATA{10,it1+1}=std_fl_mem1;
ALLDATA{11,it1+1}=tot_Fl_mem1;
ALLDATA{12,it1+1}=tot_Fl_cyt1;
ALLDATA{13,it1+1}=tot_Fl_nuc1;
ALLDATA{14,it1+1}=mean_int_per_area_C1;
ALLDATA{15,it1+1}=mean_int_per_area_N1;
ALLDATA{16,it1+1}=nuc_Vol1;
ALLDATA{17,it1+1}=cyt_Vol1;
ALLDATA{18,it1+1}=cyt_Vol_sub1;
ALLDATA{19,it1+1}=FL_Conc_T1;
ALLDATA{20,it1+1}=FL_Conc_C1;
ALLDATA{21,it1+1}=FL_Conc_N1;
ALLDATA{22,it1+1}=FL_mean_int_N_thr1;
ALLDATA{23,it1+1}=FL_mean_int_C_thr1;
ALLDATA{24,it1+1}=all_back;


end

delete(gcp('nocreate'));


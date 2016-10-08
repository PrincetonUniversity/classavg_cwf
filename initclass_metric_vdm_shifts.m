% Example code for generating class averages.
%
% The example generates 10000 noisy projection and produces 10000 class
% averages by averaging each image with its 50 nearest neighbors. 
%
% Data set:
% ---------
% 10000 noisy projections with additive white Gaussian noise at SNR/1=100.
% The 
% projections orientations are random and uniformly disributed, and they
% are shifted randomly with maximum shift of +/- 4 pixels in the x and y
% direcrtions. The clean projections are separated in 20 different defocus
% groups.
% If the file 'clean_data.mat' (in the folder ./simulation) is missing,
% please go into "ASPIRE/projections/class_average/" and run  
%   gen_simulation_data 
% to generate it, before running this example code. 
%
%
% Tuning the class averaging algorithm:
% -------------------------------------
% Variables determined by the user are the following
%
%   r_max =floor(L/2)-10; 
%       Radius of region of interest that contains the particle.
%
%   n_nbor = 50; 
%       Number of nearest neighbors for initial classification.
%
% Use new FB coefficients, new Init_class function
% Tejal April 17, 2016
run ~/aspire/initpath.m
run ~/cwf_denoise/cwf_paths.m

tic;
clear all;
K = 10000; %K is the number of images
SNR = 1/40; %SNR
use_shifted=0;

if(use_shifted)
load('/scratch/ARCHIVE_from_sdl6/tbhamre/cwf_class/clean_data_6454_65_shift3.mat'); % load clean centered projection images 
else
load('/scratch/ARCHIVE_from_sdl6/tbhamre/cwf_class/clean_data_6454_65.mat'); % load clean centered projection images 
end
disp('Loaded clean data')
%downsampleddim=65;
%sprintf('Downsampling to %dX%d grid', downsampleddim, downsampleddim)
%data.projections=cryo_downsample(data.projections,[downsampleddim downsampleddim],1);
use_CTF=1;
ndef=20; % Number of defocus groups
def1=1;
def2=4;
lambda = EWavelength(300);
B=10; % decay envelope parameter

[g_proj_CTF,CTF,defocus_group]=  add_CTF_env_v6(cfft2(data.projections(:,:,1:K)), ndef, def1,def2,B, lambda, use_CTF);
[images, noise_v_r]=addnoise_v6(icfft2(g_proj_CTF), SNR);
% Use this for test with clean data
q = data.q(:, 1:K);
%images = data.projections(:, :, 1:K);
%clear data;
L = size(images, 1);
n_nbor = 10; %number of nearest neighbors for initial classification.
n_nbor_large=50;
isrann = 0;

use_VDM=1;
k_VDM_in = n_nbor; % number of nearest neighbors for building graph for VDM.
VDM_flag = 0; % VDM using union rule
k_VDM_out = n_nbor; % number of nearest neighbors search for 
% Initial Classification with CWF
%[CWF_data, cwf_coeff_cell, denoised_coeff_ccwf, basis]=data_CWF(images, CTF, defocus_group, noise_v_r, ndef, def1, def2, B, lambda, use_CTF);
%[recon_cwf] = recon_images_FB(CWF_data.c, CWF_data.R, L, denoised_coeff_ccwf, 1, size(images,3)); % Specify range of images to reconstruct
%[mse_cwf] = calc_MSE_v6(recon_cwf, data.projections(:,:,1:K),CWF_data.R);
%
%% Optiional: Reduce coefficients corresponding to those (k,q) for which the L2 energy is low (average over all images)
%l2_thresh=0.99;
%[CWF_data_red]=reduce_coeffs(cwf_coeff_cell, CWF_data.Freqs, l2_thresh); 
%[ class, class_refl, rot, corr,  timing ] = Initial_classification_FD(CWF_data_red, n_nbor, isrann );
%%[ class, class_refl, rot, corr,  timing ] = Initial_classification_FD(CWF_data, n_nbor, isrann );
%disp('Finished initial classification...');
%% Check Classification result
%[ d, error_rot ] = check_simulation_results(class, class_refl, -rot, q); % should use minus sign for init class, no minus sign for VDM 
%[ N, X ] = hist(acosd(d), [0:180]);
%figure; bar(N); title('CWF')
%xlabel('a$\cos\langle v_i, v_j \rangle$', 'interpreter', 'latex');
%sprintf('CWF: Number of images with correlation > %f is %d',0.9, numel(find(d(d>0.9))))
%

% Initial classification with sPCA (new, fast code)
[ images_fl ] = Phase_Flip(images, defocus_group, CTF); %phase flipping 
disp('Phase flipped');
[sPCA_data, sPCA_coeff_cell, basis, recon_spca]=data_sPCA(images_fl,  noise_v_r);
[mse_spca] = calc_MSE_v6(recon_spca, data.projections(:,:,1:K),sPCA_data.R);
tic_init=tic;
[ class_f, class_refl_f, rot_f, corr_f,  timing_f ] = Initial_classification_FD(sPCA_data, n_nbor, isrann );
toc_init=toc(tic_init);
disp('Finished initial classification...');
if(use_VDM)
	tic_VDM = tic;
	[ class_VDM, class_VDM_refl, rot_f_vdm ] = VDM(class_f, ones(size(class_f)), rot_f, class_refl_f, k_VDM_in, VDM_flag, k_VDM_out);
	toc_VDM = toc(tic_VDM);
	disp('Finished VDM classification...');
	% Check Classification result
	[ d_f, error_rot_f ] = check_simulation_results(class_VDM, class_VDM_refl, rot_f_vdm, q); % should use minus sign for init class, no minus sign for VDM 
else
	[ d_f, error_rot_f ] = check_simulation_results(class_f, class_refl_f, -rot_f, q); % should use minus sign for init class, no minus sign for VDM 
end
%[ N_f, X_f ] = hist(acosd(d_f), [0:180]);
%figure; bar(N_f); title('sPCA')
%xlabel('a$\cos\langle v_i, v_j \rangle$', 'interpreter', 'latex');
%sprintf('CWF: Number of images with correlation > %f is %d',0.9, numel(find(d(d>=0.9))))
sprintf('sPCA: Number of images with correlation > %f is %d',0.9, numel(find(d_f(d_f>=0.9))))


%% With new Mahalonobis distance
% Initial classification with sPCA (new, fast code)
[ class_f_large, class_refl_f_large, rot_f_large, corr_f_large,  timing_f_large ] = Initial_classification_FD(sPCA_data, n_nbor_large, isrann );
disp('Finished initial classification...');
if(use_VDM)
	tic_VDM = tic;
	[ class_VDM_large, class_VDM_refl_large, rot_f_vdm_large ] = VDM(class_f_large, ones(size(class_f_large)), rot_f_large, class_refl_f_large, k_VDM_in, VDM_flag, k_VDM_out);
	toc_VDM = toc(tic_VDM);
	disp('Finished VDM classification...');
	% Check Classification result
	[ d_f_large, error_rot_f_large ] = check_simulation_results(class_VDM_large, class_VDM_refl_large, rot_f_vdm_large, q); % should use minus sign for init class, no minus sign for VDM 
else
	[ d_f_large, error_rot_f_large ] = check_simulation_results(class_f_large, class_refl_f_large, -rot_f_large, q); % should use minus sign for init class, no minus sign for VDM 
end
% Check Classification result
%[ N_f_large, X_f_large ] = hist(acosd(d_f_large), [0:180]);
%figure; bar(N_f_large); title('sPCA large')
%xlabel('a$\cos\langle v_i, v_j \rangle$', 'interpreter', 'latex');
%sprintf('CWF: Number of images with correlation > %f is %d',0.9, numel(find(d(d>=0.9))))
%sprintf('sPCA: Number of images with correlation > %f is %d',0.9, numel(find(d_f(d_f>=0.9))))
list_recon = [1:size(images, 3)];

if(use_shifted)
max_shift=3;
else
max_shift=0;
end

k_out=n_nbor_large;
tic_align = tic;
new_num_nn=n_nbor;
[data_cwf] =  data_cwf_metric(images, CTF, defocus_group, noise_v_r, ndef, def1, def2, B, lambda, use_CTF);
[ shifts, corr, average, norm_variance, class_m, class_refl_m, rot_m ] = align_main_cwf( images_fl,...
 rot_f_large, class_f_large, class_refl_f_large, sPCA_data, k_out, max_shift, list_recon, recon_spca, data_cwf, defocus_group, new_num_nn); % Should it be rot_f or -rot_f?

toc_align = toc(tic_align);
[ d_m, error_rot_m ] = check_simulation_results(class_m, class_refl_m, -rot_m, q); % should use minus sign for init class, no minus sign for VDM 
sprintf('sPCA + new metric: Number of images with correlation > %f is %d',0.9, numel(find(d_m(d_m>=0.9))))
%% Initial Classification with old FBsPCA
sprintf('sPCA: Number of images with correlation > %f is %d',0.9, numel(find(d_f(d_f>=0.9))))
SNR
%%[ images ] = Phase_Flip(images, defocus_group, CTF); %phase flipping 
%%disp('Phase flipped');
%r_max =sPCA_data.R; %radius of region of interest. Determined by user.
%% Check Classification result
%[ class_o, class_refl_o, rot_o, ~, FBsPCA_data, timing ] = Initial_classification(images, r_max, n_nbor, isrann );
%[ d_o, error_rot_o ] = check_simulation_results(class_o, class_refl_o, -rot_o, q); % should use minus sign for init class, no minus sign for VDM 
%[ N_o, X_o ] = hist(acosd(d_o), [0:180]);
%figure; bar(N_o); title('FBPCA')
%xlabel('a$\cos\langle v_i, v_j \rangle$', 'interpreter', 'latex');
%sprintf('FBsPCA old: Number of images with correlation > %f is %d',0.9, numel(find(d_o(d_o>0.9))))
%
%d_val=[0.8:0.05:1];
%corr_comp=[];
%for i=1:length(d_val)
%	corr_cwf(i)=numel(find(d(d>=d_val(i))));
%	corr_spca(i)=numel(find(d_f(d_f>=d_val(i))));
%	corr_oldspca(i)=numel(find(d_o(d_o>=d_val(i))));
%	corr_comp=vertcat(corr_comp,[corr_cwf(i), corr_spca(i), corr_oldspca(i)]);
%end
%h=bar(d_val,corr_comp);
%colormap(summer(5));
%grid on
%l = cell(1,3);
%l{1}='CWF'; l{2}='New sPCA'; l{3}='Old FBsPCA'    
%legend(h,l);
%xlabel('Correlation')
%ylabel('Number of pairs with correlation greater than')
%title(sprintf('%d images, %d nearest neighbors, SNR=1/%d',K , n_nbor, (1/SNR)))
%

% Compare MSE
%sprintf('MSE with CWF is %f',mse_cwf)
sprintf('MSE with sPCA is %f',mse_spca)

(toc./(60*60*24))

old_wins=zeros(10000,1);
new_wins=zeros(10000,1);

d_m1=reshape(d_m,10000,n_nbor);
d_f1=reshape(d_f,10000,n_nbor);

votes=zeros(10000,1);


corrthr=0.995

for i=1:10000
old_wins(i)=numel(d_f1(d_f1(i,:)>corrthr));
new_wins(i)=numel(d_m1(d_m1(i,:)>corrthr));
end

for i=1:10000
if(old_wins(i)<new_wins(i))
	votes(i)=1;
elseif(old_wins(i)>new_wins(i))
	votes(i)=2;
end
end

fprintf('Perc of wins for Mahalanobis %.2f',numel(votes(votes==1))/10000)
fprintf('Perc of wins for old class avg %.2f',numel(votes(votes==2))/10000)
fprintf('Perc of draws %.2f',numel(votes(votes==0))/10000)

 

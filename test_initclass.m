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

% Tejal April 17, 2016

clc;
clear all;
clf;
K = 10000; %K is the number of images
%SNR = 1/10; %SNR
data = load('clean_data.mat'); % load clean centered projection images 
disp('Loaded clean data')
%[images, defocus_group, noise, noise_spec, c, q]=create_image_wCTF(data, SNR, 'gaussian'); %create projection images with CTF and shifts
q = data.q(:, 1:K);
images = data.projections(:, :, 1:K);
clear data;
%[ images ] = Phase_Flip(images, defocus_group, c); %phase flipping 
%disp('Phase flipped');
% Low pass filtering images to make it approximately invaraint to small

% shifts
%[ images_lpf ] = low_pass_filter( images );
L = size(images, 1);
r_max =floor(L/2)-10; %radius of region of interest. Determined by user.
n_nbor = 10; %number of nearest neighbors for initial classification.
isrann = 0;
% Initial Classification
[ class, class_refl, rot, ~, FBsPCA_data, timing ] = Initial_classification(images, r_max, n_nbor, isrann );
disp('Finished initial classification...');
% Check Classification result
[ d, error_rot ] = check_simulation_results(class, class_refl, -rot, q); % should use minus sign for init class, no minus sign for VDM 
[ N, X ] = hist(acosd(d), [0:180]);
figure; bar(N);
xlabel('a$\cos\langle v_i, v_j \rangle$', 'interpreter', 'latex');

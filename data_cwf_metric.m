function [data_cwf] =  data_cwf_metric(images, CTF, index, noise_v_r, ndef, def1, def2, B, lambda, use_CTF)
% Prepare denoised coefficients and frequencies data
% images: Images in real space, with noise and CTF 

proj_CTF_noisy=cfft2(images);
%% Calculate coeffs using FFBsPCA
energy_thresh=0.99;
[ c, R ] = choose_support_v6( proj_CTF_noisy, energy_thresh); %Estimate band limit and compact support size
c=c*(0.5/floor(size(proj_CTF_noisy,1)/2)); % Rescaling between 0 and 0.5
c=0.5    
sprintf('hard coded c=%f for testing',c)
n_r = ceil(4*c*R);
tic_basis=tic;
[ basis, sample_points ] = precomp_fb( n_r, R, c );
timing.basis=toc(tic_basis)
num_pool=5;
L0=size(images,1);

regu=1;
mean_image_f=mean_LS(CTF, index, proj_CTF_noisy,regu);

[y_mu] = demean_y_v6(proj_CTF_noisy, CTF, mean_image_f, index);
[ coeff_ymu ] = coeff_demean( icfft2(y_mu) , R, basis, sample_points, num_pool);
[ coeff_mean ] = coeff_demean( icfft2(mean_image_f) , R, basis, sample_points, num_pool);

%% CTF in new basis: numerical integration

[ctf_rad_all]=  calc_CTF_rad(use_CTF, L0, index, ndef, def1,def2,B, lambda, sample_points.r*((floor(L0/2)+1)/0.5));
[ data_cwf ]  = get_cwf_metric(index, ctf_rad_all, basis, sample_points,  coeff_mean, coeff_ymu,  noise_v_r);
%energy_coeff=cumsum(cellfun(@norm,denoised_coeff_ccwf));
%cutoff_coeff=find(energy_coeff/max(energy_coeff)>0.99,1,'first');
%compressed_den_coeff=denoised_coeff_ccwf(1:cutoff_coeff);

Coeff_all=cell2mat(data_cwf.Coeff);
%data.Lmat=L;
data_cwf.Freqs=basis.ang_freqs(1:size(Coeff_all,1));
data_cwf.c=c;
data_cwf.R=R;

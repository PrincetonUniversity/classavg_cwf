function [ Coeff_b, toc_bispec ] = Bispec_2Drot_large( Coeff, Freqs )
%Description:
%This function computes rotationally invariant bispectral-like features for 2D images.
%   Input: 
%       Coeff:
%           Truncated expansion coefficients of 2D images on Fourier Bessel
%           steerable PCA basis.
%       Freqs:
%           The angular frequencies associated with each component.
% 
%   Output:
%       Coeff_b:
%           The invariant feature for original images.
%       toc_bispec:
%           Timing for bispectrum computation.
%
%   Zhizhen Zhao Aug 2013
%   Updated Jan 2015 for complex images

tic_bispec=tic;
alpha=1/3; %modify the amplitude for each component.
n_im = size(Coeff, 2);
Coeff_norm=(abs(Coeff(Freqs~=0, :))).^(alpha);
Coeff_norm=log(Coeff_norm);
check=isinf(Coeff_norm);
assert(max(check(:))~=1);
clear check
Phase=Coeff(Freqs~=0, :)./abs(Coeff(Freqs~=0, :));
Phase=atan2(imag(Phase), real(Phase));
clear Coeff
[O1, O2] = bispec_Operator(Freqs(Freqs~=0));
fprintf('\nLength of the bispectrum is %d \n', size(O1, 1));

%random sample (construct bispectrum with a subset of the the images);
shuffle_id = randperm(size(Coeff_norm, 2));
subnum = 5000; % SHould be a function of the size of the dataset
M=exp(O1*Coeff_norm(:, shuffle_id(1:min(n_im, subnum))) + sqrt(-1)*O2*Phase(:, shuffle_id(1:min(n_im, subnum)))); % Memory bottleneck

% Think about this part, not final
%If the length of bispectrum is too big, subsample, just choose components with the largest variance 
variance_bispec = var(M, [], 2);
[sorted_var, id ] = sort(variance_bispec, 'descend');
l_O1 = ceil(2*length(Freqs)*log10(length(Freqs))); 
O1 = O1(id(1:l_O1), :);
O2 = O2(id(1:l_O1), :);
P_max = floor(7*10^8/l_O1); %just define the size of the images that might be able to hold in memory
% P_max shoudl be user defined/ check user memory
% %% PCA the bispectrum

ncomp = 400;
% Should be chosen at runtime?

%dimensionality reduction using PCA, need to think about subsample in images and in bispectrum coefficients
if n_im<= P_max
    M = exp(O1*Coeff_norm+sqrt(-1)*O2*Phase));
    if l_O1<=ncomp
        Coeff_b = M;
    else
        [ ~, S, V ]=pca_Y(M, ncomp);
        bar(diag(S));
        fprintf('\nFinished PCA');
        Coeff_b = S*V';
    end;
else
    Coeff_b = zeros(min(l_O1, ncomp), n_im);   
    if l_O1<=ncomp
        for i=1:ceil(n_im/P_max)
            Coeff_b(:, (i-1)*P_max+1:min(n_im, i*P_max)) = exp(O1*Coeff_norm(:, (i-1)*P_max+1:min(n_im, i*P_max)) + sqrt(-1)*O2*Phase(:, (i-1)*P_max+1:min(n_im, i*P_max)));
        end;
    else
        M = exp(O1*(Coeff_norm(:, shuffle_id(1:min(n_im, P_max)))+sqrt(-1)*Phase(:, shuffle_id(1:min(n_im, P_max)))));
        [ U, ~, ~ ] = pca_Y(M, ncomp);
        for i=1:ceil(n_im/P_max)
            M = exp(O1*Coeff_norm(:, (i-1)*P_max+1:min(n_im, i*P_max)) + sqrt(-1)*O2*Phase(:, (i-1)*P_max+1:min(n_im, i*P_max)));
            Coeff_b(:, (i-1)*P_max+1:min(n_im, i*P_max)) = U'*M;
        end;
    end;
end;

clear O1 O2 M

%Normalize 
for i=1:size(Coeff_b, 2);
     Coeff_b(:, i)=Coeff_b(:, i)/norm(Coeff_b(:, i));
end;
toc_bispec=toc(tic_bispec);
end


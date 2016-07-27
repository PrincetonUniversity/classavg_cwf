function [ shifts, corr, average, norm_variance, class_met, class_refl_met, angle_met ] = align_main_cwf( data, angle, class_VDM, refl, FBsPCA_data, k, max_shifts, list_recon, recon_sPCA, data_cwf, def_grp, num_nn)
% Function for aligning images with its k nearest neighbors to generate
% class averages.
%   Input: 
%       data: LxLxP matrix. P projection images of size LxL pixels.
%       angle: Pxl (l>=k) matrix. Rotational alignment angle
%       class_VDM: Pxl matrix. Nearest neighbor list
%       refl: Pxl matrix. 1: no reflection. 2: reflection
%       FBsPCA_data: Fourier-Bessel steerable PCA data with r_max, UU,
%       Coeff, Mean, Freqs
%       k: number of nearest neighbors for class averages
%       max_shifts: maximum number of pixels to check for shift
%       list_recon: indices for images to compute class averages
%   Output:
%       shifts: Pxk matrix. Relative shifts for k nearest neighbors
%       corr: Pxk matrix. Normalized cross correlation of each image with
%       its k nearest neighbors
%       average: LxLxP matrix. Class averages
%       norm_variance: compute the variance of each class averages.
%
% Zhizhen Zhao Feb 2014

P=size(data, 3);
L=size(data, 1);
l=size(class_VDM, 2);

N=floor(L/2);
[x, y]=meshgrid(-N:N, -N:N);
r=sqrt(x.^2+y.^2);

r_max = FBsPCA_data.R;
UU = FBsPCA_data.U;
Coeff = FBsPCA_data.Coeff;
Mean = FBsPCA_data.Mean;
Freqs = FBsPCA_data.Freqs;
clear FBsPCA_data;
Lmat=data_cwf.Lmat;
cwf_coeff=data_cwf.Coeff;
cwf_coeff_pr=data_cwf.Coeff_pr;
clear data_cwf;

%Check if the number of nearest neighbors is too large
if l<k
    error('myApp:argChk', 'The number of nearest neighbors is too large. \nIt should be smaller or equal to %d', l)
end;

shifts=zeros(length(list_recon), k+1);
corr=zeros(length(list_recon), k+1);
average=zeros(L, L, length(list_recon));
norm_variance=zeros(length(list_recon), 1);

%generate grid. Precompute phase for shifts
range=-fix(L/2):fix(L/2);
[omega_x,omega_y]=ndgrid(range,range);
omega_x=-2*pi.*omega_x/L; omega_y=-2*pi.*omega_y/L;
omega_x=omega_x(:); omega_y=omega_y(:);
a=-max_shifts:1:max_shifts; % checking shift for every other pixel;
num=length(a);
a1=repmat(a', num, 1);
a2=kron(a', ones(num, 1));
shifts_list=[a1, a2];

phase= exp(sqrt(-1)*(omega_x*shifts_list(:, 1)'+omega_y*shifts_list(:, 2)'));

angle=round(-angle);
angle(angle<0) = angle(angle<0)+360;

angle(angle==360)=0;
angle_met=angle;
%Precompute rotation table, precision for rotation is 1 degree
M=cell(359);
for i=1:359
    M{i}=fastrotateprecomp(L, L,i);
end;


%Go concurrent
ps=parpool('local', 12);
if ps==0
  parpool('local',12);
end


parfor j=1:length(list_recon)
    
    angle_j=angle(list_recon(j), 1:k); %rotation alignment
 
    refl_j=refl(list_recon(j), 1:k);    %reflections
    
    index = class_VDM(list_recon(j), 1:k);
    
    images=data(:, :, index); % nearest neighbor images
    
    image1=data(:, :, list_recon(j));
    
    %Build denoised images from FBsPCA
    %reconstruct the images.
    %tmp = 2*real(UU(:, Freqs~=0)*Coeff(Freqs~=0, list_recon(j)));
    %tmp = tmp + UU(:, Freqs==0)*real(Coeff(Freqs==0, list_recon(j)));
    tmp=recon_sPCA(:,:,list_recon(j));
    I=tmp;
    %I = zeros(L);
    %I(r<=r_max)=tmp;  % Not needed as well, all this is taken care of while reconstructing images from sPCA
    %I = I+Mean;  % No need to add mean image separately, recon_spca already contains final reconstructed image
    %I = mask_fuzzy(I, r_max-5); % What's the need for this? Ask Jane
    
    for i=1:k
        if (refl_j(i)==2)
            images(:, :, i)=flipud(images(:, :, i));
        end;
    end;
    for i=1:k
        if (angle_j(i)~=0)
            images(:, :, i)=fastrotate(images(:, :, i), angle_j(i), M{angle_j(i)});
        end
    end;
    
    pf1 = cfft2(I);
    pf1 = pf1(:);
    
    images = cat(3, image1, images);
    pf_images=zeros(size(images));
    for i=1:k+1
        pf_images(:, :, i)=cfft2(images(:, :, i));
    end;
    pf_images=reshape(pf_images, L^2, k+1);
    
    
    pf=bsxfun(@times, phase, pf1);

    C=pf'*pf_images;
    [corr(j, :), id ] = max(C, [], 1);
    
    pf_images_shift=pf_images.*conj(phase(:, id));
    variance=var(pf_images_shift, [], 2);
    norm_variance(j)=norm(variance, 'fro');
    tmp = mean(pf_images_shift, 2);
    tmp = reshape(tmp, L, L);
    average(:, :, j) = icfft2(tmp);

    ctf_id=def_grp(:,index); % defocus indices
    metr_list=zeros(size(index));	
    cwf_coeff_nn=cwf_coeff;

    for ii=1:length(cwf_coeff)
	cwf_coeff_nn{ii}=cwf_coeff{ii}(:,index); % Coeffs for nearest neighbors from VDM
	if mod(ii-1,2)==0
		ref_factor=1;
	else
		ref_factor=-1;
	end
	for nn=1:size(cwf_coeff_nn{ii},2)
		if (refl_j(nn)==2)
	cwf_coeff_nn{ii}(:,nn)=conj(cwf_coeff_nn{ii}(:,nn))*ref_factor; % Adjusting coeffs for reflections 
		end	
	end	
        cwf_coeff_nn{ii}=bsxfun(@times,cwf_coeff_nn{ii},exp(-sqrt(-1)*(ii-1)*(-angle_j*pi/180))); % Adjusting coeffs for rotations
    end

    cwf_coeff_nn_pr=cwf_coeff_pr;

    for ii=1:length(cwf_coeff_pr)
	cwf_coeff_nn_pr{ii}=cwf_coeff_pr{ii}(:,index); % Coeffs for nearest neighbors from VDM
	if mod(ii-1,2)==0
		ref_factor=1;
	else
		ref_factor=-1;
	end
	for nn=1:size(cwf_coeff_nn_pr{ii},2)
		if (refl_j(nn)==2)
	cwf_coeff_nn_pr{ii}(:,nn)=conj(cwf_coeff_nn_pr{ii}(:,nn))*ref_factor; % Adjusting coeffs for reflections 
		end	
	end	
        cwf_coeff_nn_pr{ii}=bsxfun(@times,cwf_coeff_nn_pr{ii},exp(-sqrt(-1)*(ii-1)*(-angle_j*pi/180))); % Adjusting coeffs for rotations
    end

    ctfid_j=def_grp(:,j);		
    for jj=1:k
	ctfid_k=ctf_id(:,k);
	t1=1; t2=0;
	for tt=1:length(cwf_coeff_nn_pr)
		t2_add = (cwf_coeff_pr{tt}(:,j)-cwf_coeff_nn_pr{tt}(:,jj))'*pinv(Lmat{ctfid_j}{tt}+Lmat{ctfid_k}{tt})*(cwf_coeff_pr{tt}(:,j)-cwf_coeff_nn_pr{tt}(:,jj));
		assert(imag(t2_add)<1e-10)
		t2 = t2 + t2_add;
		if ((det(Lmat{ctfid_j}{tt}+Lmat{ctfid_k}{tt}))~=0)
			t1=t1*(det(Lmat{ctfid_j}{tt}+Lmat{ctfid_k}{tt}));
		end
	end	
    metr_list(:,jj)=-0.5*log(t1)-0.5*t2;
    end
    
    [sorted_metric,new_nn]=sort(metr_list,'descend');		
    class_met(j,:)=index(new_nn);
    class_refl_met(j,:)=refl_j(:,new_nn);
    angle_met(j,:)=angle_j(:,new_nn);
    shifts(j, :)=-shifts_list(id, 1) - sqrt(-1)*shifts_list(id, 2);

end

% Select top num_nn nearest neighbors using metric
class_met = class_met(:,1:num_nn);
class_refl_met = class_refl_met(:,1:num_nn);
angle_met = angle_met(:,1:num_nn);
delete(gcp)

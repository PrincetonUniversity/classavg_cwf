function [data]=  reduce_coeffs(data_cell, in_freqs, l2_thresh)
reduction=0;

for i=1:length(data_cell)
	l2e=cumsum(sum(abs(data_cell{i}).^2,2)/size(data_cell{i},2));
	if (max(l2e)>0) % Non zero energy in block coeffs
		cutoff=find(l2e/max(l2e)>=l2_thresh,1,'first');
	else
		cutoff=1;% The whole block can be removed, keeping only first row here	
	end
	reduction=reduction+size(data_cell{i},1)-cutoff;
	sprintf('Original size %d, final size %d', size(data_cell{i},1), cutoff)
	data_cell{i}=data_cell{i}(1:cutoff,:);
end
sprintf('Reduced rows by %d', reduction)

size_vec=zeros(1,length(data_cell));
for i=1:length(data_cell)
	size_vec(i)=size(data_cell{i},1);
end

check=(size_vec==0); % Retained at least one row in each block above, so size_vec shouldn't contain 0
assert(max(check(:))~=1);

freqs1=zeros(size(size_vec));
uniq_freq=unique(in_freqs);
so_far=0;
for i=0:max(uniq_freq)
 if size_vec(i+1)~=0
   for j=1+so_far:size_vec(i+1)+so_far
       freqs1(j)=uniq_freq(i+1);
   end
   so_far=so_far+size_vec(i+1);
 end
end


data.Coeff = cell2mat(data_cell);
data.Freqs=freqs1';

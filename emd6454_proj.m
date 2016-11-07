%gen_simulation_data
%Generate 10^4 clean projection images from TRPV1 volume on scratch. Image size
%is 129x129 pixels.
pathstr='/scratch/ARCHIVE_from_sdl6/tbhamre/cwf_class/';

K = 100000; %K is the number of images
max_shift = 0; %maximum shift range
step_size = 1;
initstate;
q = qrand(K);
shifts=round((rand(K,2)-1/2)*2*max_shift/step_size)*step_size;
vol=ReadMRC('/scratch/ARCHIVE_from_sdl6/tbhamre/EMD-6454.map');
vol=cryo_downsample(vol,65);

projections=cryo_project(vol,q,size(vol,1),'single');
data.projections=projections;
data.q=q;
data.shifts=shifts;
%save(fullfile(pathstr,'clean_data_6454'), '-v7.3', 'data')

%downsampleddim=129;
%sprintf('Downsampling to %dX%d grid', downsampleddim, downsampleddim)
%data.projections=cryo_downsample(data.projections,[downsampleddim downsampleddim],1);
%save(fullfile(pathstr,'clean_data_6454_129'), '-v7.3', 'data')

%downsampleddim=65;
%sprintf('Downsampling to %dX%d grid', downsampleddim, downsampleddim)
%data.projections=cryo_downsample(data.projections,[downsampleddim downsampleddim],1);
save(fullfile(pathstr,'clean_data_6454_65_100k'), '-v7.3', 'data')
%clear all;

run ~/aspire/initpath.m
run ~/Dropbox/asinger/80S/initpath_t.m
addpath ~/Dropbox/asinger/80S/80s_Tejal/SUGAR/sugar-master/my_scripts/ccwf_v6/
ndef=10;
def1=1;
def2=4;
lambda = EWavelength(300);
B=10; % decay envelope parameter
n_im=10;
use_CTF=1;
L=209;
hatI_curr=ones(L,L,10);
[g_proj_CTF,CTF,index]=  add_CTF_env_v6(hatI_curr(:,:,1:n_im), ndef, def1,def2,B, lambda, use_CTF);

x=[0:floor(L/2)];
figure(1);
%plot(x,CTF(55,55:end,1),'b-','LineWidth',2, x,CTF(55,55:end,2),'r-','LineWidth',2, x,CTF(55,55:end,3),'g-','LineWidth',2 );
plot(x,CTF(floor(L/2)+1,floor(L/2)+1:end,1),'b-', x,CTF(floor(L/2)+1,floor(L/2)+1:end,2),'r-', x,CTF(floor(L/2)+1,floor(L/2)+1:end,3),'g-','LineWidth',1.5 );
axis tight
legend({'def=1\mum', 'def=1.3\mum', 'def=1.6\mum'},'location','NorthWest', 'Box', 'off','FontSize',9);
xlabel('x','FontSize',11);
ylabel('Contrast','FontSize',11);
print('ctfeg_fig', '-dpng', '-r300'); %<-Save as PNG with 300 DPI

%Iterative Closest Point Differencing Exercise 
%EarthCube Advancing the Analysis of HRT
%Chelsea Scott, Arizona State University, 2018 
clear all;close all
addpath('lasdata')

%Edit the first 5 lines 
pre_dir='Datasets/pre.las';%Pre-earthquake las file
post_dir='Datasets/post.las';%Post-earthquake las file
sz=50;%Differencing window size
grd=250;%Grid spacing; make equal to sz if time allows 
margn=5;%Additional dimension of post-earthquake window. Must be larger than the expect surface displacement

pre=lasdata(pre_dir);%Read the pre-earthquake las file 
pre_x=pre.x;pre_y=pre.y;pre_z=pre.z;%Extract terms from the Matlab structure

post=lasdata(post_dir);%Read the post-earthquake las file
post_x=post.x;post_y=post.y;post_z=post.z;%Extract terms from the Matlab structure

%Construct a core point grid for differencing 
[core_x,core_y]=meshgrid([min(pre_x):grd:max(pre_x)],[min(pre_y):grd:max(pre_y)]);
core_x=core_x(:);core_y=core_y(:);

%ICP for loop 
for i=1:length(core_x)
clear q* p p_*
m=core_x(i);n=core_y(i);
 
% Select points surrounding core point 
a=find(pre_x>m-sz/2&pre_x<m+sz/2&pre_y>n-sz/2&pre_y<n+sz/2);
tz=sz+2*margn;
b=find(post_x>m-tz/2&post_x<m+tz/2&post_y>n-tz/2&post_y<n+tz/2);

%shift (0,0,0) to lie at the center of the grid 
q1=mean(pre_x(a));q2=mean(pre_y(a));q3=mean(pre_z(a));
q_trans(1,:)=pre_x(a)-q1;q_trans(2,:)=pre_y(a)-q2;q_trans(3,:)=pre_z(a)-q3;
p_trans(1,:)=post_x(b)-q1;p_trans(2,:)=post_y(b)-q2;p_trans(3,:)=post_z(b)-q3;

%Perform ICP point-to-plane differencing 
[TR, TT, ER, t] = icp(p_trans,q_trans,'Minimize','plane');

results(i,:) =[core_x(i) core_y(i) TT'];
end

%plot differencing results 
figure
quiver(results(:,1)/1e3,results(:,2)/1e3,results(:,3),results(:,4),'.k','ShowArrowHead','on','LineWidth',3);hold on 
scatter(results(:,1)/1e3,results(:,2)/1e3,45,results(:,5),'filled');
set(gca,'fontsize',14);xlabel('East (km)');ylabel('North (km)');title('ICP displacements');colorbar;colormap(jet)

%%
clear all; clc;





%%

Rt=[12.37	-16.81	-4.31
10.96	-15.64	-5.12
12.74	-15.91	-3.62
12.67	-15.50	-5.29
10.55	-14.83	-3.16
13.71	-17.30	-5.54
12.29	-16.03	-3.64
11.87	-15.87	-8.78
13.98	-13.37	-1.78
14.00	-13.18	-3.30
12.37	-12.18	-5.14
12.66	-13.17	-4.81
11.81	-16.11	-6.65
13.72	-11.16	0.59
14.35	-14.85	-1.99
11.09	-14.63	-6.86
15.33	-11.06	-0.42
12.48	-16.14	-6.52
12.36	-16.11	-5.63
14.07	-15.04	-4.81
14.66	-9.95	2.28
14.32	-15.53	-4.59
14.74	-13.61	5.02
15.41	-10.26	2.41
15.20	-14.05	1.73
13.86	-10.61	-2.26
11.12	-17.31	-5.57
13.74	-14.19	0.01
12.87	-16.89	-2.73];

Lt=[-12.94	-16.12	-4.20
-12.40	-16.21	-3.13
-12.34	-15.64	-4.34
-13.25	-16.78	-7.74
-10.84	-15.86	-3.71
-13.17	-16.87	-4.27
-12.66	-15.82	-2.73
-14.43	-12.33	-2.03
-14.64	-13.36	-1.52
-13.16	-13.05	-4.09
-12.21	-11.81	-5.89
-12.35	-12.05	-4.79
-12.08	-16.96	-7.80
-16.97	-10.84	0.83
-14.39	-16.57	-3.32
-12.99	-17.62	-6.68
-13.80	-15.16	-3.33
-17.43	-11.06	1.22
-15.10	-15.26	-5.22
-15.03	-13.54	-5.26
-17.27	-8.97	5.20
-15.68	-15.51	-1.30
-16.89	-12.12	6.29
-12.86	-16.81	-5.02
-14.02	-13.96	-0.79
-14.04	-14.57	-2.08
-14.13	-14.71	-4.09
-12.50	-15.96	-5.69
-14.52	-16.33	1.56
-13.51	-17.21	-2.37];

All=[Rt;Lt];
All(1:29,1)=All(1:29,1)*-1;

imp=[41.18
77.78
91.67
42.86
27.27
62.5
62.5
62.5
76.92
44.44
33.33
50
33.33
100
81.82
76.47
77.78
52.63
40
73.68
50
38.89
20
-23.08
78.95
25
64.29
0
69.23
50
85.71
76.47
20
37.5
57.14
55.56
50
100
20
55
83.33
81.82
64.29
83.33
60
90.91
70.59
37.5
42.86
-40
71.43
75
81.82
22.22
22.22
77.78
73.33
23.81
69.23];

%%
nii0=ea_load_nii('nii0.nii');

min_x=min(All(:,1));
max_x=max(All(:,1));

min_y=min(All(:,2));
max_y=max(All(:,2));

min_z=min(All(:,3));
max_z=max(All(:,3));
%%
min_xyz=ea_mm2vox([(min_x-3) (min_y-3) (min_z-3)],'/Users/A/Desktop/lead_23/templates/space/MNI_ICBM_2009b_NLIN_ASYM/t1.nii');
max_xyz=ea_mm2vox([(max_x+3) (max_y+3) (max_z+3)],'/Users/A/Desktop/lead_23/templates/space/MNI_ICBM_2009b_NLIN_ASYM/t1.nii');

nii0.img((min_xyz(:,1):max_xyz(:,1)),min_xyz(:,2):max_xyz(:,2),min_xyz(:,3):max_xyz(:,3))=1;
nii0.fname='nii_mask.nii';
ea_write_nii(nii0)
%%
nii=ea_load_nii('vtas.nii');
c=find(nii.img>0);
[x y z]=ind2sub(size(nii.img),c);
a=ea_vox2mm([x y z],'/Users/A/Desktop/lead_23/templates/space/MNI_ICBM_2009b_NLIN_ASYM/t1.nii');

for m=1:length(All)
    for v=1:length(a)
        X=[[All(m,1) All(m,2) All(m,3)];a(v,:)];
        dist2vox(m,v)=pdist(X,'euclidean');  
    end
    disp(m)
end


%%

nii=ea_load_nii('nii_mask.nii');
c=find(nii.img>0);
[x y z]=ind2sub(size(nii.img),c);
a=ea_vox2mm([x y z],'/Users/A/Desktop/lead_23/templates/space/MNI_ICBM_2009b_NLIN_ASYM/t1.nii');

for m=1:length(All)
    for v=1:length(a)
        X=[[All(m,1) All(m,2) All(m,3)];a(v,:)];
        dist2vox(m,v)=pdist(X,'euclidean');  
    end
    disp(m)
end

%%
for i=1:length(a)
     [r(i) p(i)]=corr((dist2vox(:,i)).*-1,imp);
%     [r(i) p(i)]=corr((dist2vox(:,i)).*-1,imp,'type','spearman');
%# do permutations
n_iter = 1000; %# number of permutations
% rperm = zeros(n_iter,1); %# preallocate the vector
% pperm = zeros(n_iter,1); %# preallocate the vector
for k = 1:n_iter
    ind = randperm(59); %# vector of random permutations
    [rperm(i,k) pperm(i,k)] = corr((dist2vox(:,i)).*-1,imp(ind));
end

%# calculate statistics
cc_mean = mean(rperm(i,:));
cc_std = std(rperm(i,:));
zval = r(i) - cc_mean ./ cc_std;
%# probability that the real cc belongs to the same distribution as cc from permuted data
pv(i) = 2 * normcdf(-abs(zval),cc_mean,cc_std); 

disp(i)

end

%%

for pt=1:length(All)
    others=1:length(All);
    others(pt)=[];
  for i=1:length(a)
      dist(others,i)=(dist2vox(others,i)).*-1;
     [r(i) p(i)]=corr(dist(others,i),imp(others));
  end
  for i=1:length(c)
      nii.img(c(i))=r(i); 
  end
    nii.fname=['vox_rmap',num2str(pt),'.nii'];
    ea_write_nii(nii)
    
end

%%
for i=1:59
    
end
%%
for i=1:length(c)
   nii.img(c(i))=r(i); 
end
nii.fname='vox_rmap_nii_mask.nii';
ea_write_nii(nii)

for i=1:length(c)
   nii.img(c(i))=p(i); 
end
nii.fname='vox_pmap_nii_mask.nii';
ea_write_nii(nii)

nii.img(nii.img>0.05)=0;
nii.fname='sig_vox_pmap_nii_mask.nii';
ea_write_nii(nii)

nii.img(nii.img>0.01)=0;
nii.fname='hsig_vox_pmap_nii_mask.nii';
ea_write_nii(nii)

%%
for i=1:length(c)
   nii.img(c(i))=pv(i); 
end
nii.fname='perm_vox_pmap1k_nii_mask.nii';
ea_write_nii(nii);
%%
nii.img(nii.img>0.01)=0;
nii.fname='survived_p_image1k_nii_mask.nii';
ea_write_nii(nii);

%%
cog=[[-14.5000000000000,-16.5000000000000,-3.50000000000000]];
for m=1:length(All)
        COG=[[All(m,1) All(m,2) All(m,3)];cog];
        dist2cog(m)=pdist(COG,'euclidean');  
end
    [r_cog p_cog]=corr(dist2cog',imp)
    
%%
     for pt=1:59
         nii1=ea_load_nii(['vox_rmap',num2str(pt),'.nii']);

    [maxval(pt),maxind(pt)]=max(nii1.img(:));
     [xx(pt) yy(pt) zz(pt)]=ind2sub(size(nii1.img),maxind(pt));
      aa=ea_vox2mm([xx(pt) yy(pt) zz(pt)],'/Users/A/Desktop/lead/templates/space/MNI_ICBM_2009b_NLIN_ASYM/t1.nii');
      XX=[aa;All(pt,:)];
       dist2pred(pt)=pdist(XX,'euclidean');  
     end
     %%
     [r_pred p_pred]=corr(dist2pred',imp)
        %%
cog=[-14.5000000000000,-16.5000000000000,-4];
for m=1:length(All)
        COG=[[All(m,1) All(m,2) All(m,3)];cog];
        dist2cog(m)=pdist(COG,'euclidean');  
end
    [r_cog p_cog]=corr(dist2cog',imp)

%%
dist=(dist2vox(:,:)).*-1;
for i=1:length(a)
     mdl{i}=fitlm([dist(:,i),All(:,1),All(:,2),All(:,3)],imp);
disp(i)

end
%%
for i=1:length(mdl)
   T=anova(mdl{i},'summary');
   
   p(i)=table2array(T(2,5));

end
%%

for i=1:length(c)
   nii.img(c(i))=p(i,1); 
end
nii.fname='vox_pmap_fitlm.nii';
ea_write_nii(nii)


%%
for i=1:length(c)
   nii.img(c(i))=p(i,1); 
end
nii.img(nii.img>0.01)=0;
nii.fname='sig_vox_pmap_fitlm.nii';
ea_write_nii(nii)

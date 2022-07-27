
clear; close all;

%% add DG Laplace coords (not present seince template shape injection was not run!)

% add source/sink
seg = '/export03/data/BigBrain/hippunfold_v1.0.0_noTemplate/work/sub-bbhist/anat/sub-bbhist_hemi-R_space-corobl_dseg.nii.gz';
labelmap = niftiread(seg);
sheader = niftiinfo(seg);
se = strel('sphere',2);
source = imdilate(labelmap==8,se) & imdilate(labelmap==2 | labelmap==4,se) & labelmap==1;
source = imdilate(source,se);
conn = bwconncomp(source);
for i = 1:length(conn.PixelIdxList); e(i) = length(conn.PixelIdxList{i}); end
[~,i] = max(e);
source = false(size(source)); source(conn.PixelIdxList{i}) = true;
labelmap(source & labelmap~=8) = 9;

sink = imdilate(labelmap==8,se) & imdilate(labelmap==0,se) & labelmap==1;
sink = imdilate(sink,se);
e=[];
conn = bwconncomp(sink);
for i = 1:length(conn.PixelIdxList); e(i) = length(conn.PixelIdxList{i}); end
[~,i] = max(e);
sink = false(size(sink)); sink(conn.PixelIdxList{i}) = true;
labelmap(sink & labelmap~=8) = 10;

conn = bwconncomp(labelmap==8);
e=[];
for i = 1:length(conn.PixelIdxList); e(i) = length(conn.PixelIdxList{i}); end
[~,i] = max(e);
labelmap(labelmap==8) = 0;
labelmap(conn.PixelIdxList{i}) = 8;
%niftiwrite(labelmap,'tmp',sheader);

% solve Laplace
se = strel('sphere',25); % HATA has a gap to DG
Laplace_PD = laplace_solver(labelmap==8,labelmap==9,labelmap==10,100);
%Laplace_AP = laplace_solver(labelmap==8,imdilate(labelmap==5,se),imdilate(labelmap==5,se),0);
ap = niftiread('/export03/data/BigBrain/hippunfold_v1.0.0_noTemplate/work/sub-bbhist/coords/sub-bbhist_dir-AP_hemi-R_space-corobl_label-hipp_desc-laplace_coords.nii.gz');
[x,y,z] = ind2sub(size(ap),find(ap>0));
interp='linear';
extrap='nearest';
scattInterp=scatteredInterpolant(x,y,z,double(ap(ap>0)),interp,extrap);
[x,y,z] = ind2sub(size(labelmap),find(labelmap==8));
Laplace_AP = scattInterp(x,y,z);
Laplace_IO = laplace_solver(labelmap==8,(labelmap==1 | labelmap>8),(labelmap~=1 & labelmap<8),10);
Laplace_PD = double(Laplace_PD);
Laplace_AP = double(Laplace_AP);
Laplace_IO = double(Laplace_IO);

% generate surface
[x,y,z] = ind2sub(size(labelmap),find(labelmap==8));
native_coords_mat = [x, y, z, ones(size(x))]';
native_coords_mat = sheader.Transform.T * native_coords_mat;
APres = 254; PDres=30; 
APsamp = [1:APres]/(APres+1);
PDsamp = [1:PDres]/(PDres+1);
IOsamp = 0.5;
[v,u,w] = meshgrid(PDsamp,APsamp,IOsamp); % have to switch AP and PD because matlab sucks
Vuvw = [u(:),v(:),w(:)];
interp='natural';
extrap='nearest';
scattInterp=scatteredInterpolant(Laplace_AP,Laplace_PD,Laplace_IO,native_coords_mat(1,:)',interp,extrap);
x = scattInterp(Vuvw(:,1),Vuvw(:,2),Vuvw(:,3));
scattInterp=scatteredInterpolant(Laplace_AP,Laplace_PD,Laplace_IO,native_coords_mat(2,:)',interp,extrap);
y = scattInterp(Vuvw(:,1),Vuvw(:,2),Vuvw(:,3));
scattInterp=scatteredInterpolant(Laplace_AP,Laplace_PD,Laplace_IO,native_coords_mat(3,:)',interp,extrap);
z = scattInterp(Vuvw(:,1),Vuvw(:,2),Vuvw(:,3));
clear scattInterp
%%
v = [x(:) y(:) z(:)];
v =  v + sheader.Transform.T(4,1:3);
midthick_DG = gifti('/export03/data/opt/hippunfold/hippunfold/resources/unfold_template_dentate/tpl-avg_space-unfold_den-unfoldiso_midthickness.surf.gii');
midthick_DG.vertices = v;

template = '/data/mica3/BigBrain/hippunfold/sub-bbhist/surf/sub-bbhist_hemi-R_space-corobl_den-unfoldiso_label-hipp_midthickness.surf.gii';
midthick_HP = gifti(template);

figure;
p1 = plot_gifti(midthick_DG);
p1.LineStyle = 'none';
p1.FaceColor = 'r';
hold on;
p2 = plot_gifti(midthick_HP);
p2.LineStyle = 'none';
p2.FaceColor = 'b';
light;

%% rotate to corobl and translate

v = midthick_HP.vertices';

sagrot = 40; % degrees
tform = [1 0 0; 
    0 cos(deg2rad(sagrot)) -sin(deg2rad(sagrot)); 
    0 sin(deg2rad(sagrot)) cos(deg2rad(sagrot))]; % 30deg sagittal
v = tform*v;
axrot = 10; % degrees
tform = [cos(deg2rad(axrot)) -sin(deg2rad(axrot)) 0;
    sin(deg2rad(axrot)) cos(deg2rad(axrot)) 0;
    0 0 1]; % 
v = tform*v;
% translate to origin
translat = -min(v');
v = v'+translat;
midthick_HP.vertices = v;


% repeat for DG
v = midthick_DG.vertices';

sagrot = 40; % degrees
tform = [1 0 0; 
    0 cos(deg2rad(sagrot)) -sin(deg2rad(sagrot)); 
    0 sin(deg2rad(sagrot)) cos(deg2rad(sagrot))]; % 30deg sagittal
v = tform*v;
axrot = 10; % degrees
tform = [cos(deg2rad(axrot)) -sin(deg2rad(axrot)) 0;
    sin(deg2rad(axrot)) cos(deg2rad(axrot)) 0;
    0 0 1]; % 
v = tform*v;
% translate to origin
v = v'+translat;
midthick_DG.vertices = v;

%% reconstruct as a smooth version

v = reshape(midthick_HP.vertices,[254,126,3]);
vRec_HP = CosineRep_2Dsurf(v,24,0.01);
midthick_HP.vertices = vRec_HP;

v = reshape(midthick_DG.vertices,[254,30,3]);
vRec_DG = CosineRep_2Dsurf(v,24,0.01);
midthick_DG.vertices = vRec_DG;

figure;
p1 = plot_gifti(midthick_DG);
p1.LineStyle = 'none';
p1.FaceColor = 'r';
hold on;
p2 = plot_gifti(midthick_HP);
p2.LineStyle = 'none';
p2.FaceColor = 'b';
light;

save(midthick_HP,'tpl-avg_space-canonical_den-unfoldiso_label-hipp_midthickness.surf.gii');
save(midthick_DG,'tpl-avg_space-canonical_den-unfoldiso_label-dentate_midthickness.surf.gii');

%% downsample to match various surfaces

densities = {'2mm', '1mm', '0p5mm'};
template32 = gifti('/export03/data/opt/hippunfold/hippunfold/resources/unfold_template_hipp/tpl-avg_space-unfold_den-unfoldiso_midthickness.surf.gii');
for i = 1:length(densities)
    den = densities{i};
    template = gifti(['/export03/data/opt/hippunfold/hippunfold/resources/unfold_template_hipp/tpl-avg_space-unfold_den-' den '_midthickness.surf.gii']);
    ind = dsearchn(template32.vertices,template.vertices);
    template.vertices = vRec_HP(ind,:);
    figure; plot_gifti(template);
    save(template,['tpl-avg_space-canonical_den-' den '_label-hipp_midthickness.surf.gii']);
end

template32 = gifti('/export03/data/opt/hippunfold/hippunfold/resources/unfold_template_dentate/tpl-avg_space-unfold_den-unfoldiso_midthickness.surf.gii');
for i = 1:length(densities)
    den = densities{i};
    template = gifti(['/export03/data/opt/hippunfold/hippunfold/resources/unfold_template_dentate/tpl-avg_space-unfold_den-' den '_midthickness.surf.gii']);
    ind = dsearchn(template32.vertices,template.vertices);
    template.vertices = vRec_DG(ind,:);
    save(template,['tpl-avg_space-canonical_den-' den '_label-dentate_midthickness.surf.gii']);
end

%% patch for https://github.com/khanlab/hippunfold/issues/202

template = gifti(['/export03/data/opt/hippunfold/hippunfold/resources/unfold_template_dentate/tpl-avg_space-unfold_den-2mm_midthickness.surf.gii']);

figure;
p1 = plot_gifti(template);
p1.LineStyle = 'none';
p1.FaceColor = 'r';

i = find(template.vertices(:,2) > -195.3);
[m,ii] = min(template.vertices(i,1));
template.vertices(i(ii),1) = min(template.vertices(:,1));
[m,ii] = max(template.vertices(i,1));
template.vertices(i(ii),1) = max(template.vertices(:,1));

figure;
p1 = plot_gifti(template);
p1.LineStyle = 'none';
p1.FaceColor = 'r';

save(template,['/export03/data/opt/hippunfold/hippunfold/resources/unfold_template_dentate/tpl-avg_space-unfold_den-2mm_midthickness.surf.gii']);
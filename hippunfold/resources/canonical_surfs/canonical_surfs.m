
clear; close all;

hippunfolddir = '/data/mica3/BIDS_BigBrain/derivatives/hippunfold_v0.6.0/';
sub = 'bbhist';
hemi = 'R';
midthick = gifti([hippunfolddir '/results/sub-' sub '/surf_cropseg/sub-' sub '_hemi-' hemi '_space-corobl_den-32k_midthickness.surf.gii']);

%% rotate to corobl and translate

sagrot = 40; % degrees
tform = [1 0 0; 
    0 cos(deg2rad(sagrot)) -sin(deg2rad(sagrot)); 
    0 sin(deg2rad(sagrot)) cos(deg2rad(sagrot))]; % 30deg sagittal
v2 = midthick;
v2.vertices = tform*v2.vertices';
v2.vertices = v2.vertices';

axrot = 10; % degrees
tform = [cos(deg2rad(axrot)) -sin(deg2rad(axrot)) 0;
    sin(deg2rad(axrot)) cos(deg2rad(axrot)) 0;
    0 0 1]; % 
v2.vertices = tform*v2.vertices';
v2.vertices = v2.vertices';


% translate to origin
v2.vertices = v2.vertices-min(v2.vertices);

midthick.vertices = v2.vertices;

%% reconstruct as a smooth version

vertices = reshape(midthick.vertices,[254,126,3]);
vRec = CosineRep_2Dsurf(vertices,24,0.01);

midthick.vertices = vRec;
plot_gifti(midthick);
save(midthick,'tpl-avg_space-canonical_den-32k_midthickness.surf.gii');

%% downsample to match various surfaces

resources = '/host/percy/local_raid/jordand/opt/hippunfold/hippunfold/resources/unfold_template/';
template32 = gifti([resources '/tpl-avg_space-unfold_den-32k_midthickness.surf.gii']);
desnities = {'400', '2k', '7k'};
for d = 1:length(desnities)
    den = desnities{d};
    template = gifti([resources '/tpl-avg_space-unfold_den-' den '_midthickness.surf.gii']);
    ind = dsearchn(template32.vertices,template.vertices);
    template.vertices = vRec(ind,:);
    figure; plot_gifti(template);
    save(template,['tpl-avg_space-canonical_den-' den '_midthickness.surf.gii']);
end

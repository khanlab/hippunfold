LBL = niftiread('tpl-upenn_desc-hipptissue_dseg.nii.gz');
hdr = niftiinfo('tpl-upenn_dir-AP_coords.nii.gz');

idxgm = (LBL==1 | LBL==8);
AP = laplace_solver(idxgm,LBL==5,LBL==6,10000);
PD = laplace_solver(idxgm,LBL==3,LBL==8,10000);
IO = laplace_solver(idxgm,(LBL==2 | LBL==4 | LBL==7),LBL==0,10000);
ap = zeros(size(LBL)); ap(idxgm) = AP;
pd = zeros(size(LBL)); pd(idxgm) = PD;
io = zeros(size(LBL)); io(idxgm) = IO;
niftiwrite(single(ap),'tpl-upenn_dir-AP_coords',hdr,'compressed',true);
niftiwrite(single(pd),'tpl-upenn_dir-PD_coords',hdr,'compressed',true);
niftiwrite(single(io),'tpl-upenn_dir-IO_coords',hdr,'compressed',true);


LBLDG = niftiread('tpl-upenn_desc-DGtissue_dseg.nii.gz');
ap = zeros(size(LBL));
ap(idxgm) = AP;
AP = ap(LBLDG==8);
PD = laplace_solver(LBLDG==8,LBLDG==1,LBLDG==2,5000);
IO = laplace_solver(LBLDG==8,LBL==1,(LBL==2 | LBL==7 | LBL==0),1000);
ap = zeros(size(LBL)); ap(LBLDG==8) = AP;
pd = zeros(size(LBL)); pd(LBLDG==8) = PD;
io = zeros(size(LBL)); io(LBLDG==8) = IO;
niftiwrite(single(ap),'tpl-upenn_dir-AP_coords-DG',hdr,'compressed',true);
niftiwrite(single(pd),'tpl-upenn_dir-PD_coords-DG',hdr,'compressed',true);
niftiwrite(single(io),'tpl-upenn_dir-IO_coords-DG',hdr,'compressed',true);

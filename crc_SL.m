%
% Small script to hack through some PRT structure (from PRoNTo) in order to
% proceed with a "search light" approach.
%
% The idea is simple:
% 1/ proceed with your whole brain PRoNTo classification
% 2/ use the full prt structure to as a template for "search light" and
%    loop over the voxels in the 1st level mask
%       - create a small spherical volume
%       - re-estimate the kernel from the subset of voxels
%       - launch the model estimation
%       - collect tha accuracies
%__________________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by Christophe Phillips
% University of Liège, Belgium


R = 10; % search light radius in mm
Pprt = spm_select(1,'mat','Select the PRT.mat file');
[pth,nam,ext,num] = spm_fileparts(Pprt);
load(Pprt)

i_model = 1; % Model index to use -> adjust according to your need!!!
PRTorig = PRT;

Pmsk = PRT.masks.fname ; % should be the updated 1st level mask
Vmsk = spm_vol(Pmsk);

%-Get space details
M         = Vmsk.mat;                          %-voxels to mm matrix
iM        = inv(M);                             %-mm to voxels matrix
DIM       = Vmsk.dim;                          %-image dimensions
[x,y,z]   = ndgrid(1:DIM(1),1:DIM(2),1:DIM(3));
XYZ       = [x(:)';y(:)';z(:)']; clear x y z    %-voxel coordinates {vx}
XYZmm     = M(1:3,:)*[XYZ; ones(1,size(XYZ,2))];%-voxel coordinates {mm}
XYZmm_cpy = XYZmm;                              %-copy without masking

% List of mask voxels>1
lVx = PRT.fas.idfeat_img;
nVx = numel(lVx);

%-Search volume (from mask)
MM     = spm_read_vols(Vmsk);
MM     = logical(MM);
XYZmm  = XYZmm(:,MM(:));
XYZ    = XYZ(:,MM(:));

% sanity check
ll = find(MM(:));
any(lVx-ll)

%-Searchlight options (clique definition)
searchopt.def = 'sphere';
searchopt.spec = R;

%-Build local clique
xY      = searchopt;
xY.xyz  = [NaN NaN NaN];
xY.rej  = {'cluster','mask'};
xY      = spm_ROI(xY);

%-Build local clique, in voxels
c            = round(DIM(:)/2);
xY.xyz       = M(1:3,:) * [c;1];
[xY, clique] = spm_ROI(xY,XYZmm_cpy);
clique       = round(iM(1:3,:) * [clique;ones(1,size(clique,2))]);
clique       = bsxfun(@minus, clique, c);
dc           = (max(clique,[],2) - min(clique,[],2) + 1)'; % size of clique radii, in voxels

%% Get the input ready before looping
PRTw = rmfield(PRTorig,'model'); % PRT before model estimation, still need to adjust feature selection!
PRTw.model = rmfield(PRTorig.model,'output');

% Use filearray
% fs_whole = PRTw.fas.dat;
% or load in memory
fs_whole = PRTw.fas.dat(:,:);

%% Loops over all voxels from the 1st level mask & collect accuracies
SLres(nVx+1) = PRTorig.model(i_model).output.stats; % initialier SL results structure array and include at N+1 the full mask results
kern_file = fullfile(pth,PRTorig.fs.fs_name);
h_wb = waitbar(0,'Voxel counts');

tic
for ii=1:nVx
    PRT = PRTw;
    i_vx = lVx(ii);
    % 1/ update the list of voxels for 2nd level mask in PRT
%     xyz_cl = bsxfun(@plus,XYZ(:,i_vx),clique);
    xyz_cl = bsxfun(@plus,XYZ(:,ii),clique);
    % remove stuff outside image volume
    l_2rm = logical(sum(bsxfun(@lt, xyz_cl, ones(3,1))));
    xyz_cl(:,l_2rm) = [];
    l_2rm = logical(sum(bsxfun(@gt, xyz_cl, DIM')));
    xyz_cl(:,l_2rm) = [];
    % create list
    Lcl = xyz_cl(1,:) + (xyz_cl(2,:)-1)*DIM(1) + (xyz_cl(3,:)-1)*DIM(1)*DIM(2);
    Lres = find(any(bsxfun(@eq, lVx', Lcl')));
    % List of voxels to use
    PRT.fs.modality.idfeat_fas = Lres';
    
    % 2/ Rebuild kernel & save it
    datSL = fs_whole(:,Lres);
    Phim = datSL*datSL';
    [~,idmax] = max(Phim);
    [~,idmin] = min(Phim);
    min_max = find(idmax==idmin);
    if isempty(min_max) || unique(Phim(:,min_max))~=0 %Kernel does not contain a whole line of zeros
        igd = true;
        Phi{1} = Phim;
        save(kern_file,'Phi');
    else
        igd = false;
        fprintf('\n Skipping voxel %d.',ii);
    end
    
    if igd
        % 3/ launch estimate
        %     PRT = prt_model(PRT,mod_w);
        in.fname      = Pprt;
        in.model_name = PRT.model(i_model).model_name;
        in.savePRT = false;
        PRT = prt_cv_model(PRT, in);
        % 4/ collect accuracies
        SLres(ii) = PRT.model(i_model).output.stats;
    end
    clear PRT
    waitbar(ii/nVx,h_wb)
end
toc

close(h_wb)

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
%
% NOTE: 
% PRoNTo must be initialized before using the script as it relies on
% PRoNTo's machinery.
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
for ivx=1:nVx
    PRT = PRTw;
    i_vx = lVx(ivx);
    % 1/ update the list of voxels for 2nd level mask in PRT
    xyz_cl = bsxfun(@plus,XYZ(:,ivx),clique);
    % remove stuff outside image volume
    l_2rm = logical(sum(bsxfun(@lt, xyz_cl, ones(3,1))));
    xyz_cl(:,l_2rm) = [];
    l_2rm = logical(sum(bsxfun(@gt, xyz_cl, DIM')));
    xyz_cl(:,l_2rm) = [];
    % create list
    Lcl = sub2ind(DIM,xyz_cl(1,:),xyz_cl(2,:),xyz_cl(3,:));
    Lres = find(any(bsxfun(@eq, lVx', Lcl')));
    % List of voxels to use
    PRT.fs.modality.idfeat_fas = lVx(Lres);
    
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
        fprintf('\n Skipping voxel %d.',ivx);
    end
    
    if igd
        % 3/ launch estimate
        %     PRT = prt_model(PRT,mod_w);
        in.fname      = Pprt;
        in.model_name = PRT.model(i_model).model_name;
        in.savePRT = false;
        PRT = prt_cv_model(PRT, in);
        % 4/ collect accuracies
        SLres(ivx) = PRT.model(i_model).output.stats;
        
        % TO DO
        % 5/ permuation testing to get a p-value
        % -> only do this if it's worth (check the accuracy) !
    end
    clear PRT
    waitbar(ivx/nVx,h_wb)
end
toc
close(h_wb)

%% Save results into images.
% Results structure for classification & regression
%   Classification:
% stats.con_mat: Confusion matrix (nClasses x nClasses matrix, pred x true)
% stats.acc:     Accuracy (scalar)
% stats.b_acc:   Balanced accuracy (nClasses x 1 vector)
% stats.c_acc:   Accuracy by class (nClasses x 1 vector)
% stats.c_pv:    Predictive value for each class (nClasses x 1 vector)
% stats.acc_lb:  \_ Lower/upper 5% confidence interval bounds for a  
% stats.acc_ub:  /  binomial distribution using Wilson's 'score interval'
%
%   Regression:
% stats.mse:     Mean square error between test and prediction
% stats.corr:    Correlation between test and prediction
% stats.r2:      Squared correlation

% TO DO, when the p-value is available, save it as an image too!

switch PRT.model(i_model).input.type
    case 'classification'
        % save acc/b_acc/c_acc only, maybe add others or leave option later on...
        % 1/ prepre files: acc, b_acc, 
        Vacc = Vmsk;
        Vacc.fname = fullfile(pth,['SLacc_',PRT.model(i_model).model_name,'.nii']);
        Vacc.dt(1) = 16; % save as float
        Vacc.descrip = 'PRoNTo Search Light accuracy';
        Vacc = spm_create_vol(Vacc);
        Vbacc = Vacc;
        Vbacc.fname = fullfile(pth,['SLbacc_',PRT.model(i_model).model_name,'.nii']);
        Vbacc.descrip = 'PRoNTo Search Light balanced accuracy';
        Vbacc = spm_create_vol(Vbacc);
        nClasses = numel(SLres(end).c_acc);
        Vcacc(nClasses) = Vacc;
        for ii=1:nClasses
            Vcacc(ii) = Vacc;
            Vcacc(ii).fname = fullfile(pth,['SLcacc',num2str(ii),'_',PRT.model(i_model).model_name,'.nii']);
            Vcacc(ii).descrip = ['PRoNTo Search Light class ',num2str(ii),' accuracy'];
            Vcacc(ii) = spm_create_vol(Vcacc(ii));
        end
        % 2/ prepare & collect values
        val_acc = zeros(prod(DIM),1)+NaN;
        val_bacc = zeros(prod(DIM),1)+NaN;
        val_cacc = zeros(prod(DIM),nClasses)+NaN;
        for ivx=1:nVx
            if ~isempty(SLres(ivx).acc)
                i_vx = lVx(ivx);
                val_acc(i_vx) = SLres(ivx).acc;
                val_bacc(i_vx) = SLres(ivx).b_acc;
                val_cacc(i_vx,:) = SLres(ivx).c_acc;
            end
        end
        val_acc = reshape(val_acc,DIM);
        val_bacc = reshape(val_bacc,DIM);
        val_cacc = reshape(val_cacc,[DIM nClasses]);
        % 3/ save values
        Vacc = spm_write_vol(Vacc,val_acc);
        Vbacc = spm_write_vol(Vbacc,val_bacc);
        for ii=1:nClasses
            Vcacc(ii) = spm_write_vol(Vcacc(ii),squeeze(val_cacc(:,:,:,ii)));
        end
    case 'regression'
        % Nothing here yet!
    otherwise
        fprintf('\nUnknown model type!\n')
end

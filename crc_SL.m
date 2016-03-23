function [SLres,Pout] = crc_SL(Pprt,opt)

% Function to hack through some PRT structure (from PRoNTo) in order to
% proceed with a "search light" approach.
%
% The idea is simple:
% 1/ proceed with your whole brain PRoNTo classification
% 2/ use the full prt structure as a template for "search light" and
%    loop over the voxels in the *1st level mask*
%       - create a small spherical volume
%       - re-estimate the kernel from the subset of voxels
%       - launch the model estimation
%       - collect tha accuracies
% 3/ save accuracies (full, balanced & per class) into images
%
% FORMAT
% [SLres,Pout] = crc_SL(Pprt,opt)
%
% INPUT:
% - Pprt    Filename (with path) of PRT.mat file to use
% - opt     some options (more could be added...)
%       .R        radius in mm of the spherical searchlight [10 by def.]
%       .i_model  index of PRoNTo model to use [1 by def.]
%       .loadF    load all features or not [true by def.]
%       .savImg   save results in image format or not [true by def.]
%       .permStat assess accuracy through permuation [false by def.]
%
% OUTPUT:
% - SLres   Searchlight results structure array. There is 1 structure per
%           voxel in the mask + 1 (last one) with the original whole mask
%           results
% - Pout    Filenames of generated stuff
%
%--------------------------------------------------------------------------
% NOTE:
% - PRoNTo must be initialized before using the script as it relies on
%   PRoNTo's machinery.
% - other clique format could be used, for example: cubes in images or
%   planes/lines for other types of data (connectivity matrix, ERP's, etc.)
% - no inference at the moment but this could be added...
%__________________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by Christophe Phillips
% University of Liège, Belgium

% Initialize PRoNTo if needed
if ~exist('prt_checkAlphaNumUnder','file')
    prt('startup','nogui')
end

% Deal with inputs
opt_def = struct( ...
    'R',10, ...
    'i_model',1, ...
    'loadF',true, ...
    'savImg',true, ...
    'permStat',false);
if nargin<2
    opt = opt_def;
else
    opt = prt_check_flag(opt_def,opt);
end
R = opt.R; % search light radius in mm
i_model = opt.i_model; % Model index to use
loadF = opt.loadF;
savImg = opt.savImg;
permStat = opt.permStat;

if nargin<1 || isempty(Pprt)
    Pprt = spm_select(1,'mat','Select the PRT.mat file');
end
[pth,nam,ext,num] = spm_fileparts(Pprt); %#ok<*NASGU,*ASGLU>
load(Pprt)
PRTorig = PRT; %#ok<*NODEF>

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
MM(isnan(MM)) = 0; % turn NaN's into 0's
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

if loadF
    % Load all features in memory
    fs_whole = PRTw.fas.dat(:,:);
else
    % or use filearray -> slower but less memory hungry!
    fs_whole = PRTw.fas.dat;
end

% find feature set index: i_fs
fs_names = cellstr(char(PRTorig.fs(:).fs_name));
i_fs = find(strcmp(PRTorig.model(i_model).input.fs.fs_name,fs_names));
if isempty(i_fs)
    error('prtSL:fs_model','No matching feature set for model #%d.\n',i_model)
end

%% Loops over all voxels from the 1st level mask & collect accuracies
% initialier SL results structure array and
% include at N+1 the full mask results
tmp_stats = PRTorig.model(i_model).output.stats;
if ~permStat && isfield(PRTorig.model(i_model).output.stats,'permutation')
    tmp_stats = rmfield(tmp_stats,'permutation');
end
SLres(nVx+1) = tmp_stats;
kern_file = fullfile(pth,PRTorig.fs(i_fs).k_file);
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
    PRT.fs(i_fs).modality.idfeat_fas = lVx(Lres);
    
    % 2/ Rebuild kernel & save it
    datSL = fs_whole(:,Lres);
    Phim = datSL*datSL';
    [~,idmax] = max(Phim);
    [~,idmin] = min(Phim);
    min_max = find(idmax==idmin);
    if isempty(min_max) || unique(Phim(:,min_max))~=0 %Kernel does not contain a whole line of zeros
        igd = true; % -> proceed with model estimation
        Phi{1} = Phim;
        save(kern_file,'Phi');
    else
        igd = false; % -> no need to estimate the model
        fprintf('\n Skipping voxel %d.',ivx);
    end
    
    if igd % Proceed if the index is "good"
        % 3/ launch estimate
        %     PRT = prt_model(PRT,mod_w);
        in.fname      = Pprt;
        in.model_name = PRT.model(i_model).model_name;
        in.savePRT = false;
        PRT = prt_cv_model(PRT, in);
        SLres_ii = PRT.model(i_model).output.stats;

        % TO DO
        % 4/ permuation testing to get a p-value
        % -> only do this if it's worth (check the accuracy) !
        if permStat
            % do the stats
            SLres_ii.permutation = [];
        end
        
        % 5/ collect results
        SLres(ivx) = SLres_ii;
    end
    clear PRT
    waitbar(ivx/nVx,h_wb)
end
toc
close(h_wb)

P_SLresults = spm_file(Pprt,'filename','SLresults.mat');
save(P_SLresults,'SLres')

Pout{1} = P_SLresults;

%% Save results into images, if requested
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
if savImg
    switch PRTw.model(i_model).input.type
        case 'classification'
            % save acc/bacc/cacc only, maybe add others or leave option later on...
            % 1/ prepre files: Vacc, Vbacc, Vcacc(1/2/...)
            Vacc = Vmsk;
            Vacc.fname = fullfile(pth,['SLacc_',PRTw.model(i_model).model_name,'.nii']);
            Vacc.dt(1) = 16; % save as float
            Vacc.descrip = 'PRoNTo Search Light accuracy';
            Vacc = spm_create_vol(Vacc);
            Vbacc = Vacc;
            Vbacc.fname = fullfile(pth,['SLbacc_',PRTw.model(i_model).model_name,'.nii']);
            Vbacc.descrip = 'PRoNTo Search Light balanced accuracy';
            Vbacc = spm_create_vol(Vbacc);
            nClasses = numel(SLres(end).c_acc);
            Vcacc(nClasses) = Vacc;
            for ii=1:nClasses
                Vcacc(ii) = Vacc;
                Vcacc(ii).fname = fullfile(pth,['SLcacc',num2str(ii),'_',PRTw.model(i_model).model_name,'.nii']);
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
            % 3/ save values into images
            Vacc = spm_write_vol(Vacc,val_acc);
            Vbacc = spm_write_vol(Vbacc,val_bacc);
            for ii=1:nClasses
                Vcacc(ii) = spm_write_vol(Vcacc(ii),squeeze(val_cacc(:,:,:,ii)));
            end
            % 4/ save filenames
            Pout{2} = Vacc.fname;
            Pout{3} = Vbacc.fname;
            for ii=1:nClasses
                Pout{3+ii} = Vcacc(ii).fname; %#ok<*AGROW>
            end
        case 'regression'
            % save mse, corr, r2
            % 1/ prepre files: Vmse, Vcorr, Vr2
            Vmse = Vmsk;
            Vmse.fname = fullfile(pth,['SLmse_',PRTw.model(i_model).model_name,'.nii']);
            Vmse.dt(1) = 16; % save as float
            Vmse.descrip = 'PRoNTo Search Light MSE';
            Vmse = spm_create_vol(Vmse);
            Vcorr = Vmse;
            Vcorr.fname = fullfile(pth,['SLcorr_',PRTw.model(i_model).model_name,'.nii']);
            Vcorr.descrip = 'PRoNTo Search Light Correlation';
            Vcorr = spm_create_vol(Vcorr);
            Vr2 = Vmse;
            Vr2.fname = fullfile(pth,['SLr2_',PRTw.model(i_model).model_name,'.nii']);
            Vr2.descrip = 'PRoNTo Search Light R2';
            Vr2 = spm_create_vol(Vr2);
            
            % 2/ prepare & collect values
            val_mse = zeros(prod(DIM),1)+NaN;
            val_corr = zeros(prod(DIM),1)+NaN;
            val_r2 = zeros(prod(DIM),1)+NaN;
            for ivx=1:nVx
                if ~isempty(SLres(ivx).mse)
                    i_vx = lVx(ivx);
                    val_mse(i_vx) = SLres(ivx).mse;
                    val_corr(i_vx) = SLres(ivx).corr;
                    val_r2(i_vx) = SLres(ivx).r2;
                end
            end
            val_mse = reshape(val_mse,DIM);
            val_corr = reshape(val_corr,DIM);
            val_r2 = reshape(val_r2,DIM);
            % 3/ save values
            Vmse = spm_write_vol(Vmse,val_mse);
            Vcorr = spm_write_vol(Vcorr,val_corr);
            Vr2 = spm_write_vol(Vr2,val_r2);
            Pout{2} = Vmse.fname;
            Pout{3} = Vcorr.fname;
            Pout{4} = Vr2.fname;
        otherwise
            fprintf('\nUnknown model type!\n')
    end
end

end
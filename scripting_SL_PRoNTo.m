% Test variation of accuracy with region size.
% Usee radii of increasing size: 5, 10, 15, 20, 25

Rlist = 5:5:25;
Pin = 'C:\Dropbox\Work\3_data\PRoNTo\Trav_SL\results5_allSt_br\PRT.mat';
for ii = 1:numel(Rlist)
    opt.R = Rlist(ii);
    [SLres{ii},Pout{ii}] = crc_parSL(Pin,opt); %#ok<*SAGROW>
end


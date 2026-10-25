function [sweeps, seeds, schemes, metrics] = sweep_defs()
%SWEEP_DEFS  Sweeps, seeds, schemes and metrics shared by all run scripts.
seeds   = 1:10;
schemes = {'LAMS-noAlt','LAMS-noAlt-noCal','LAMS','TH','RC','MH','EPI','SNW'};
metrics = {'pdr','delayMean','delayP95','txPerDel','ctrlPerDel','occMean', ...
           'shareSCF','shareAlt','switchPerMin','dropTTL','dropBuf','dropNoRoute'};
% name, parameter field, values
sweeps = {
    'N',      'N',      [20 30 40 50 60];
    'speed',  'vMax',   [15 25 35];
    'load',   'lambda', [0.5 1 2 4];
    'shadow', 'sigmaS', [2 4 8];
};
end

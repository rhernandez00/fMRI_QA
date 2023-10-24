clear

thr_fwd = 0.5; %fwd threshold
thr_dvars = 0.5; %dvars threshold
tableName = [pwd,'/table.xlsx']; %path of the output table 
filesPath = [pwd,'/input']; %folder of files to analyze
savePath = [pwd,'/txt']; %output folder


cfg.radius = 50; %some option needed for fwd, we never change this
fileList = {dir([filesPath,'/*.*']).name};
fileList(1:2) = [];

primaryKeyList = cell(1,size(fileList,2));
for nFile = 1:size(fileList,2)
    fileName = fileList{nFile};
    dot = strfind(fileName,'.');
    dot = dot(end);
    extention = fileName(dot+1:end);
    switch extention
        case 'par' %This is a movement file. Running fwd
            primaryKey = fileName(1:end-4);
        case 'txt' %This is a movement file. Running fwd
            primaryKey = fileName(1:end-4);
        case 'nii' %This is a nii file. Running dvars
            primaryKey = fileName(1:end-4);
        case 'gz' %This is a nii file. Running dvars
            primaryKey = fileName(1:end-7);
        otherwise
            error(['Unrecognized extention: ', extention, ' for file: ',fileName]);
    end
    primaryKeyList{nFile} = primaryKey;
end

keyList = unique(primaryKeyList)';
totalKeys = length(keyList);
preproSuiteList = cell(numel(keyList),1);
preproSuiteList(:) ={'unkown'};
tableOut = table(keyList,preproSuiteList,zeros(size(keyList)),...
    zeros(size(keyList)),zeros(size(keyList)),...
    'VariableNames',{'primaryKey','PreproSuite','fwd','fwd_clean','volumesLostAfter_fwd'});

%This runs fwd for each primaryKey if movement files are found
indxfwd = cell(1,numel(keyList));
for nKey = 1:numel(keyList) 
    primaryKey = keyList{nKey};
    if exist([filesPath,'/',primaryKey,'.txt'],'file') && exist([filesPath,'/',primaryKey,'.par'],'file')
        error(['There are txt and par files for: ',primaryKey,' there should only be one']);
    elseif exist([filesPath,'/',primaryKey,'.txt'],'file')
        PreproSuite = 'spm';
        cfg.prepro_suite =  'spm';
        cfg.motionparam = [filesPath,'/',primaryKey,'.txt'];
    elseif exist([filesPath,'/',primaryKey,'.par'],'file')
        PreproSuite = 'fsl';
        cfg.prepro_suite =  'fsl-fs';
        cfg.motionparam = [filesPath,'/',primaryKey,'.par'];
    else %no movement files for that primaryKey
        disp(['no movement files for: ', primaryKey]);
        continue
    end
    
    fwd = bramila_framewiseDisplacement(cfg);
    indx = fwd > thr_fwd;
    indxNums = find(indx);
    totalVolumes = numel(indx);
    lostVolumes = numel(indxNums);
    tableOut.fwd(nKey) = mean(fwd);
    tableOut.fwd_clean(nKey) = mean(fwd(fwd < thr_fwd));
    tableOut.lostVolumesAfterfwd(nKey) = lostVolumes;
    
    tableOut.totalVolumes(nKey) = totalVolumes;
    

    tableOut.PreproSuite{nKey} = PreproSuite;
    indxfwd{nKey} = indx;
end
clear cfg

%This runs dvars
indxDvars = cell(1,numel(keyList));

for nKey = 1:numel(keyList)
    primaryKey = keyList{nKey};
    if exist([filesPath,'/',primaryKey,'.nii.gz'],'file') && exist([filesPath,'/',primaryKey,'.nii'],'file')
        error(['There are two nifti files for: ',primaryKey,' there should only be one']);
    elseif exist([filesPath,'/',primaryKey,'.nii'],'file')
        vol = load_untouch_nii([filesPath,'/',primaryKey,'.nii']);
    elseif exist([filesPath,'/',primaryKey,'.nii.gz'],'file')
        vol = load_untouch_nii([filesPath,'/',primaryKey,'.nii.gz']);
    else
        continue
    end
    cfg.vol = vol.img; %assigns the nii file
    [dvars]=bramila_dvars(cfg); %runs dvars
    dvarsZ = abs([0;zscore(dvars(2:end))]);%transforms dvars to zscore
    indx = dvarsZ > thr_dvars; %
    indxNums = find(indx);
    totalVolumes = numel(indx);
    lostVolumes = length(indxNums);
    if tableOut.totalVolumes(nKey) == 0
        tableOut.totalVolumes(nKey) = totalVolumes;
    else
        if tableOut.totalVolumes(nKey) ~= totalVolumes
            error(['The volumes in the movement file and the nifti file',...
                'do not match for: ', primaryKey, ' volumes in nii: ',...
                num2str(totalVolumes), ' volumes in movement file: ',...
                num2str(tableOut.totalVolumes(nKey))]);
        end
    end
    tableOut.dvars(nKey) = mean(dvarsZ);
    tableOut.dvars_clean(nKey) = mean(dvarsZ(dvarsZ < thr_fwd));
    tableOut.lostVolumesAfterdvars(nKey) = lostVolumes;
    indxDvars{nKey} = indx;
end

for nKey = 1:numel(keyList)
    primaryKey = keyList{nKey};
    if isempty(indxDvars{nKey})
        indxDvars{nKey} = zeros(tableOut.totalVolumes(nKey),1);
    end
    if isempty(indxfwd{nKey})
        indxfwd{nKey} = zeros(tableOut.totalVolumes(nKey),1);
    end
    indx1 = indxfwd{nKey};
    indx2 = indxDvars{nKey};
    indx = sum([indx1,indx2],2)>0;
    indx = logical(indx);
    writeTxt([savePath,'/',primaryKey,'.txt'],indx); %
    
end


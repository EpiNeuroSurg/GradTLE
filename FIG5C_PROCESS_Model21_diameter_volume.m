%% S16 PROCESS Model 21: diameter-volume vascular burden
% X = raw ipsilateral mean diameter x raw ipsilateral volume, or
%     diameter asymmetry x volume asymmetry; W = duration;
% M = -1*thickness; V = -1*MidFLAIR; Y = PET.
% Territories: ACTP, MCAT, PCAT, and an equal-weight ALL composite
% Flat script: no local functions and no file/figure export.

clear; clc; close all;
rng(20260914,'twister');

%% Settings
rootDir = '/host/verges/tank/data/mo/BIDS_20251201_TLE_BTTH1';
imagingFile = '/host/verges/tank/data/mo/demo/02Voxel_feats_hmz_flip_asym_wPDB.mat';
maskDir = '/host/verges/tank/data/mo/atlas/arterial/fsLR32k';
nBoot = 5000;

territoryNames = {'ACTP','MCAT','PCAT','ALL'};
maskFiles = {
    fullfile(maskDir,'ACTPL_hemi-L_fsLR32k.func.gii')
    fullfile(maskDir,'MCATL_hemi-L_fsLR32k.func.gii')
    fullfile(maskDir,'PCATL_hemi-L_fsLR32k.func.gii')};
leftLabels  = {'LAChA','LMCA1','LPCA1'};
rightLabels = {'RAChA','RMCA1','RPCA1'};

combinationNames = {
    'Raw diameter x raw volume'
    'Diameter asymmetry x volume asymmetry'};

%% 1. Clinical variables and exact imaging row order
run('/host/verges/tank/data/mo/load_demo.m');

imagingID = upper(strtrim(string(Demo_disPTs.ID)));
imagingID = regexprep(imagingID,'^SUB-','');
duration227 = double(Demo_disPTs.Duration);
side227 = upper(strtrim(string(Demo_disPTs.Side)));
side227(startsWith(side227,'LEFT')) = "L";
side227(startsWith(side227,'RIGHT')) = "R";

if numel(imagingID) ~= 227
    error('Demo_disPTs contains %d subjects; expected 227.',numel(imagingID));
end

%% 2. Extract three mediators and PET from arterial territories
S = load(imagingFile,'Voxel_PTs_neofeats_hmz_flip_asym_wPDB');
D = S.Voxel_PTs_neofeats_hmz_flip_asym_wPDB;

CorThick = double(D.CorThick);
MidFLAIR = double(D.MidFLAIR);
Myelin = double(D.Myelin);
MidPET = double(D.MidPET);

if any([size(CorThick,1),size(MidFLAIR,1),size(Myelin,1),size(MidPET,1)] ~= 227)
    error('CorThick, MidFLAIR, Myelin, and MidPET must each contain 227 rows.');
end
if any([size(CorThick,2),size(MidFLAIR,2),size(Myelin,2),size(MidPET,2)] < 32492)
    error('An imaging matrix has fewer than 32,492 cortical vertices.');
end

% First 32,492 vertices are ipsilateral after hemisphere flipping.
CorThick = CorThick(:,1:32492);
MidFLAIR = MidFLAIR(:,1:32492);
Myelin = Myelin(:,1:32492);
MidPET = MidPET(:,1:32492);

thickness227 = nan(227,3);
flair227 = nan(227,3);
myelin227 = nan(227,3);
pet227 = nan(227,3);

for t = 1:3
    G = gifti(maskFiles{t});
    maskValue = double(G.cdata(:));
    if numel(maskValue) ~= 32492
        error('%s contains %d vertices; expected 32,492.',maskFiles{t},numel(maskValue));
    end
    idx = isfinite(maskValue) & maskValue > 0;
    thickness227(:,t) = mean(CorThick(:,idx),2,'omitnan');
    flair227(:,t) = mean(MidFLAIR(:,idx),2,'omitnan');
    myelin227(:,t) = mean(Myelin(:,idx),2,'omitnan');
    pet227(:,t) = mean(MidPET(:,idx),2,'omitnan');
end

clear CorThick MidFLAIR Myelin MidPET D S G maskValue idx

%% 3. Extract bilateral raw volumes and calculate ipsilateral asymmetry
eFiles = dir(fullfile(rootDir,'derivatives','eICAB','sub-*','ses-pre', ...
    'anat','metric_space','*_CircleOfWillis','*_SUMMARY.csv'));
if isempty(eFiles)
    error('No eICAB *_SUMMARY.csv files were found.');
end

eicabID = strings(numel(eFiles),1);
ipsiDiameter = nan(numel(eFiles),3);
contraDiameter = nan(numel(eFiles),3);
diameterAsymmetry = nan(numel(eFiles),3);
ipsiVolume = nan(numel(eFiles),3);
contraVolume = nan(numel(eFiles),3);
volumeAsymmetry = nan(numel(eFiles),3);

for k = 1:numel(eFiles)
    summaryFile = fullfile(eFiles(k).folder,eFiles(k).name);
    idToken = regexp(summaryFile,'[\\/]sub-([^\\/]+)[\\/]','tokens','once');
    if isempty(idToken); continue; end

    currentID = upper(strtrim(string(idToken{1})));
    currentID = regexprep(currentID,'^SUB-','');
    eicabID(k) = currentID;

    patientRow = find(imagingID==currentID,1);
    if isempty(patientRow); continue; end
    currentSide = side227(patientRow);

    raw = readcell(summaryFile);
    textRaw = strings(size(raw));
    for rr = 1:size(raw,1)
        for cc = 1:size(raw,2)
            try
                textRaw(rr,cc) = strtrim(string(raw{rr,cc}));
            catch
                textRaw(rr,cc) = "";
            end
        end
    end

    [volumeRows,~] = find(strcmpi(textRaw,'Volume'));
    if isempty(volumeRows); continue; end
    volumeRow = volumeRows(1);

    [diameterRows,~] = find(strcmpi(textRaw,'Diameter mean'));
    if isempty(diameterRows); continue; end
    diameterRow = diameterRows(1);

    for t = 1:3
        if currentSide=="L"
            ipsiLabel = leftLabels{t};
            contraLabel = rightLabels{t};
        elseif currentSide=="R"
            ipsiLabel = rightLabels{t};
            contraLabel = leftLabels{t};
        else
            continue
        end

        [~,ipsiColumns] = find(strcmpi(textRaw,ipsiLabel));
        [~,contraColumns] = find(strcmpi(textRaw,contraLabel));
        ipsiValue = NaN; contraValue = NaN;

        if ~isempty(ipsiColumns)
            ipsiCell = raw{volumeRow,ipsiColumns(1)};
            if isnumeric(ipsiCell)
                ipsiValue = double(ipsiCell);
            else
                ipsiValue = str2double(strrep(string(ipsiCell),',',''));
            end
        end
        if ~isempty(contraColumns)
            contraCell = raw{volumeRow,contraColumns(1)};
            if isnumeric(contraCell)
                contraValue = double(contraCell);
            else
                contraValue = str2double(strrep(string(contraCell),',',''));
            end
        end

        % Preserve raw ipsilateral volume exactly as in the earlier script.
        if isscalar(ipsiValue) && isfinite(ipsiValue) && ipsiValue>0
            ipsiVolume(k,t) = ipsiValue;
        end
        if isscalar(contraValue) && isfinite(contraValue) && contraValue>0
            contraVolume(k,t) = contraValue;
        end
        % Asymmetry legitimately requires both sides.
        if isfinite(ipsiVolume(k,t)) && isfinite(contraVolume(k,t))
            volumeAsymmetry(k,t) = (ipsiValue-contraValue) / ...
                (0.5*(abs(ipsiValue)+abs(contraValue)));
        end

        ipsiDiameterValue = NaN;
        contraDiameterValue = NaN;
        if ~isempty(ipsiColumns)
            ipsiCell = raw{diameterRow,ipsiColumns(1)};
            if isnumeric(ipsiCell)
                ipsiDiameterValue = double(ipsiCell);
            else
                ipsiDiameterValue = str2double(strrep(string(ipsiCell),',',''));
            end
        end
        if ~isempty(contraColumns)
            contraCell = raw{diameterRow,contraColumns(1)};
            if isnumeric(contraCell)
                contraDiameterValue = double(contraCell);
            else
                contraDiameterValue = str2double(strrep(string(contraCell),',',''));
            end
        end

        if isscalar(ipsiDiameterValue) && isfinite(ipsiDiameterValue) && ...
                ipsiDiameterValue>0
            ipsiDiameter(k,t) = ipsiDiameterValue;
        end
        if isscalar(contraDiameterValue) && isfinite(contraDiameterValue) && ...
                contraDiameterValue>0
            contraDiameter(k,t) = contraDiameterValue;
        end
        if isfinite(ipsiDiameter(k,t)) && isfinite(contraDiameter(k,t))
            diameterAsymmetry(k,t) = ...
                (ipsiDiameterValue-contraDiameterValue) / ...
                (0.5*(abs(ipsiDiameterValue)+abs(contraDiameterValue)));
        end
    end
end

%% 4. Match all sources by ID
[matchedID,imagingIndex,eicabIndex] = intersect(imagingID,eicabID,'stable');
matchedSide = side227(imagingIndex);
matchedDuration = duration227(imagingIndex);
matchedIpsiDiameter = ipsiDiameter(eicabIndex,:);
matchedContraDiameter = contraDiameter(eicabIndex,:);
matchedDiameterAsymmetry = diameterAsymmetry(eicabIndex,:);
matchedIpsiVolume = ipsiVolume(eicabIndex,:);
matchedContraVolume = contraVolume(eicabIndex,:);
matchedVolumeAsymmetry = volumeAsymmetry(eicabIndex,:);
matchedRawBurden = matchedIpsiDiameter .* matchedIpsiVolume;
matchedAsymmetryBurden = ...
    matchedDiameterAsymmetry .* matchedVolumeAsymmetry;
matchedThickness = thickness227(imagingIndex,:);
matchedFLAIR = flair227(imagingIndex,:);
matchedMyelin = myelin227(imagingIndex,:);
matchedPET = pet227(imagingIndex,:);

fprintf('\n=============================================================\n');
fprintf('eICAB files found:           %d\n',numel(eFiles));
fprintf('IDs matched to imaging:      %d\n',numel(matchedID));
fprintf('Models to fit:               8 (2 exposures x 4 regions)\n');
fprintf('=============================================================\n');

% Construct equal-weight ALL composites. A composite is retained only if
% all three territorial measurements are available.
Z = (matchedRawBurden-mean(matchedRawBurden,1,'omitnan')) ./ ...
    std(matchedRawBurden,0,1,'omitnan');
allRawBurden = mean(Z,2,'omitnan');
allRawBurden(any(~isfinite(matchedRawBurden),2)) = NaN;

Z = (matchedAsymmetryBurden-mean(matchedAsymmetryBurden,1,'omitnan')) ./ ...
    std(matchedAsymmetryBurden,0,1,'omitnan');
allAsymmetryBurden = mean(Z,2,'omitnan');
allAsymmetryBurden(any(~isfinite(matchedAsymmetryBurden),2)) = NaN;

Z = (matchedThickness-mean(matchedThickness,1,'omitnan')) ./ ...
    std(matchedThickness,0,1,'omitnan');
allThickness = mean(Z,2,'omitnan');
allThickness(any(~isfinite(matchedThickness),2)) = NaN;

Z = (matchedFLAIR-mean(matchedFLAIR,1,'omitnan')) ./ ...
    std(matchedFLAIR,0,1,'omitnan');
allFLAIR = mean(Z,2,'omitnan');
allFLAIR(any(~isfinite(matchedFLAIR),2)) = NaN;

Z = (matchedMyelin-mean(matchedMyelin,1,'omitnan')) ./ ...
    std(matchedMyelin,0,1,'omitnan');
allMyelin = mean(Z,2,'omitnan');
allMyelin(any(~isfinite(matchedMyelin),2)) = NaN;

Z = (matchedPET-mean(matchedPET,1,'omitnan')) ./ ...
    std(matchedPET,0,1,'omitnan');
allPET = mean(Z,2,'omitnan');
allPET(any(~isfinite(matchedPET),2)) = NaN;
clear Z

matchedRawBurden(:,4) = allRawBurden;
matchedAsymmetryBurden(:,4) = allAsymmetryBurden;
matchedThickness(:,4) = allThickness;
matchedFLAIR(:,4) = allFLAIR;
matchedMyelin(:,4) = allMyelin;
matchedPET(:,4) = allPET;

% Auditable data remain in the MATLAB workspace.
MatchedData = table(matchedID,matchedSide,matchedDuration, ...
    matchedIpsiDiameter(:,1),matchedIpsiVolume(:,1), ...
    matchedDiameterAsymmetry(:,1),matchedVolumeAsymmetry(:,1), ...
    matchedRawBurden(:,1),matchedAsymmetryBurden(:,1), ...
    matchedIpsiDiameter(:,2),matchedIpsiVolume(:,2), ...
    matchedDiameterAsymmetry(:,2),matchedVolumeAsymmetry(:,2), ...
    matchedRawBurden(:,2),matchedAsymmetryBurden(:,2), ...
    matchedIpsiDiameter(:,3),matchedIpsiVolume(:,3), ...
    matchedDiameterAsymmetry(:,3),matchedVolumeAsymmetry(:,3), ...
    matchedRawBurden(:,3),matchedAsymmetryBurden(:,3), ...
    'VariableNames',{'ID','Side','Duration', ...
    'ACTP_IpsiDiameter','ACTP_IpsiVolume','ACTP_DiameterAsymmetry', ...
    'ACTP_VolumeAsymmetry','ACTP_RawBurden','ACTP_AsymmetryBurden', ...
    'MCAT_IpsiDiameter','MCAT_IpsiVolume','MCAT_DiameterAsymmetry', ...
    'MCAT_VolumeAsymmetry','MCAT_RawBurden','MCAT_AsymmetryBurden', ...
    'PCAT_IpsiDiameter','PCAT_IpsiVolume','PCAT_DiameterAsymmetry', ...
    'PCAT_VolumeAsymmetry','PCAT_RawBurden','PCAT_AsymmetryBurden'});

%% 5. PROCESS Model 21; complete cases only, no outlier removal
% W moderates X->M; V moderates M->Y.
% M = a1*X + a2*W + a3*XW
% Y = cprime*X + b1*M + b2*V + b3*MV

moderatorLevels = [-1 0 1];

SampleSummary = table();
PathResults = table();
ConditionalResults = table();
IndexResults = table();
indexEffect = nan(2,4);
indexSignificant = false(2,4);

for combination = 1:2
    if combination==1
        Xmatrix = matchedRawBurden;
        xName = "Raw mean diameter x raw volume";
    else
        Xmatrix = matchedAsymmetryBurden;
        xName = "Diameter asymmetry x volume asymmetry";
    end
    model21Name = string(combinationNames{combination}) + ...
        " + (-Thickness) + (-MidFLAIR)";
    Mmatrix = matchedThickness;
    Vmatrix = matchedFLAIR;
    vName = "-MidFLAIR";

    for t = 1:4
        Xraw = Xmatrix(:,t);
        Wraw = matchedDuration;
        Mraw = -1*Mmatrix(:,t);
        Vraw = -1*Vmatrix(:,t);
        Yraw = matchedPET(:,t);
        if combination==1 && t<=3
            validX = isfinite(Xraw) & Xraw>0;
        else
            validX = isfinite(Xraw);
        end
        keep = validX & isfinite(Wraw) & isfinite(Mraw) & ...
            isfinite(Vraw) & isfinite(Yraw);

        X = zscore(Xraw(keep));
        W = zscore(Wraw(keep));
        M = zscore(Mraw(keep));
        V = zscore(Vraw(keep));
        Y = zscore(Yraw(keep));
        XW = X.*W;
        MV = M.*V;
        n = numel(X);
        SampleSummary = [SampleSummary;table(combination,model21Name, ...
            xName,vName,string(territoryNames{t}),n,'VariableNames', ...
            {'Combination','CombinationName','X','V','Territory','N'})];
        if n<30; continue; end

        mediatorModel = fitlm(table(X,W,XW,M),'M ~ X + W + XW');
        outcomeModel = fitlm(table(X,M,V,MV,Y),'Y ~ X + M + V + MV');
        a1 = mediatorModel.Coefficients{'X','Estimate'};
        a2 = mediatorModel.Coefficients{'W','Estimate'};
        a3 = mediatorModel.Coefficients{'XW','Estimate'};
        b1 = outcomeModel.Coefficients{'M','Estimate'};
        b2 = outcomeModel.Coefficients{'V','Estimate'};
        b3 = outcomeModel.Coefficients{'MV','Estimate'};
        cp = outcomeModel.Coefficients{'X','Estimate'};

        pathNames = ["a1: X to M";"a2: W to M";"a3: XW to M"; ...
            "b1: M to Y";"b2: V to Y";"b3: MV to Y";"cprime: X to Y"];
        pathBeta = [a1;a2;a3;b1;b2;b3;cp];
        pathP = [mediatorModel.Coefficients{'X','pValue'}; ...
            mediatorModel.Coefficients{'W','pValue'}; ...
            mediatorModel.Coefficients{'XW','pValue'}; ...
            outcomeModel.Coefficients{'M','pValue'}; ...
            outcomeModel.Coefficients{'V','pValue'}; ...
            outcomeModel.Coefficients{'MV','pValue'}; ...
            outcomeModel.Coefficients{'X','pValue'}];
        PathResults = [PathResults;table(repmat(combination,7,1), ...
            repmat(model21Name,7,1),repmat(xName,7,1), ...
            repmat(vName,7,1),repmat(string(territoryNames{t}),7,1), ...
            repmat(n,7,1),pathNames,pathBeta,pathP,pathP<0.05, ...
            'VariableNames',{'Combination','CombinationName','X','V','Territory', ...
            'N','PathName','Beta','PValue','Significant'})];

        observedConditional = nan(9,1);
        conditionW = nan(9,1);
        conditionV = nan(9,1);
        q = 0;
        for wi = 1:3
            for vi = 1:3
                q = q+1;
                conditionW(q) = moderatorLevels(wi);
                conditionV(q) = moderatorLevels(vi);
                observedConditional(q) = (a1+a3*conditionW(q))* ...
                    (b1+b3*conditionV(q));
            end
        end
        observedIndexMM = a3*b3;
        observedIndexW = a3*(b1+b3*moderatorLevels(:));
        observedIndexV = b3*(a1+a3*moderatorLevels(:));

        bootConditional = nan(nBoot,9);
        bootIndexMM = nan(nBoot,1);
        bootIndexW = nan(nBoot,3);
        bootIndexV = nan(nBoot,3);
        A_M = [ones(n,1),X,W,XW];
        A_Y = [ones(n,1),X,M,V,MV];
        rng(20260914,'twister');
        for r = 1:nBoot
            ii = randi(n,n,1);
            bm = A_M(ii,:)\M(ii);
            by = A_Y(ii,:)\Y(ii);
            ba1 = bm(2); ba3 = bm(4);
            bb1 = by(3); bb3 = by(5);
            q = 0;
            for wi = 1:3
                for vi = 1:3
                    q = q+1;
                    bootConditional(r,q) = ...
                        (ba1+ba3*moderatorLevels(wi))* ...
                        (bb1+bb3*moderatorLevels(vi));
                end
            end
            bootIndexMM(r) = ba3*bb3;
            bootIndexW(r,:) = ba3*(bb1+bb3*moderatorLevels);
            bootIndexV(r,:) = bb3*(ba1+ba3*moderatorLevels);
        end

        conditionalCI = prctile(bootConditional,[2.5 97.5],1)';
        conditionalSig = conditionalCI(:,1)>0 | conditionalCI(:,2)<0;
        ConditionalResults = [ConditionalResults;table(repmat(combination,9,1), ...
            repmat(model21Name,9,1),repmat(xName,9,1), ...
            repmat(vName,9,1),repmat(string(territoryNames{t}),9,1), ...
            repmat(n,9,1),conditionW,conditionV,observedConditional, ...
            conditionalCI(:,1),conditionalCI(:,2),conditionalSig, ...
            'VariableNames',{'Combination','CombinationName','X','V','Territory','N', ...
            'W_z','V_z','IndirectEffect','BootCI_L','BootCI_U','Significant'})];

        combinedObservedIndex = [observedIndexMM;observedIndexW;observedIndexV];
        combinedBootIndex = [bootIndexMM,bootIndexW,bootIndexV];
        indexCI = prctile(combinedBootIndex,[2.5 97.5],1)';
        indexSig = indexCI(:,1)>0 | indexCI(:,2)<0;
        indexNames = ["Index moderated moderated mediation"; ...
            "Index W at V=-1SD";"Index W at V=mean";"Index W at V=+1SD"; ...
            "Index V at W=-1SD";"Index V at W=mean";"Index V at W=+1SD"];
        IndexResults = [IndexResults;table(repmat(combination,7,1), ...
            repmat(model21Name,7,1),repmat(xName,7,1), ...
            repmat(vName,7,1),repmat(string(territoryNames{t}),7,1), ...
            repmat(n,7,1),indexNames,combinedObservedIndex,indexCI(:,1), ...
            indexCI(:,2),indexSig,'VariableNames',{'Combination','CombinationName', ...
            'X','V','Territory','N','IndexName','Effect','BootCI_L','BootCI_U','Significant'})];
        indexEffect(combination,t) = observedIndexMM;
        indexSignificant(combination,t) = indexSig(1);
    end
end

%% 6. Display results for both vascular-burden definitions
fprintf('\n================ MODEL 21 SAMPLE SUMMARY ================\n');
disp(SampleSummary);
fprintf('\n================ SIGNIFICANT PATHS (P < 0.05) ================\n');
disp(PathResults(PathResults.Significant,:));
fprintf('\n================ MODERATED-MODERATED MEDIATION INDEX ================\n');
disp(IndexResults(IndexResults.IndexName== ...
    "Index moderated moderated mediation",:));

figure('Color','w','Position',[30 40 1900 920]);
tiledlayout(2,4,'TileSpacing','compact','Padding','compact');

for combination = 1:2
    for t = 1:4
        ax = nexttile;
        hold(ax,'on'); axis(ax,[0 10 0 6]); axis(ax,'off');

        pr = PathResults( ...
            PathResults.Combination==combination & ...
            PathResults.Territory==string(territoryNames{t}),:);
        ir = IndexResults( ...
            IndexResults.Combination==combination & ...
            IndexResults.Territory==string(territoryNames{t}) & ...
            IndexResults.IndexName=="Index moderated moderated mediation",:);

        if isempty(pr) || isempty(ir)
            text(5,3,'Insufficient complete cases', ...
                'HorizontalAlignment','center');
            title(sprintf('%s: %s',territoryNames{t}, ...
                combinationNames{combination}),'FontWeight','bold');
            continue
        end

        beta = pr.Beta;
        pValue = pr.PValue;
        pathColor = repmat({[0.15 0.15 0.15]},7,1);
        for k = 1:7
            if pValue(k)<0.05
                if beta(k)>0
                    pathColor{k}=[0.80 0.10 0.10];
                else
                    pathColor{k}=[0.10 0.35 0.80];
                end
            end
        end

        rectangle('Position',[0.4 1.7 1.8 0.7], ...
            'FaceColor',[0.90 0.95 1.00]);
        rectangle('Position',[4.1 3.2 1.8 0.7], ...
            'FaceColor',[0.92 1.00 0.92]);
        rectangle('Position',[7.8 1.7 1.8 0.7], ...
            'FaceColor',[1.00 0.92 0.92]);
        rectangle('Position',[0.6 4.7 1.5 0.7], ...
            'FaceColor',[1.00 0.97 0.86]);
        rectangle('Position',[7.9 4.7 1.5 0.7], ...
            'FaceColor',[0.96 0.92 1.00]);

        if combination==1
            xBoxText = {'X','Raw D x V'};
        else
            xBoxText = {'X','Asym D x V'};
        end
        text(1.3,2.05,xBoxText,'HorizontalAlignment','center', ...
            'FontWeight','bold');
        text(5.0,3.55,{'M','Thickness'}, ...
            'HorizontalAlignment','center','FontWeight','bold');
        text(8.7,2.05,{'Y','PET'}, ...
            'HorizontalAlignment','center','FontWeight','bold');
        text(1.35,5.05,{'W','Duration'}, ...
            'HorizontalAlignment','center','FontWeight','bold');
        text(8.65,5.05,{'V','MidFLAIR'}, ...
            'HorizontalAlignment','center','FontWeight','bold');

        quiver(2.2,2.10,1.9,1.25,0,'k', ...
            'LineWidth',1.3,'MaxHeadSize',0.25);
        quiver(5.9,3.35,1.9,-1.25,0,'k', ...
            'LineWidth',1.3,'MaxHeadSize',0.25);
        quiver(2.2,1.92,5.6,0,0,'k', ...
            'LineWidth',1.3,'MaxHeadSize',0.12);
        quiver(1.35,4.70,1.85,-1.75,0, ...
            'Color',[0.45 0.45 0.45], ...
            'LineStyle','--','LineWidth',1.1,'MaxHeadSize',0.20);
        quiver(8.65,4.70,-1.85,-1.75,0, ...
            'Color',[0.45 0.45 0.45], ...
            'LineStyle','--','LineWidth',1.1,'MaxHeadSize',0.20);

        text(2.55,3.20,sprintf('a_1 = %.3f, P = %.3f', ...
            beta(1),pValue(1)),'Color',pathColor{1},'FontSize',8);
        text(1.55,4.30,sprintf('a_2 = %.3f, P = %.3f', ...
            beta(2),pValue(2)),'Color',pathColor{2},'FontSize',8);
        text(2.05,3.85,sprintf('a_3(XW) = %.3f, P = %.3f', ...
            beta(3),pValue(3)),'Color',pathColor{3},'FontSize',8);
        text(6.00,3.20,sprintf('b_1 = %.3f, P = %.3f', ...
            beta(4),pValue(4)),'Color',pathColor{4},'FontSize',8);
        text(7.15,4.30,sprintf('b_2 = %.3f, P = %.3f', ...
            beta(5),pValue(5)),'Color',pathColor{5},'FontSize',8);
        text(6.10,3.85,sprintf('b_3(MV) = %.3f, P = %.3f', ...
            beta(6),pValue(6)),'Color',pathColor{6},'FontSize',8);
        text(4.05,1.50,sprintf('c'' = %.3f, P = %.3f', ...
            beta(7),pValue(7)),'Color',pathColor{7},'FontSize',8);

        if ir.Significant
            if ir.Effect>0
                indexColor=[0.80 0.10 0.10];
            else
                indexColor=[0.10 0.35 0.80];
            end
        else
            indexColor=[0.15 0.15 0.15];
        end

        text(5.0,0.65,sprintf( ...
            'Index a_3b_3 = %.3f, bootstrap 95%% CI [%.3f, %.3f]', ...
            ir.Effect,ir.BootCI_L,ir.BootCI_U), ...
            'HorizontalAlignment','center', ...
            'Color',indexColor,'FontWeight','bold','FontSize',8);

        title(sprintf('%s | %s | N=%d',territoryNames{t}, ...
            combinationNames{combination},pr.N(1)), ...
            'FontWeight','bold','FontSize',9);
    end
end

sgtitle(['Model 21: diameter-volume burden, -thickness, ' ...
    '-MidFLAIR, and PET'],'FontWeight','bold','FontSize',13);

return

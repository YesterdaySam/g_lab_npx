
function [ripples,envStd] = get_ripples(root,chan,sess,lThresh,hThresh,dThresh)
%% Returns an updated root
% Based off Jadhav lab methods and Buzsaki lab code
% https://www.cell.com/neuron/fulltext/S0896-6273(19)30785-8
% https://github.com/buzsakilab/buzcode/blob/master/analysis/SharpWaveRipples/bz_FindRipples.m
% 
% Filter 150-250 Hz
% Get envelope of signal with Hilbert Transform
% Detect ripple peaks at >3SD of mean on that electrode
% Start/stop as times when envelope exceeds mean near each event
% Optionally discard excessively short/long ripples
% Optionally discard events during runs >4cm/s
%
% Inputs:
% root = root object. Must have root.tssync and root.tsb fields
% chan = LFP channel ID
% sess = session struct from importBhvr
%
% Outputs:
% ripples = Nx4 array of [start peak end normpwr] indices for N ripples after thresholding
% envStd  = standard deviation of the analytic signal envelope
%
% Created 1/17/25 LKW; Grienberger Lab; Brandeis University
%--------------------------------------------------------------------------

arguments
    root            %struct containing neural info
    chan {double}   %LFP channel ID
    sess            %session struct
    lThresh = 3     %First threshold for ripple envelope
    hThresh = 5     %Second threshold for ripple peak power
    dThresh = [15 250] % Min and Max ripple duration
end

% Get Analytic siganl envelope
rawlf = root.lfp(chan,:); 
riplf = bandpass(rawlf, [150, 250], root.fs_lfp);
ripEnv = abs(hilbert(riplf));

% Set thresholds
envStd = std(ripEnv);
lowThresh = lThresh*envStd;
hiThresh = hThresh*envStd;
durThresh = dThresh;
runThresh = 4;

%% Find candidate start/stops
ripThreshed = ripEnv > lowThresh;
stt = find(diff(ripThreshed)>0);
stp = find(diff(ripThreshed)<0);

% Exclude last ripple if it is incomplete
if length(stp) == length(stt)-1
	stt = stt(1:end-1);
end
% Exclude first ripple if it is incomplete
if length(stp)-1 == length(stt)
    stp = stp(2:end);
end
% Correct special case when both first and last ripples are incomplete
if stt(1) > stp(1)
	stp(1) = [];
	stt(end) = [];
end
firstPass = [stt',stp'];
if isempty(firstPass)
	disp('Detection by thresholding failed');
	return
else
	disp(['After detection by thresholding: ' num2str(length(firstPass)) ' events.']);
end

%% Merge ripples that are nearby
iriThresh = 30/1000*root.fs_lfp; %30 msec inter-ripple-interval
secondPass = [];
rips = firstPass(1,:);
for i = 2:size(firstPass,1)
	if firstPass(i,1) - rips(2) < iriThresh
		% Merge
		rips = [rips(1) firstPass(i,2)];
	else
		secondPass = [secondPass ; rips];
		rips = firstPass(i,:);
	end
end
secondPass = [secondPass ; rips];
if isempty(secondPass)
	disp('Ripple merge failed');
	return
else
	disp(['After ripple merge: ' num2str(length(secondPass)) ' events.']);
end

%% Discard ripples with a peak power < highThresholdFactor
thirdPass = [];
pkNormPwr = [];
for i = 1:size(secondPass,1)
	[maxValue,maxIndex] = max(ripEnv([secondPass(i,1):secondPass(i,2)]));
	if maxValue > hiThresh
		thirdPass = [thirdPass ; secondPass(i,:)];
		pkNormPwr = [pkNormPwr ; maxValue];
	end
end
if isempty(thirdPass)
	disp('Peak thresholding failed.');
	return
else
	disp(['After peak thresholding: ' num2str(length(thirdPass)) ' events.']);
end

%% Detect negative peak position for each ripple
pkPos = zeros(size(thirdPass,1),1);
for i=1:size(thirdPass,1)
	[minValue,minIndex] = min(riplf(thirdPass(i,1):thirdPass(i,2)));
	pkPos(i) = minIndex + thirdPass(i,1) - 1;
end

%% Discard ripples that are way too long
ripples = [root.lfp_tsb(thirdPass(:,1))' root.lfp_tsb(pkPos)' ...
           root.lfp_tsb(thirdPass(:,2))' pkNormPwr];
dur = ripples(:,3)-ripples(:,1);
ripples(dur > durThresh(2) / 1000 * root.fs_lfp, :) = NaN;
disp(['After Long duration test: ' num2str(size(ripples,1)) ' events.']);

%% Discard ripples that are too short
ripples(dur < durThresh(1) / 1000 * root.fs_lfp, :) = NaN;
ripples = ripples((all((~isnan(ripples)),2)),:);

disp(['After Short duration test: ' num2str(size(ripples,1)) ' events.']);

%% Discard ripples during runs
hiVel = sess.velshft(ripples(:,2)) > runThresh;
ripples(hiVel,:) = NaN;
ripples = ripples((all((~isnan(ripples)),2)),:);

disp(['After Run Velocity threshold test: ' num2str(size(ripples,1)) ' events.']);

end
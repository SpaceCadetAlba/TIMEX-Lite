% Void Function TIMEX LITE
function timexLite(noiseThresh1, noiseThresh2, peakThresh, zeroMaskSize, peakSearchWindowSize)

% Read Audio Input---------------------------------------------------------
audioFile = uigetfile;
[audio, fs] = audioread(audioFile);
lengthAudio = length(audio);

audioOrig = audio; % Store original audio

% Plot Audio Input---------------------------------------------------------
figure;
subplot(8,1,1)
plot(audioOrig)
xlim([0, lengthAudio]);
ylabel('Amplitude')
xlabel('Sample Number')
title('Audio')

% Noise Gate 01 -----------------------------------------------------------
% This gate is used for note partitioning

noiseThresh = noiseThresh1; %Get gate threshold

% check all audio samples. If value is less than gate threshold then zero
i = 1; 
while i < (lengthAudio + 1)
   if (audio(i, 1)) < noiseThresh && (audio(i, 1)) > (-1*noiseThresh)
        audio(i, 1) = 0;
    end
    i = i + 1;
end
clearvars i noiseThresh;

audioGated = audio; % store gated audio

% Plot Gated Audio --------------------------------------------------------
subplot(8,1,2)
plot(audioGated)
xlim([0, lengthAudio]);
ylabel('Amplitude')
xlabel('Sample Number')
title('Gated Audio')

% restoring audio vector to original audio for further processing
audio = audioOrig;

% Noise Gate 02 -----------------------------------------------------------
%This gate is used for removing noise for pitch tracking

noiseThresh = noiseThresh2; % Get gate threshold

% use original audio
audio = audioOrig;
% check all audio samples. If value is less than gate threshold then zero
i = 1;
while i < (lengthAudio + 1)
   if (audio(i, 1)) < noiseThresh && (audio(i, 1)) > (-1*noiseThresh)
        audio(i, 1) = 0;
    end
    i = i + 1;
end
clearvars i noiseThresh;

% Estimate Pitch of Audio Input--------------------------------------------

% NCF pitch estimation 
%use audio from by gate2 as we are only interested in 'gate:open' frames
windowSize = round(0.052*fs);
overlapSize = round (0.042*fs);
[fo,location] = pitch(audio,fs, 'WindowLength', windowSize, 'OverlapLength', overlapSize, 'Method', 'NCF', 'MedianFilterLength', 3);

% store fo vector
foOrig = fo;

% Find NAN Values (rest) and replace with zero
nanCheck = isnan(fo);
fo(nanCheck) = 0;

% Zero pitch below 20Hz
i = 1;
while i < ((length(fo))+ 1)
    if fo(i,1)< 20
        fo(i, 1) = 0;
    else
        fo(i,1) = fo(i,1);
    end
    i = i+1;
end
clearvars i;

%Plot unprocessed Pitch----------------------------------------------------
subplot(8,1,3)
plot(location,foOrig)
xlim([0, lengthAudio]);
ylabel('Pitch (Hz)')
xlabel('Sample Number')
title('Raw Pitch')

% Vibrato suppression -----------------------------------------------------

% Define Musical Notes from C0 to B8---------------------------------------
noteMatrix = zeros(2, 108);
% Place Note Numbers in noteMatrix
i = -57:50;
noteMatrix(1,:) = i;
clearvars i;
% Calculate Frequencies for each note number
i = 1;
while i < 109
    noteMatrix(2, i) = 440*(2^((noteMatrix(1, i))/12));
    i = i + 1;
end
clearvars i;

% Pitch Snap---------------------------------------------------------------
% Get Length
pitchLength = length(fo);
% Loop
i = 1;
while i < (pitchLength + 1)
    nanCheck = fo(i, 1);
    if nanCheck == 0
        fo(i, 1) = fo(i, 1);
    else
        n = round((-1)*12*(log2(440/(fo(i, 1)))));
        fo(i, 1) = noteMatrix(2, (n+58));
    end
    i = i+1;
end
clearvars i;

%Plot processed Pitch -----------------------------------------------------
subplot(8,1,4)
plot(location,fo)
xlim([0, lengthAudio]);
ylabel('Pitch (Hz)')
xlabel('Sample Number')
title('Snapped Pitch')

% Find Pitch Slope---------------------------------------------------------
foSlope = diff(fo);

% Make location vector
i = 1;
pitchLength1 = length(fo);
locationSlope = zeros ((pitchLength1-1), 1);
while i < pitchLength1 - 2
    x = (location((i+1), 1)) - (location(i, 1));
    x = x/2;
    x = (location(i, 1)) + x;
    locationSlope(i, 1) = x;
    i = i + 1;
end
clearvars i x;

% Smooth pitch slope
foSlope = smoothdata(foSlope, 'sgolay');

% Threshold slope
i = 1;
slopeThreshold = peakThresh;
while i<((length(foSlope))+1)
    if foSlope(i,1) < slopeThreshold
        foSlope(i,1) = 0;
    else
        foSlope(i,1) = foSlope(i,1);
    end
    i = i+1;
end
clearvars i;

% Plot Pitch Slope---------------------------------------------------------
subplot(8,1,5)
plot(locationSlope,foSlope)
xlim([0, lengthAudio]);
ylabel('Pitch Slope')
xlabel('Sample Number')
title('Pitch Slope')

% Find Peaks --------------------------------------------------------------
[onsetPeaks, onsetLocation]=findpeaks(foSlope);

onsetPeakOrig = onsetPeaks; %store original peak values

%translate peak locations to locations in audio file
%peak locations correlate to a variable value in location matrix for the
%slope
i = 1;
lengthPeaks = length(onsetLocation);
while i < (lengthPeaks + 1)
    n = onsetLocation(i, 1);
    onsetLocation(i, 1) = locationSlope(n, 1);
    i = i+1;
end
clearvars i n;

onsetLocationOrig = onsetLocation; %store unsorted peak locations

%Plot Unsorted Peaks ------------------------------------------------------
subplot(8,1,6)
scatter(onsetLocationOrig,onsetPeakOrig)
xlim([0, lengthAudio]);
ylabel('Onset Peaks')
xlabel('Sample Number')
title('Peaks(Unsorted)')

% Window Each Note --------------------------------------------------------
audio = audioGated; % Get audio from note partitioning gate
audio = audio'; %being lazy
maskTime = zeroMaskSize; % get length of zero mask
masksize = fs * maskTime; %length of mask in samples
check = zeros(1, masksize); % Make mask of zeros to identify gaps between notes(regions where gate trigger is 'off)

b = find(~audio); %return vector, indices of zero value samples
z = lengthAudio;
i = 1;
count = length(b);

% for each zo check if the length of zero block is size of mask block
while i < count && (b(1, i) + (masksize-1)) < z
    x = b(1, i);
    window = audio(x:(x + masksize - 1));
    truewindow = isequal(window, check);
    if truewindow ~= 1
        b(1, i) = 0;
    end
    i = i+ 1;
end

% make vector of value 1 where signal is on and zero where signal is off
c = b(find(b~=0));
d = ones(1, z);
i = 1;
count = length(c);
while i < count
    x = c(1, i);
    d(1, c) = 0;
    i = i+1;
end
gateValue = d; % Store Gate switch

% Find start and end of 'On' blocks
startSamples = strfind([0, gateValue==1],[0 1]);
endSamples = strfind([gateValue==1, 0],[1 0]);

clearvars count x i d c z b truewindow %tidy

% Plot note partitioning---------------------------------------------------
subplot(8,1,7)
plot(gateValue)
xlim([0, lengthAudio]);
ylabel('Amplitude')
xlabel('Sample Number')
title('Onset Detection Gate')

%Peak Sorting Loop---------------------------------------------------------
% Find peak which occurs closest to the gate trigger while gate is open

% For each partitioned note
% Look at first 3 peaks in partition
% Ensure these peaks are within the search window (are suitably close to
% beginning of note partition)
% Select the biggest of the peaks which satisfy these conditions
i = 1;
sorted = 0;
sortedValue = 0;
frame = fs*peakSearchWindowSize; % get peak search window size

while i<(length(startSamples)+1)
    x = startSamples(1, i);%Find where gate opens
    z = onsetLocation - x; %find value closest
    z = abs(z);
    [y, index] = min(z);
    n = 0;
    a = onsetPeaks(index, 1);
    if index<(lengthPeaks-1)
        if (onsetLocation(index+1, 1) - (onsetLocation(index, 1)))<frame
            a = [a onsetPeaks(index+1, 1)];
            n = n+1;
        end
        if (onsetLocation(index+2, 1) - (onsetLocation(index, 1)))<frame
            a = [a onsetPeaks(index+2, 1)];
            n = n+1;
        end
    elseif index == (lengthPeaks-1)
        if (onsetLocation(index+1, 1) - (onsetLocation(index, 1)))<frame
            a = [a onsetPeaks(index+1, 1)];
            n = n+1;
        end
    elseif index == lengthPeaks
        a = a;
    end
    [y2, index2] = max(a);
    indexOrig = length(a)-n;
    indexShift = index2-indexOrig;
    index = index + indexShift;
    sorted = [sorted (onsetLocation(index, 1))]; %Store Values
    sortedValue = [sortedValue (onsetPeaks(index, 1))];
    i= i+1;
end

% Store sorted peaks
sortedIndex = sorted(sorted~=0);
sortedPeaks = sortedValue(sortedValue~=0);

clearvars i a x z y zPos index sorted sortedValue; %Tidy

%Plot Sorted Peaks --------------------------------------------------------
subplot(8,1,8)
scatter(sortedIndex,sortedPeaks)
xlim([0, lengthAudio]);
ylabel('Onset Peaks')
xlabel('Sample Number')
title('Peaks(Sorted)')

%Convert to time ----------------------------------------------------------
sortedIndex = sortedIndex/fs;

%Write to text file--------------------------------------------------------
dlmwrite('Onsets.txt', sortedIndex);



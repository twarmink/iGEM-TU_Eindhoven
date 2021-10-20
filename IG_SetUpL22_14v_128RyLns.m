% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use
%
% File name: SetUpL22_14v128RyLns.m.m - Example of scanline imaging with focused transmits at a single focal depth.
%
% Description:
%   Sequence programming file for L22-14v Linear array, using 128 ray lines
%   (focus transmits) and 128 receive acquisitions. Of the 128 transmit
%   channels, the active transmit aperture is limited based on user-entered
%   transmit focus and f-number. All 128 receive channels are active for
%   each acquisition. This script uses 4X sampling with A/D sample rate of
%   62.5 MHz for a 15.625 MHz processing center frequency.  Transmit is at
%   17.8 MHz and receive bandpass filter has been shifted to 18 MHz center
%   frequency, 13.9MHz -3 dB bandwidth to support the 12 MHz bandwidth of
%   the L22-14v (18MHz center frequency, 67% bandwidth). Processing is
%   asynchronous with respect to acquisition.
%
% Last update:
% 12/07/2015 - modified for SW 3.0

% imagesc(ImgData{1}(:,:,1,1)) 
% imgDataPre = ImgData{1}(:,:,1,1);
% imgDataPost = ImgData{1}(:,:,1,1);
% imagesc(imgDataPost-imgDataPre) 

clear all
PGain = 10.0;
P.startDepth = 0;
P.endDepth = 110;            % Acquisition depth in wavelengths
P.txFocus = 64.3;  % Initial transmit focus.
P.numRays = 128;              % no. of Rays

% Specify system parameters.
Resource.Parameters.connector = 1;
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

% Specify Trans structure array.
Trans.name = 'L22-14v';
Trans.frequency = 19;
Trans.units = 'wavelengths';

Trans = computeTrans(Trans);

RcvProfile.LnaZinSel = 31;

% Specify PData structure array.
PData.PDelta = [0.2*Trans.spacing, 0, 0.2];
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3));
PData.Size(2) = ceil((Trans.numelements*Trans.spacing)/PData.PDelta(1));
PData.Size(3) = 1;      % single image page
PData.Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr
% - specify 128 Region structures.
PData(1).Region = repmat(struct('Shape',struct( ...
                    'Name','Rectangle',...
                    'Position',[0,0,P.startDepth],...
                    'width',Trans.spacing,...
                    'height',P.endDepth-P.startDepth)),1,128);
% - set position of regions to correspond to beam spacing.
for i = 1:128
    PData(1).Region(i).Shape.Position(1) = (-63.5 + (i-1))*Trans.spacing;
end
PData(1).Region = computeRegions(PData(1));

% Specify Media object.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
Resource.RcvBuffer.datatype = 'int16';
Resource.RcvBuffer.rowsPerFrame = 2*184320; % this should be larger than 128*Receive.endDepth*4 for max depth (doubled for 4X sampling)
Resource.RcvBuffer.colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer.numFrames = 10;
Resource.InterBuffer(1).numFrames = 2;  % one intermediate buffer needed.
Resource.ImageBuffer.numFrames = 50;
Resource.DisplayWindow.Title = 'L22-14v128RyLns 4X sampling at 62.5 MHz';
Resource.DisplayWindow.pdelta = 0.25;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 250;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow.Colormap = gray(256);

% Specify TW structure array.
TW = struct('type','parametric',...
            'Parameters',[19,.67,2,1]);

% Specify P.numRays TX structure arrays. Transmit centered on element n in the array for event n.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, P.numRays);

% Determine TX aperture based on focal point and desired f number.
txFNum = 3;  % set to desired f-number value for transmit (range: 1.0 - 20)
P.numTx = round((P.txFocus/txFNum)/Trans.spacing); % no. of elements in 1/2 aperture.
txNumEl = floor(P.numTx/2);
if txNumEl > (Trans.numelements/2 - 1), txNumEl = floor(Trans.numelements/2 - 1); end
% txNumEl is the number of elements to include on each side of the
% center element, for the specified focus and sensitivity cutoff.
% Thus the full transmit aperture will be 2*txNumEl + 1 elements.
%display('Number of elements in transmit aperture:');
%disp(2*txNumEl+1);

% - Set event specific TX attributes.
for n = 1:128   % 128 transmit events
    % Set transmit Origins to positions of elements.
    TX(n).Origin = [(-63.5 + (n-1))*Trans.spacing, 0.0, 0.0];
    % Set transmit Apodization so (1 + 2*TXnumel) transmitters are active.
    lft = n - txNumEl;
    if lft < 1, lft = 1; end
    rt = n + txNumEl;
    if rt > Trans.numelements, rt = Trans.numelements; end
    TX(n).Apod(lft:rt) = 1.0;
    TX(n).Delay = computeTXDelays(TX(n));
end

% Specify Receive structure arrays.
% - We need P.numRays Receives for every frame.

% sampling center frequency is 15.625, but we want the bandpass filter
% centered on the actual transducer center frequency of 18 MHz with 67%
% bandwidth, or 12 to 24 MHz.  Coefficients below were set using
% "G3_BPFdevelopment" with normalized cf=1.15 (18 MHz), bw=0.85,
% xsn wdth=0.41 resulting in -3 dB 0.71 to 1.6 (11.1 to 25 MHz), and
% -20 dB 0.57 to 1.74 (8.9 to 27.2 MHz)
%
BPF1 = [ -0.00009 -0.00128 +0.00104 +0.00085 +0.00159 +0.00244 -0.00955 ...
         +0.00079 -0.00476 +0.01108 +0.02103 -0.01892 +0.00281 -0.05206 ...
         +0.01358 +0.06165 +0.00735 +0.09698 -0.27612 -0.10144 +0.48608 ];

maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
wlsPer128 = 128/(2*4); % wavelengths in 128 samples for 4 samplesPerWave
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', P.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128), ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'BS67BW',...
                        'InputFilter', BPF1, ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, P.numRays*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = P.numRays*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:P.numRays
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [1023,1023,1023,1023,1023,1023,1023,1023]; %[300,511,716,920,1023,1023,1023,1023]; %[0,138,260,287,385,593,674,810];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

TPC.hv = 2.0;

% Specify Recon structure array.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [0,0], ...
               'ImgBufDest', [1,-1], ...
               'RINums', 1:P.numRays);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'replaceIntensity', ...  % replace intensity data
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 0), 1, P.numRays);
% - Set specific ReconInfo attributes.
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j;
end

% Specify Process structure array.
pers = 90;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',PGain,...            % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'timeToNextAcq';
SeqControl(1).argument = 200;  % 200 usec between ray lines
SeqControl(2).command = 'timeToNextAcq';
SeqControl(2).argument = round(44400 - 127*SeqControl(1).argument); % 22.5 frames per second
SeqControl(3).command = 'returnToMatlab';
SeqControl(4).command = 'jump'; % Jump back to start.
SeqControl(4).argument = 1;
nsc = 5;

% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:P.numRays                 % Acquire all ray lines for frame
        Event(n).info = 'Acquire ray line';
        Event(n).tx = j;
        Event(n).rcv = P.numRays*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 1;
        n = n+1;
    end
    % Replace last events SeqControl with inter-frame timeToNextAcq and transfer to host.
    Event(n-1).seqControl = [2,nsc];
       SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
       nsc = nsc+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    if floor(i/4) == i/4     % Exit to Matlab every 4th frame reconstructed
        Event(n).seqControl = 3;
    else
        Event(n).seqControl = 0;
    end
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 4;

% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoffCallback');

% - Range Change
wls2mm = 1;
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        wls2mm = Resource.Parameters.speedOfSound/1000/Trans.frequency;
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',[64,300,P.endDepth]*wls2mm,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');

% - Transmit focus change
UI(3).Control = {'UserB4','Style','VsSlider','Label',['TX Focus (',AxesUnit,')'],...
                 'SliderMinMaxVal',[20,320,P.txFocus]*wls2mm,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(3).Callback = text2cell('%TxFocusCallback');

% - F number change
UI(4).Control = {'UserB3','Style','VsSlider','Label','F Number',...
                 'SliderMinMaxVal',[0,10,round(P.txFocus/(P.numTx*Trans.spacing))],'SliderStep',[0.05,0.1],'ValueFormat','%2.0f'};
UI(4).Callback = text2cell('%FNumCallback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 4;

% Save all the structures to a .mat file.
filename = 'iGEM/MatFiles/L22-14v_128RyLns';
save(filename);
% filename = ('L22-14v_128RyLns'); % 
VSX;    % permits immediately running VSX without specifying the matfile name
return

% **** Callback routines to be converted by text2cell function. ****
%SensCutoffCallback - Sensitivity cutoff change
ReconL = evalin('base', 'Recon');
for i = 1:size(ReconL,2)
    ReconL(i).senscutoff = UIValue;
end
assignin('base','Recon',ReconL);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Recon'};
assignin('base','Control', Control);
return
%SensCutoffCallback

%RangeChangeCallback - Range change
simMode = evalin('base','Resource.Parameters.simulateMode');
% No range change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.endDepth'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
P.endDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P.endDepth = UIValue*scaleToWvl;
    end
end
assignin('base','P',P);

PData = evalin('base','PData');
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
PData(1).Region = repmat(struct('Shape',struct( ...
                    'Name','Rectangle',...
                    'Position',[0,0,P.startDepth],...
                    'width',Trans.spacing,...
                    'height',P.endDepth-P.startDepth)),1,128);
% - set position of regions to correspond to beam spacing.
for i = 1:128
    PData(1).Region(i).Shape.Position(1) = (-63.5 + (i-1))*Trans.spacing;
end
assignin('base','PData',PData);
evalin('base','PData(1).Region = computeRegions(PData(1));');
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback

%TxFocusCallback - TX focus changel
simMode = evalin('base','Resource.Parameters.simulateMode');
% No focus change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.txFocus'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
P.txFocus = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P.txFocus = UIValue*scaleToWvl;
    end
end
assignin('base','P',P);

TX = evalin('base', 'TX');
for n = 1:128   % 128 transmit events
    TX(n).focus = P.txFocus;
    TX(n).Delay = computeTXDelays(TX(n));
end
assignin('base','TX', TX);

% Update Fnumber based on new P.txFocus
evalin('base','set(UI(4).handle(2),''Value'',round(P.txFocus/(P.numTx*Trans.spacing)));');
evalin('base','set(UI(4).handle(3),''String'',num2str(round(P.txFocus/(P.numTx*Trans.spacing))));');
% Set Control command to update TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
return
%TxFocusCallback

%FNumCallback - F number change
simMode = evalin('base','Resource.Parameters.simulateMode');
P = evalin('base','P');
Trans = evalin('base','Trans');
% No F number change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',round(P.txFocus/(P.numTx*Trans.spacing)));
    return
end
P.txFNum = UIValue;
P.numTx = round(P.txFocus/(P.txFNum*Trans.spacing));
assignin('base','P',P);
% - Redefine event specific TX attributes for the new P.numTx.
TX = evalin('base', 'TX');
for n = 1:128   % 128 transmit events
    % Set transmit Apodization.
    lft = n - floor(P.numTx/2);
    if lft < 1, lft = 1; end;
    rt = n + floor(P.numTx/2);
    if rt > Trans.numelements, rt = Trans.numelements; end;
    TX(n).Apod = zeros(1,Trans.numelements);
    TX(n).Apod(lft:rt) = 1.0;
    TX(n).Delay = computeTXDelays(TX(n));
end
assignin('base','TX', TX);
% Set Control command to update TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
return
%FNumCallback

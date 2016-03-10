
downSampleFreq=100000;
readParams.startTime = 1;
readParams.endTime = 9.7;
% readParams.endTime = 4.7;
readParams.highCut = 5000;
readParams.roundTimesVal = 1;
readParams.numThreshSweeps = 10;
readParams.doDownSample = 1;

readParams.doPlotSweeps = true;
readParams.doMarkEvents = false;
readParams.doPlotEvents = false;

detectParams.timesRMS = 2;

... minimum allowable amplitude for alpha functions (in units of samples)
    detectParams.minAmp = -1000;
... maximum allowable amplitude
    detectParams.maxAmp = -5;
... minimum allowable tau for alpha functions (in units of samples)
    detectParams.minTau = 0;%original 500e-6,50e-6
... maximum allowable tau
    detectParams.maxTau = inf;%original 5e-3,worked 50e-3
... minimum allowable yOffset for alpha functions (in units of mV)
    detectParams.minYOffset = -2000;
... maximum allowable yOffset
    detectParams.maxYOffset = 2000;
... minimum allowable decay tau
    detectParams.minDecay = 0; % CHANGE FROM 1E-3 ON 11-05-19 (Scott had set to 3e-3)
... maximum allowable decay tau
    detectParams.maxDecay = inf; % CHANGE FROM 50E-3 ON 11-05-19  (Scott had set to 20e-3)
... threshold used to determine if the change of derivative is great
    ... enough to separate out two alpha functions
    detectParams.derThresh = 5;%original 10,5 worked for cell4, derthresh=5 for 3 and 4
... second EPSP is not recognized if it is closer than this value to
    ... the previous one (in units of samples)
    detectParams.closestEPSPs = 1E-3;
% threshold for standard error above which a fit with multiple alphas is
% attempted
detectParams.errThresh = 0.08;
% 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
detectParams.dataFilterType = 3;
% 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
detectParams.derFilterType = 3;
% length of data filter
detectParams.dataFilterLength = 5E-3;
% length of derivative filter
detectParams.derFilterLength = 2E-3;
% if set to 1 then debugging figures appear
detectParams.debugging = 0;
% index of first data point
detectParams.dataStart = 1;
% forces a graphical output even if other outputs are taken
detectParams.forceDisplay = 0;
% turns off best fitting of start time and decay Tau when 1
detectParams.noFit = 0;

detectParams.threshVal=-1;




filepath='/Users/ero/Documents/Axon/';
savefile='/Users/ero/Tsien_Lab/Analysis/mEPSCs_Matlab/LPI_cell';
numcell=1;


SF=10000;%sampling frequency
for i=1:numcell
    
    if i==1
        filepath2=strcat(filepath,'16-1-24/16124');
        range=22:46;
        exclude=[31;32;43;45];
        savefile2=strcat(savefile,num2str(i),'.mat');
        checkthresh=1;
        PSPsDown=1;
        timerange=1.15e4;
    elseif i==2
        
    elseif i==3
        
    elseif i==4
        
    elseif i==5
        
        
    elseif i==7
        
    end
    results={};
    pp=1;%dummy counter to iterate through results
    for j=1:length(range)
        if ~ismember(range(j),exclude)
            if range(j)<10
                filepath3=strcat(filepath2,'00',num2str(range(j)),'.abf');
            elseif range(j)>9&&range(j)<100
                filepath3=strcat(filepath2,'0',num2str(range(j)),'.abf');
            elseif range(j)>99
                filepath3=strcat(filepath2,num2str(range(j)),'.abf');
            end
            x=abfload(filepath3);
            x=squeeze(x(timerange:end,1,:));
            
            for k=1:length(x(1,:))
                [b,a]=butter(6,500/SF,'low');
                filtData=filtfilt(b,a,x(:,k));
                if checkthresh
                    for iii=1:9
                        filtdatacut=filtData(((iii-1)*20000+1):iii*20000);
                        detectPSPs(filtdatacut, PSPsDown, ...
                            'minAmp', detectParams.minAmp, ...
                            ... minimum allowable amplitude for alpha functions (in units of samples)
                            'maxAmp', detectParams.threshVal, ...
                            ... maximum allowable amplitude
                            'minTau', (detectParams.minTau * downSampleFreq), ...
                            ... minimum allowable tau for alpha functions (in units of samples)
                            'maxTau', (detectParams.maxTau * downSampleFreq), ...
                            ... maximum allowable tau
                            'minYOffset', detectParams.minYOffset, ...
                            ... minimum allowable yOffset for alpha functions (in units of mV)
                            'maxYOffset', detectParams.maxYOffset, ...
                            ... maximum allowable yOffset
                            'minDecay', (detectParams.minDecay * downSampleFreq), ...
                            ... minimum allowable decay tau
                            'maxDecay', (detectParams.maxDecay * downSampleFreq), ...
                            ... maximum allowable decay tau
                            'derThresh', (detectParams.derThresh), ...
                            ... threshold used to determine if the change of derivative is ...
                            ... great enough to separate out two alpha functions
                            'closestEPSPs', (detectParams.closestEPSPs * downSampleFreq), ...
                            ... second EPSP is not recognized if it is closer than this ...
                            ... value to the previous one (in units of samples)
                            'errThresh', (detectParams.errThresh), ...
                            ... threshold for standard error above which a fit with ...
                            ... multiple alphas is attempted
                            'dataFilterType', detectParams.dataFilterType, ...
                            ... 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
                            'derFilterType', detectParams.derFilterType, ...
                            ... 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
                            'dataFilterLength', (detectParams.dataFilterLength * downSampleFreq), ...
                            ... length of data filter
                            'derFilterLength', (detectParams.derFilterLength * downSampleFreq), ...
                            ... length of derivative filter
                            'debugging', detectParams.debugging, ...
                            ... if set to 1 then debugging figures appear
                            'dataStart', detectParams.dataStart, ...
                            ... index of first data point
                            'forceDisplay', detectParams.forceDisplay,...
                            ... forces a graphical output even if other outputs are taken
                            'noFit', detectParams.noFit)
                        
                        pause();
                        close;
                    end
                else
                    results{pp}=detectPSPs(filtData, PSPsDown, ...
                        'minAmp', detectParams.minAmp, ...
                        ... minimum allowable amplitude for alpha functions (in units of samples)
                        'maxAmp', detectParams.threshVal, ...
                        ... maximum allowable amplitude
                        'minTau', (detectParams.minTau * downSampleFreq), ...
                        ... minimum allowable tau for alpha functions (in units of samples)
                        'maxTau', (detectParams.maxTau * downSampleFreq), ...
                        ... maximum allowable tau
                        'minYOffset', detectParams.minYOffset, ...
                        ... minimum allowable yOffset for alpha functions (in units of mV)
                        'maxYOffset', detectParams.maxYOffset, ...
                        ... maximum allowable yOffset
                        'minDecay', (detectParams.minDecay * downSampleFreq), ...
                        ... minimum allowable decay tau
                        'maxDecay', (detectParams.maxDecay * downSampleFreq), ...
                        ... maximum allowable decay tau
                        'derThresh', (detectParams.derThresh), ...
                        ... threshold used to determine if the change of derivative is ...
                        ... great enough to separate out two alpha functions
                        'closestEPSPs', (detectParams.closestEPSPs * downSampleFreq), ...
                        ... second EPSP is not recognized if it is closer than this ...
                        ... value to the previous one (in units of samples)
                        'errThresh', (detectParams.errThresh), ...
                        ... threshold for standard error above which a fit with ...
                        ... multiple alphas is attempted
                        'dataFilterType', detectParams.dataFilterType, ...
                        ... 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
                        'derFilterType', detectParams.derFilterType, ...
                        ... 0 = none 1 = windowFilter, 2 = medianFilter, 3 = savitsky-golay
                        'dataFilterLength', (detectParams.dataFilterLength * downSampleFreq), ...
                        ... length of data filter
                        'derFilterLength', (detectParams.derFilterLength * downSampleFreq), ...
                        ... length of derivative filter
                        'debugging', detectParams.debugging, ...
                        ... if set to 1 then debugging figures appear
                        'dataStart', detectParams.dataStart, ...
                        ... index of first data point
                        'forceDisplay', detectParams.forceDisplay,...
                        ... forces a graphical output even if other outputs are taken
                        'noFit', detectParams.noFit);
                    pp=pp+1;
                    
                    
                end
                
            end
            
          
        else
            results{pp}=NaN;
            pp=pp+1;
            
        end
    end
    save(savefile,'results');
    
    
    
    
    
    
    
    
end
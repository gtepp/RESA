
% Detector for Repeating Event/Slow Clap Sequences
%
% --------------------
%
% Work flow:
%
% - STA-LTA event detector
%
% - Check events with templates (add templates if no match)
% 
% - If match, stack event with template
%
% - At end of each time window, merge any highly correlated templates
%
% - When enough matching events have accumulated, declare a sequence
%
% - When enough stations have an active sequence (or drop below), send alert
%
% --------------------
%
% - Can read CSS, SAC, and generic binary (data only) files or get data
%       from a Winston waveserver or irisFetch
%
% Required functions:
%
%   - readcss.m (if using CSS files)
%   - opencss.m (if using CSS files)
%   - extractdatairis (if data will be retrieved via irisFetch)
%
% - Uses the GISMO Correlation Toolbox and Waveform Suite
%
%       By: Michael West and Celso Reyes (respectively)
%           Geophysical Institute, Alaska Volcano Observatory, U. Alaska Fairbanks
%       Currently maintained by Glenn Thompson
%           University of South Florida
%
%       Website: https://geoscience-community-codes.github.io/GISMO/
%
% - Make sure the toolbox is installed. It is freely avaiable at the above website.
%
% - Requires updated version of adjusttrig.m from GISMO suite to properly
%       align traces
%
% - Note: Currently can only use one channel per station name
%
% - If using files, a file database array is needed with the following format:
%        {filename start_time end_time sample_rate start_time_vec}

% --------------------

% By: Gabrielle Tepp, USGS AVO
% Created: 3/15/2017
% Last updated: 12/18/2017

%--------------------------------------------------------------------------%

clear variables global; % clear workspace
warning off % don't display warnings


%% Hard-wired info

% Set start & end of desired data range

params.time.start.str = '15-Apr 2017 18:00:00'; % start time in a string format
params.time.end.str = '16-Apr 2017 06:00:00'; % end time in a string format

% If you want to set SCNL for convenience, else comment out & code will prompt

params.net.str = 'AV'; % currently assumes only one network or array (could easily be changed)
params.sta.str = 'MAPS,MGOD,MSW ,OKFG,OKNC'; % list of stations being used (names must be same length - whitespaces removed later)
params.cha.str = 'BHZ,BHZ,EHZ,BHZ,BHZ'; % need same # as stations in same (pair) order

directory = './'; % directory to save files to

logfile = 'test.txt';

% General data parameters

dbtype = 'WIN'; % database type (Winston WS ('WIN', default),irisFetch ('IRIS'),'files')
ftype = 'SAC'; % file type if using files (options: CSS, BIN (binary), SAC)

mySource = datasource('winston',IPaddress,Port); % data source object for Winston WS

datadir = '/directory/data/'; % data file directory if using files
fileinfo = 'Volc_fileinfo.mat'; % list of file information if using files

twin = 10.8; % in minutes; data window to search
% twin includes the time window overlap (set = l_lta + evwin to avoid duplicate events)

% Filter parameters (bandpass)

poles = 4;

lf = 2; % in Hz
hf = 7; % in Hz

% Turn on waveform-event plotting (for testing)

evtspl = 0; % 1 = on, 0 = off (default)

% STA-LTA Parameters

min_sep = 10; % in seconds
min_dur = 6; % in seconds
l_sta = 1; % in seconds (measured at the right end of the LTA window)
l_lta = 20; % in seconds

threson = 2.3; % threshold for event on
thresoff = 0.6; % threshold for event off

% Event Waveform & X-correlation parameters

mincc = 0.625; % minimum cross-correlation value
maxlag = 3; % in sec; max allowable lag time for match (note: may be good to set as 2*l_sta)
xcst = 1; % in sec; cross-correlation start time
xcend = 16; % in sec; cross-correlation end time (must be <= evwin)

buff = 2; % in sec; time of window start before event onset
evwin = 28; % in sec; length of event waveform window

tempmax = 15; % max # of active template events to hold in memory

keepseqtemp = 1; % 1=yes, 0=no (default); keep templates of all sequences in memory (for post-script analysis)

% Sequence parameters

minev = 5; % minimimum # of events required for sequence to turn on
seqT = 90; % in minutes; time period within which sequence can be declared
minevoff = 3; % minimimum # of events required for sequence to remain on (will turn off when # drops below this)
seqToff = seqT; % in minutes; time period to use for sequence off

l2chk = 1; % check for increasing (level 2) sequence (1 for on, 0 for off (default))
minevl2 = 5; % minimimum # of events required for increased (lv. 2) sequence
seqTl2 = 90; % in minutes; time period within which increased (lv. 2) sequence can be declared
minevoffl2 = 3; % minimimum # of events required for (lv. 2) sequence to remain on (will turn off when # drops below this)
seqToffl2 = seqTl2; % in minutes; time period to use for (lv. 2) sequence off

minsta = 2; % min # of stations with seq on needed to send alert
minreq = 1; % min # of stations required from preferred list to send alert
reqstastr = 'MAPS,MGOD,OKRE,OKER'; % list of preferred/required stations
offwait = 0; % # of extra time windows with no sequences on before level 1 "off" alert is sent (set to 0 for none)
offwaitl2 = 0; % # of extra time windows with no sequences on before level 2 "off" alert is sent (set to 0 for none)
seqoffwin = 0; % # of extra time windows to wait before station-level sequence turns "off" (set to 0 for none)
alertont = 0; % in minutes; amount of time required before alert status can change to "off" (set to 0 for none)
alertontl2 = 0; % in minutes; amount of time required before level 2 alert status can change to "off" (set to 0 for none)

% Make plots

pltl = 1; % timeline of sequences
splrows = 4; % # of subplot rows per figure for seq timeline (each row is <= 1 week of time)
plstatl = 0; % timeline of sequences on each station
plseqs = 0; % individual sequences with events and on/off times

% These files are optional to add more info to plots (in format [datenum_start datenum_end])
% scfile = 'Volc_sequences.mat'; % file with start/end times of known sequences
% expfile = 'Volc_explosions.mat'; % file with start/end times of known explosions/eruptions
% staoutfile = 'Volc_station_outages.mat'; % file with start/end times of known station outages

%sendtoonl2 = {'email@host.com'}; % to test alert notifications

% Write parameters to log file

if ~exist([directory logfile],'file')
    
    fID = fopen([directory logfile],'a'); % open log file for writing (append)
    
    % write parameters
    
    fprintf(fID,strcat('General Parameters\ntwin = %7.2f\nlf = %7.2f\nhf = %7.2f\npoles = %7.2f\ntempmax = %7.2f\n',...
        '\nSTA-LTA Parameters\nmin_sep = %7.2f\nmin_dur = %7.2f\nl_sta = %7.2f\nl_lta = %7.2f',...
        '\nthreson = %7.2f\nthresoff = %7.2f\n\nCross-correlation Parameters\nmincc = %8.3f\nmaxlag = %7.2f\n',...
        'xcst = %7.2f\nxcend = %7.2f\nbuff = %7.2f\nevwin = %7.2f\n\nSequence Parameters\n',...
        'minev = %7.2f\nseqT = %7.2f\nminevoff = %7.2f\nseqToff = %7.2f\nl2chk = %7.2f\nminevl2 = %7.2f\n',...
        'seqTl2 = %7.2f\nminevoffl2 = %7.2f\nseqToffl2 = %7.2f\nminsta = %7.2f\nminreq = %7.2f\nreqstastr = %s\n\n'),...
        twin,lf,hf,poles,tempmax,min_sep,min_dur,l_sta,l_lta,threson,thresoff,mincc,maxlag,xcst,xcend,buff,evwin,minev,...
        seqT,minevoff,seqToff,minevl2,seqTl2,minevoffl2,seqToffl2,minsta,minreq,reqstastr);
    
end

    
%% Get Information

% Get array, station, channel, and network info

disp(' ')

if ~exist('params','var') && isfield(params,'net') && isfield(params.net,'str') % if no channel is set, ask for one
    
    params.net.str = input('Enter the network code(s) (* for all): ','s');
    
end

if ~exist('params','var') && isfield(params,'sta') && isfield(params.sta,'str') % if no station is set, ask for one
    
    params.sta.str = input('Enter a station(s) (* for all): ', 's');
    
end

if ~exist('params','var') && isfield(params,'cha') && isfield(params.cha,'str') % if no channel is set, ask for one
    
    params.cha.str = input('Enter a channel(s) (* for all, same # and pair order as stations): ','s');
    
end

% Put into useable lists

params.net.list = makelist(params.net.str,'str');

params.sta.list = makelist(params.sta.str,'str');

params.cha.list = makelist(params.cha.str,'str');

% Make sure there are station-channel pairs

while size(params.sta.list,1) ~= size(params.cha.list,1)
    
    disp(' ')
    disp('Number of stations and channels do not match. Please try again.')
    disp(' ')
    params.sta.str = input('Enter a station(s) (* for all): ', 's');
    params.cha.str = input('Enter a channel(s) (* for all): ','s');
    
    params.sta.list = makelist(params.sta.str,'str');
    params.cha.list = makelist(params.cha.str,'str');
    
end

% Get times to examine

if ~exist('params','var') && isfield(params,'time') && isfield(params.sta,'start') % if no start time is set, ask for one
    
    params.time.start.str = input('Enter a starting time: ', 's');
    
end

if ~exist('params','var') && isfield(params,'time') && isfield(params.cha,'end') % if no end time is set, ask for one
    
    params.time.end.str = input('Enter an ending time: ','s');
    
end

% Make preferred station list

if ~exist('reqstastr','var') % if no stations given
    
    reqsta = ''; % empty list
    
else

    reqsta = makelist(reqstastr,'str'); % make usable list of preferred stations

end

% Make sure minimum required stations = 0 if none are provided

if isempty(reqsta) == 1
    
   minreq = 0; 
    
end

ovlp = (l_lta + evwin) / 60; % set overlap of time windows; in min

while twin <= ovlp
    
   disp(' ')
   disp('Warning! Time window is shorter than overlap time.')
   disp('Please try again.')
   disp(' ')
   twin = input('How long is the data time window to analyze? (minutes): ');
    
end

% Convert times to days

twin = twin / 60 / 24;
ovlp = ovlp / 60 / 24;
seqT = seqT / 60 / 24;
seqTl2 = seqTl2 / 60 / 24;
seqToff = seqToff / 60 / 24;
seqToffl2 = seqToffl2 / 60 / 24;
alertont = alertont / 60 / 24;
alertontl2 = alertontl2 / 60 / 24;

% Check that parameters are good

if xcend > evwin
    
    xcend = evwin;
    
end

% If using files, read in file database

if ~isempty(dbtype) && strcmpi(dbtype,'files') == 1
    
    finfo = importdata(fileinfo); % read in file info
    
end


%% Start Detector

tstrt = datenum(params.time.start.str); % convert to datenum
tend = datenum(params.time.end.str); % convert to datenum

if tstrt > tend % if start is after end, switch them

    disp(' ')
    disp('Oops! Start time is after end time. Switching....')
    
    tempt = tstrt;
    tempts = params.time.start.str;
   
    tstrt = tend;
    params.time.start.str = params.time.end.str;
    
    tend = tempt;
    params.time.end.str = tempts;
    
    clear tempt tempts; % clear temporary variables
    
end

% Print out headers for sequences

disp(' ')
disp('      1st Event              Most Recent Event      Total # of Events    Template ID    On/Off    Station')
disp('----------------------    ----------------------    -----------------    -----------    ------    -------')
disp(' ')

% Start new log file if it doesn't yet exist

fID = fopen([directory logfile],'a'); % open log file for writing (append)

% write headers

fprintf(fID,'%s\n%s\n\n','      1st Event              Most Recent Event      Total # of Events    Template ID    On/Off    Station',...
    '----------------------    ----------------------    -----------------    -----------    ------    -------');

fclose(fID);

% Initiate variables/fields

tcheck = tstrt;

for ss = 1:size(params.sta.list,1)
    
    sta = strcat(params.sta.list(ss,:));
    
    data.(sta).seqon = [];
    
    if l2chk == 1
        
        data.(sta).seqonl2 = []; % lv 2 on list
        
    end
    
    data.(sta).gapon = 0;
    
    data.(sta).gap = [];
    
end

alertt.lv1.on = [];
alertt.lv1.off = [];

curstat = 0; % start current sequence status at off

if l2chk == 1
    
    alertt.lv2.on = [];
    alertt.lv2.off = [];
    
    curstatl2 = 0; % start current level 2 sequence status at off
    
end

disp(['Starting at ',datestr(tcheck)])

fID = fopen([directory logfile],'a'); % open log file for writing (append)

% Apply detector to each time window

for t = tstrt:(twin - ovlp):(tend - twin)
    
    % Display time update approx every 5 hours
    
    if ((t - tcheck) * 24) > (5) 
        
        tcheck = t;
        
        disp(' ')
        disp(['Now at ' datestr(t)])
        
    end
    
    % If using files, find files that have data from this time window
    
    if ~isempty(dbtype) && strcmpi(dbtype,'files') == 1
    
        tfiles = find((cell2mat(finfo(:,2)) <= t & cell2mat(finfo(:,3)) > t) | (cell2mat(finfo(:,2)) < t+twin & cell2mat(finfo(:,3)) >= t+twin));
    
    end
    
    for s = 1:size(params.sta.list,1)
        
        %% Get data
        
        sta = strcat(params.sta.list(s,:)); % get station name without any whitespaces
        cha = params.cha.list(s,:); % get channel
            
        if ~isempty(dbtype) && strcmpi(dbtype,'IRIS') == 1
            
            % Get data from IRIS
            
            tempw = irisFetch.Traces(params.net.list,sta,'--',cha,t,(t+twin));
            
            if isempty(tempw) == 1 % if no data retrieved
                
                data.(sta).wf = []; % make empty waveform
                
            else % put data into waveform object
                
                tempd = extractdatairis(tempw,tempw(1).sampleRate,t,(t+twin),NaN); % combine data "chunks"
                
                data.(sta).wf = waveform(sta,cha,tempw(1).sampleRate,t,tempd); % put into waveform object
                
            end
            
        elseif ~isempty(dbtype) && strcmpi(dbtype,'files') == 1
            
            data.(sta).wf = []; % make sure waveform starts empty

            % If full time window is not in current waveform (or no current wf exists), get needed data
            
            if isempty(data.(sta).dwf) || (get(data.(sta).dwf,'start') >= t || get(data.(sta).dwf,'end') <= (t+twin))
                
                % Find files for this station-channel for this time window
                
                goodf = []; % clear variable
                
                for onf = 1:size(tfiles,1)
                    
                    if ((strcmpi(ftype,'CSS') == 1 && strcmpi(finfo{tfiles(onf),5},params.net.list(s,:)) == 1) ||...
                            (strcmpi(finfo{tfiles(onf),5},sta) == 1 && strcmpi(ftype,'CSS') == 0)) && strcmpi(finfo{tfiles(onf),6},cha) == 1
                        
                        goodf = [goodf;tfiles(onf)];
                        
                    end
                    
                end
                
                if size(goodf,1) == 1 % if time window in only one file
                    
                    % Get data from file
                    
                    if strcmpi(ftype,'CSS') == 1
                        
                        cssdata = opencss(datadir,'station',sta,'channel',cha,'files',finfo{goodf,1});
                        
                        stt = (cssdata.file1.(sta).(cha).begtime / 24 / 3600) + datenum([1970 1 1 0 0 0]); % in days
                        
                        data.(sta).dwf = waveform(sta,cha,cssdata.file1.(sta).(cha).samprate,stt,cssdata.file1.(sta).(cha).amp);
                        
                    elseif strcmpi(ftype,'BIN') == 1
                        
                        % Find file(s) with correct
                        
                        binid = fopen([datadir,finfo{goodf,1}]);
                        bindata.(sta).(cha).amp = fread(binid,'single');
                        fclose(binid);
                        
                        % if max value is some huge number, there's missing data - replace with NaN
                        
                        if max(abs(bindata.(sta).(cha).amp)) > 1*10^7
                            
                            bindata.(sta).(cha).amp(abs(bindata.(sta).(cha).amp) > 1*10^7) = NaN;
                            
                        end
                        
                        bindata.(sta).(cha).samprate = length(bindata.(sta).(cha).amp) / 3600 / 24; % find the sampling rate (assuming file is 1 day long)
                        
                        data.(sta).dwf = waveform(sta,cha,bindata.(sta).(cha).samprate,finfo{goodf,2},bindata.(sta).(cha).amp);
                        
                    elseif strcmpi(ftype,'SAC') == 1
                        
                        ds = datasource('sac',[datadir,finfo{goodf,1}]);
                        
                        scnl = scnlobject(sta,cha,params.net.list(1,:),'--');
                        
                        data.(sta).dwf = waveform(ds,scnl,finfo{goodf,2},finfo{goodf,3});
                        
                    end
                    
                    % Extract current data window
                    
                    data.(sta).wf = extract(data.(sta).dwf,'TIME',t,(t+twin));
                    
                elseif size(goodf,1) > 1 % overlaps two or more day-files
                    
                    % If window start time is within current file, will only need to get later file(s)
                    
                    if ~isempty(data.(sta).dwf) && get(data.(sta).dwf,'start') <= t && get(data.(sta).dwf,'end') >= t

                        % Clear any earlier "extra" files since no longer needed
                        
                        if isfield(data.(sta),'extradata') == 1
                            
                            data.(sta) = rmfield(data.(sta),'extradata');
                            
                        end
                        
                        % Hold onto the part of the current file that's still needed (now and for future windows)
                        
                        data.(sta).extradata.wf1 = extract(data.(sta).dwf,'TIME',t,get(data.(sta).dwf,'end')); % from window start to file end
                        
                        % Get other needed files
                        
                        for onf = 2:size(goodf,1) % assume first matching file is current one
                            
                            if strcmpi(ftype,'CSS') == 1
                                
                                cssdata = opencss(datadir,'station',sta,'channel',cha,'files',finfo{goodf(onf),1});
                                
                                stt = (cssdata.file1.(sta).(cha).begtime / 24 / 3600) + datenum([1970 1 1 0 0 0]); % in days
                                
                                if onf == size(goodf,1) % last file is "current" file
                                    disp('Opening last needed file')
                                    data.(sta).dwf = waveform(sta,cha,cssdata.file1.(sta).(cha).samprate,stt,cssdata.file1.(sta).(cha).amp);
                                    
                                else
                                    disp(['Opening file ' num2str(onf)])
                                    data.(sta).extradata.(['wf',num2str(onf)]) =...
                                        waveform(sta,cha,cssdata.file1.(sta).(cha).samprate,stt,cssdata.file1.(sta).(cha).amp);
                                    
                                end
                                
                            elseif strcmpi(ftype,'BIN') == 1
                                
                                binid = fopen([datadir,finfo{goodf(onf),1}]);
                                bindata.(sta).(cha).amp = fread(binid,'single');
                                fclose(binid);
                                
                                % if min value is some huge number, there's missing data - replace with NaN
                                
                                if max(abs(bindata.(sta).(cha).amp)) > 1*10^7
                                    
                                    bindata.(sta).(cha).amp(abs(bindata.(sta).(cha).amp) > 1*10^7) = NaN;
                                    
                                end
                                
                                bindata.(sta).(cha).samprate = length(bindata.(sta).(cha).amp) / 3600 / 24; % find the sampling rate (assuming file is 1 day long)
                                
                                if onf == size(goodf,1) % last file is "current" file
                                    
                                    data.(sta).dwf = waveform(sta,cha,bindata.(sta).(cha).samprate,finfo{goodf(onf),2},bindata.(sta).(cha).amp);
                                    
                                else
                                    
                                    data.(sta).extradata.(['wf',num2str(onf)]) =...
                                        waveform(sta,cha,bindata.(sta).(cha).samprate,finfo{goodf(onf),2},bindata.(sta).(cha).amp);
                                    
                                end
                                
                            elseif strcmpi(ftype,'SAC') == 1
                                
                                ds = datasource('sac',[datadir,finfo{goodf(onf),1}]);
                                
                                scnl = scnlobject(sta,cha,params.net.list(1,:),'--');
                                
                                if onf == size(goodf,1) % last file is "current" file
                                    
                                    data.(sta).dwf = waveform(ds,scnl,finfo{goodf(onf),2},finfo{goodf(onf),3});
                                    
                                else
                                    
                                    data.(sta).extradata.(['wf',num2str(onf)]) = waveform(ds,scnl,finfo{goodf(onf),2},finfo{goodf(onf),3});
                                    
                                end
                                
                            end
                            
                        end
                        
                        % Extract data for current window
                        
                        wfdata = []; % initiate empty vector
                        
                        for ed = 1:size(goodf,1)
                            
                            if ed == size(goodf,1) % last "file" is current waveform
                                
                                temp = extract(data.(sta).dwf,'TIME',get(data.(sta).dwf,'start'),(t+twin));
                                
                            elseif ed == 1 % oldest waveform is wf1
                                
                                temp = extract(data.(sta).extradata.wf1,'TIME',t,get(data.(sta).dwf,'end'));
                                
                            else % for any "middle" files, will entire waveform
                                
                                temp = extract(data.(sta).extradata.(['wf',num2str(ed)]),'TIME',get(data.(sta).dwf,'start'),get(data.(sta).dwf,'end'));
                                
                            end
                            
                            wfdata = [wfdata;get(temp,'Data')];
                            
                        end
                        
                        data.(sta).wf = waveform(sta,cha,get(temp,'freq'),t,wfdata); % assume all files have same sampling rate
                        
                        clear temp wfdata; % clear temporary variables
                        
                        % If window end time is in current file, check for extra earlier data & get from file if needed
                        
                    elseif ~isempty(data.(sta).dwf) && get(data.(sta).dwf,'start') <= (t+twin) && get(data.(sta).dwf,'end') >= (t+twin)
                        
                        if isfield(data.(sta),'extradata') == 1 % if earlier data already exists, check it
                            
                            extranum = size(fieldnames(data.(sta).extradata),1); % find how many extra data waveforms exist
                            
                            wfdata = []; % initiate empty vector
                            
                            rmewf = {}; % empty vector for waveforms to remove
                            
                            for ed = 1:extranum
                                
                                ewf = ['wf',num2str(ed)];
                                
                                % if start time is within data "file"
                                
                                if get(data.(sta).extradata.(ewf),'start') >= t && get(data.(sta).extradata.(ewf),'end') <= t
                                    
                                    temp = extract(data.(sta).extradata.(ewf),'TIME',t,get(data.(sta).dwf,'end'));
                                    
                                    % if start time is beyond window, make note to remove it (no longer needed)
                                    
                                elseif get(data.(sta).extradata.(ewf),'end') >= t
                                    
                                    rmewf = [rmewf;ewf];
                                    
                                    temp = []; % no data to keep
                                    
                                    % if entire waveform "file" is within window, get it
                                    
                                elseif get(data.(sta).extradata.(ewf),'start') >= t && get(data.(sta).extradata.(ewf),'end') <= (t+twin)
                                    
                                    temp = extract(data.(sta).extradata.(ewf),'TIME',get(data.(sta).dwf,'start'),get(data.(sta).dwf,'end'));
                                    
                                else % in case none of those match, make sure temp is empty variable
                                    
                                    temp = [];
                                    
                                end
                                
                                wfdata = [wfdata;get(temp,'Data')];
                                
                            end
                            
                            % Lastly, add on end from current file
                            
                            temp = extract(data.(sta).dwf,'TIME',get(data.(sta).dwf,'start'),(t+twin));
                            
                            wfdata = [wfdata;get(temp,'Data')];
                            
                            % Remove unneeded extra waveforms and rename any remaining ones
                            
                            data.(sta).extradata = rmfield(data.(sta).extradata,rmewf);
                            
                            if ~isempty(fieldnames(data.(sta).extradata)) && ~isempty(rmewf) % if fields remain (and any were removed), rename them
                                
                                % make list of remaining "files"
                                
                                fieldns = fieldnames(data.(sta).extradata);
                                
                                fnums = nan(size(fieldns,1),1);
                                
                                for f = 1:size(fieldns,1)
                                    
                                    fnums(f,1) = str2double(fieldns{f,1}(3:end)); % get file/wf #
                                    
                                end
                                
                                fnums = sort(fnums); % put list in numeric order
                                
                                % loop through, rename, and remove old names
                                
                                for f = 1:size(fieldns,1)
                                    
                                    data.sta.extradata.(['wf',num2str(f)]) = data.sta.extradata.(['wf',num2str(fieldns(f,1))]);
                                    
                                    clear data.sta.extradata.(['wf',num2str(fieldns(f,1))]); % remove the "old" name
                                    
                                end
                                
                            else % remove extradata field
                                
                                data.(sta) = rmfield(data.(sta),'extradata');
                                
                            end
                            
                        else % get needed data
                            
                            for onf = 1:(size(goodf,1)-1) % all except last file
                                
                                % Get data from file
                                
                                if strcmpi(ftype,'CSS') == 1
                                    
                                    cssdata = opencss(datadir,'station',sta,'channel',cha,'files',finfo{goodf(onf),1});
                                    
                                    stt = (cssdata.file1.(sta).(cha).begtime / 24 / 3600) + datenum([1970 1 1 0 0 0]); % in days
                                    
                                    if onf == size(goodf,1) % last file is "current" file
                                        
                                        data.(sta).dwf = waveform(sta,cha,cssdata.file1.(sta).(cha).samprate,stt,cssdata.file1.(sta).(cha).amp);
                                        
                                    else
                                        
                                        data.(sta).extradata.(['wf',num2str(onf)]) =...
                                            waveform(sta,cha,cssdata.file1.(sta).(cha).samprate,stt,cssdata.file1.(sta).(cha).amp);
                                        
                                    end
                                    
                                elseif strcmpi(ftype,'BIN') == 1
                                    
                                    binid = fopen([datadir,finfo{goodf(onf),1}]);
                                    bindata.(sta).(cha).amp = fread(binid,'single');
                                    fclose(binid);
                                    
                                    % if min value is some huge number, there's missing data - replace with NaN
                                    
                                    if max(abs(bindata.(sta).(cha).amp)) > 1*10^7
                                        
                                        bindata.(sta).(cha).amp(abs(bindata.(sta).(cha).amp) > 1*10^7) = NaN;
                                        
                                    end
                                    
                                    bindata.(sta).(cha).samprate = length(bindata.(sta).(cha).amp) / 3600 / 24; % find the sampling rate (assuming file is 1 day long)
                                    
                                    if onf == size(goodf,1) % last file is "current" file
                                        
                                        data.(sta).dwf = waveform(sta,cha,bindata.(sta).(cha).samprate,finfo{goodf(onf),2},bindata.(sta).(cha).amp);
                                        
                                    else
                                        
                                        data.(sta).extradata.(['wf',num2str(onf)]) =...
                                            waveform(sta,cha,bindata.(sta).(cha).samprate,finfo{goodf(onf),2},bindata.(sta).(cha).amp);
                                        
                                    end
                                    
                                elseif strcmpi(ftype,'SAC') == 1
                                    
                                    ds = datasource('sac',[datadir,finfo{goodf(onf),1}]);
                                    
                                    scnl = scnlobject(sta,cha,params.net.list(1,:),'--');
                                    
                                    if onf == size(goodf,1) % last file is "current" file
                                        
                                        data.(sta).dwf = waveform(ds,scnl,finfo{goodf(onf),2},finfo{goodf(onf),3});
                                        
                                    else
                                        
                                        data.(sta).extradata.(['wf',num2str(onf)]) = waveform(ds,scnl,finfo{goodf(onf),2},finfo{goodf(onf),3});
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            % Extract data for current window
                            
                            wfdata = []; % initiate empty vector
                            
                            for ed = 1:size(goodf,1)
                                
                                if ed == size(goodf,1) % last "file" is current waveform
                                    
                                    temp = extract(data.(sta).dwf,'TIME',get(data.(sta).dwf,'start'),(t+twin));
                                    
                                elseif ed == 1 % oldest waveform is wf1
                                    
                                    temp = extract(data.(sta).extradata.wf1,'TIME',t,get(data.(sta).dwf,'end'));
                                    
                                else % for any "middle" files, will entire waveform
                                    
                                    temp = extract(data.(sta).extradata.(['wf',num2str(ed)]),'TIME',get(data.(sta).dwf,'start'),get(data.(sta).dwf,'end'));
                                    
                                end
                                
                                wfdata = [wfdata;get(temp,'Data')];
                                
                            end
                            
                            data.(sta).wf = waveform(sta,cha,get(temp,'freq'),t,wfdata); % assume all files have same sampling rate
                            
                            clear temp wfdata; % clear temporary variables
                            
                        end
                        
                    else % will need to open all files
                        
                        for onf = 1:size(goodf,1)
                            
                            % Get data from file
                            
                            if strcmpi(ftype,'CSS') == 1
                                
                                cssdata = opencss(datadir,'station',sta,'channel',cha,'files',finfo{goodf(onf),1});
                                
                                stt = (cssdata.file1.(sta).(cha).begtime / 24 / 3600) + datenum([1970 1 1 0 0 0]); % in days
                                
                                if onf == size(goodf,1) % last file is "current" file
                                    
                                    data.(sta).dwf = waveform(sta,cha,cssdata.file1.(sta).(cha).samprate,stt,cssdata.file1.(sta).(cha).amp);
                                    
                                else
                                    
                                    data.(sta).extradata.(['wf',num2str(onf)]) =...
                                        waveform(sta,cha,cssdata.file1.(sta).(cha).samprate,stt,cssdata.file1.(sta).(cha).amp);
                                    
                                end
                                
                            elseif strcmpi(ftype,'BIN') == 1
                                
                                binid = fopen([datadir,finfo{goodf(onf),1}]);
                                bindata.(sta).(cha).amp = fread(binid,'single');
                                fclose(binid);
                                
                                % if min value is some huge number, there's missing data - replace with NaN
                                
                                if max(abs(bindata.(sta).(cha).amp)) > 1*10^7
                                    
                                    bindata.(sta).(cha).amp(abs(bindata.(sta).(cha).amp) > 1*10^7) = NaN;
                                    
                                end
                                
                                bindata.(sta).(cha).samprate = length(bindata.(sta).(cha).amp) / 3600 / 24; % find the sampling rate (assuming file is 1 day long)
                                
                                if onf == size(goodf,1) % last file is "current" file
                                    
                                    data.(sta).dwf = waveform(sta,cha,bindata.(sta).(cha).samprate,finfo{goodf(onf),2},bindata.(sta).(cha).amp);
                                    
                                else
                                    
                                    data.(sta).extradata.(['wf',num2str(onf)]) =...
                                        waveform(sta,cha,bindata.(sta).(cha).samprate,finfo{goodf(onf),2},bindata.(sta).(cha).amp);
                                    
                                end
                                
                            elseif strcmpi(ftype,'SAC') == 1
                                
                                ds = datasource('sac',[datadir,finfo{goodf(onf),1}]);
                                
                                scnl = scnlobject(sta,cha,params.net.list(1,:),'--');
                                
                                if onf == size(goodf,1) % last file is "current" file
                                    
                                    data.(sta).dwf = waveform(ds,scnl,finfo{goodf(onf),2},finfo{goodf(onf),3});
                                    
                                else
                                    
                                    data.(sta).extradata.(['wf',num2str(onf)]) = waveform(ds,scnl,finfo{goodf(onf),2},finfo{goodf(onf),3});
                                    
                                end
                                
                            end
                            
                        end
                        
                        % Extract data for current window
                        
                        wfdata = []; % initiate empty vector
                        
                        for ed = 1:size(goodf,1)
                            
                            if ed == size(goodf,1) % last "file" is current waveform
                                
                                temp = extract(data.(sta).dwf,'TIME',get(data.(sta).dwf,'start'),(t+twin));
                                
                            elseif ed == 1 % oldest waveform is wf1
                                
                                temp = extract(data.(sta).extradata.wf1,'TIME',t,get(data.(sta).dwf,'end'));
                                
                            else % for any "middle" files, will entire waveform
                                
                                temp = extract(data.(sta).extradata.(['wf',num2str(ed)]),'TIME',get(data.(sta).dwf,'start'),get(data.(sta).dwf,'end'));
                                
                            end
                            
                            wfdata = [wfdata;get(temp,'Data')];
                            
                        end
                        
                        data.(sta).wf = waveform(sta,cha,get(temp,'freq'),t,wfdata); % assume all files have same sampling rate
                        
                        clear temp wfdata; % clear temporary variables
                        
                    end
                    
                elseif isempty(goodf) == 1 % no files found
                    
                    data.(sta).wf = []; % make sure waveform is empty
                    
                end
                
                % get window from current data waveform
                
            elseif ~isempty(data.(sta).dwf) && (get(data.(sta).dwf,'start') <= t && get(data.(sta).dwf,'end') >= (t+twin))
                
                data.(sta).wf = extract(data.(sta).dwf,'TIME',t,(t+twin));
                
                % Clear any earlier "extra" files since no longer needed
                
                if isfield(data.(sta),'extradata') == 1
                    
                    data.(sta) = rmfield(data.(sta),'extradata');
                    
                end
                
            else % something went wrong....
                
                disp('Error! Can''t find data.')
                
                data.(sta).wf = []; % make sure waveform is empty
                
            end
            
            % If not enough non-NaN data, make waveform empty
            
            if length(find(~isnan(get(data.(sta).wf,'data')),1)) > 24 %(at least 24 non-NaN samples)
                
                disp(['No data for ' sta])
                
                data.(sta).wf = [];
                
            end
            
        else
            
            % Get data from Winston waveserver
            
            scnlList = scnlobject(sta,cha,params.net.list,'--');
            data.(sta).wf = waveform(mySource,scnlList,t,(t+twin));
            
        end


        %% Run Detector if There's Data

        if ~isempty(data.(sta).wf)

            if data.(sta).gapon == 1 % if data was in long gap, turn off and record
                
                data.(sta).gap(end,2) = t; % start time of current window is gap end
                
                data.(sta).gapon = 0;
                
            end
            
            %% Filter

            data.(sta).wf = zero2nan(data.(sta).wf,10); % change zeros (min of 10 consecutive) to NaNs
            
            data.(sta).wf = demean(data.(sta).wf); % remove mean

            data.(sta).wf = fillgaps(data.(sta).wf,0); % replace NaNs with 0s before filtering

            % Make filter object

            f = filterobject('b', [lf hf], poles); % band-pass filter

            data.(sta).wf = filtfilt(f,data.(sta).wf);


            %% Do STA-LTA Filtering

            %         disp(' ')
            %         disp('Applying STA-LTA filter....')
            
            % Call STA-LTA function
            % [l_sta l_lta th_on th_off min_sep min_dur] (times in s)

            data.(sta).evcat = sta_lta(data.(sta).wf,'edp',[l_sta l_lta threson thresoff min_sep min_dur]);%,'eot','sst');

            % if event(s) found, go to next step

            if isempty(data.(sta).evcat) == 0

                % Plot event times on waveform for testing
                
                if exist('evtspl','var') == 1 && ~isempty(evtspl) && evtspl == 1
                    
                    figure;
                    
                    plot((get(data.(sta).wf,'timevector')),get(data.(sta).wf,'data'));
                    
                    hold all;
                    
                    scatter(data.(sta).evcat(:,1),ones(size(data.(sta).evcat(:,1)))*(max(get(data.(sta).wf,'data'))/2),49,'k','filled');
%                     scatter(data.(sta).evcat(:,2),ones(size(data.(sta).evcat(:,2)))*(max(get(data.(sta).wf,'data'))/2),49,'k');
                    
                    title(sta);
                    
                end
                
                
                %% Check for Matching Waveforms

                for onev = 1:size(data.(sta).evcat,1)

                    wfst = data.(sta).evcat(onev,1) - (buff / 3600 / 24); % in days
                    wfend = wfst + (evwin / 3600 / 24); % in days
                    
                    % Make sure waveform fits in current data window 
                    % Note that the STA/LTA won't pick up events within the first LTA window of
                    % the waveform, but that's enforced here anyway to avoid potential repeats
                    
                    if wfend <= get(data.(sta).wf,'end') && wfst >= get(data.(sta).wf,'start') &&...
                        wfst > (get(data.(sta).wf,'start')+(l_lta/3600/24))

                        evwf = extract(data.(sta).wf,'TIME',wfst,wfend); % get event waveform
                        
                        evwf = demean(evwf); % demean again just to make sure all event waveforms centered around zero
                        
                        % If no templates currently exist, start list

                        if ~isfield(data.(sta),'templates')

                            data.(sta).templates = evwf;

                            data.(sta).tempdata = [1 1 0 0]; % [template_ID# total_#_matches sequence_declared lv2_declared]

                            tID.(sta) = 1; % on template 1

                            data.(sta).temp1.evdata = data.(sta).evcat(onev,1); % keep event onset time

                        else % Check to see if event(s) matches any templates

                            % Put event and templates into correlation object

                            trigs = vertcat(get(data.(sta).templates,'start'),get(evwf,'start')); % trigger on waveform start times
                            
                            cobj = correlation(vertcat(data.(sta).templates,evwf),trigs);

                            % Cross-correlate

                            cobj = xcorr(cobj,[xcst xcend],'row',get(cobj,'traces')); % only xcorr with new waveform

                            % Find max correlation value (that's not with self)

                            xcmat = get(cobj,'corr'); % get correlation matrix
                            [maxcorr,tnum] = max(xcmat(length(xcmat),(1:(length(xcmat)-1))));

                            lagmat = get(cobj,'lag'); % get lag matrix

                            % if match found, stack waveform with template and update info
                            % also make sure not a duplicate

                            if maxcorr >= mincc && abs(lagmat(length(lagmat),tnum)) <= maxlag &&...
                                    isempty(find(data.(sta).evcat(onev,1)==...
                                    data.(sta).(['temp',num2str(data.(sta).tempdata(tnum,1))]).evdata(:,1),1)) == 1

                                %disp(['Match found with template ' num2str(data.(sta).tempdata(tnum,1))])

                                % Make corr object with just new waveform & matching template

                                mcobj = correlation(waveform(cobj,[tnum length(lagmat)]),trigs([tnum length(lagmat)],1));
                                
                                % Adjust trigger to matching template

                                mcobj = xcorr(mcobj,[xcst xcend]); % need to repopulate lag matrix
                                
                                mcobj = adjusttrig(mcobj,'INDEX',1);

                                % Crop to matching template window

                                mcobj = crop(mcobj,0,evwin);

                                stackcobj = waveform(stack(mcobj)); % stack matching template and new event

                                nonmatches = find((1:length(data.(sta).templates) ~= tnum)); % list of non-matching template #s

                                % Add event onset to template event list

                                tname = strcat('temp',num2str(data.(sta).tempdata(tnum,1)));

                                data.(sta).(tname).evdata = vertcat(data.(sta).(tname).evdata,data.(sta).evcat(onev,1));

                                if keepseqtemp == 1
                                   
                                    data.(sta).(tname).template = stackcobj(length(stackcobj)); % keep "new" stacked template

                                end
                                
                                % update total match count

                                data.(sta).tempdata(tnum,2) = data.(sta).tempdata(tnum,2) + 1;

                                % update template and info

                                data.(sta).templates = vertcat(stackcobj(length(stackcobj)),data.(sta).templates(nonmatches));
                                data.(sta).tempdata = vertcat(data.(sta).tempdata(tnum,:),data.(sta).tempdata(nonmatches,:)); % re-order tempdata

                            % no match found and not a duplicate, add to list
                                 
                            elseif isempty(data.(sta).evcat(onev,1)==data.(sta).(['temp',num2str(data.(sta).tempdata(tnum,1))]).evdata(:,1)) == 0

                                if size(data.(sta).templates,1) < tempmax % if still below max allowed template #

                                    tID.(sta) = tID.(sta) + 1; % move to next ID number

                                    data.(sta).templates = vertcat(evwf,data.(sta).templates);

                                    data.(sta).tempdata = vertcat([tID.(sta) 1 0 0],data.(sta).tempdata);

                                    tname = strcat('temp',num2str(tID.(sta)));

                                    data.(sta).(tname).evdata = data.(sta).evcat(onev,1); % keep event onset time

                                else % only keep (tempmax - 1) most recent templates + newest addition (and any in sequence)

                                    tID.(sta) = tID.(sta) + 1; % move to next ID number

                                    % Find last listed template without sequence "on"
                                    
                                    rmtnum = find((data.(sta).tempdata(:,3) == 0),1,'last');
                                    
                                    % remove field entry of disappearing template

                                    if data.(sta).tempdata(rmtnum,2) < minev

                                        tname = strcat('temp',num2str(data.(sta).tempdata(rmtnum,1)));
                                        data.(sta) = rmfield(data.(sta),tname);

                                    end

                                    keeptemp = find((1:length(data.(sta).templates)) ~= rmtnum); % make list of templates to keep
                                    
                                    % add new template to top of list and corresponding info

                                    data.(sta).templates = vertcat(evwf,data.(sta).templates(keeptemp,:));

                                    data.(sta).tempdata = vertcat([tID.(sta) 1 0 0],data.(sta).tempdata(keeptemp,:));

                                    tname = strcat('temp',num2str(tID.(sta)));

                                    data.(sta).(tname).evdata = data.(sta).evcat(onev,1); % keep event onset time
                                    
                                end

                            end

                        end

                    end

                end

            end


            %% Merge Highly Correlated Templates

            if isfield(data.(sta),'templates') == 1 % if there are templates

                checktemp = 1;

                while checktemp == 1

                    % x-correlate current template events

                    cc = correlation(data.(sta).templates,get(data.(sta).templates,'start'));
                    
                    cc = xcorr(cc,[xcst xcend]);

                    ccmat = get(cc,'corr'); % get correlation matrix

                    % Find any cc values above the minimum

                    ccmat = ccmat - eye(length(ccmat)); % subtract off diagonal to avoid auto-correlation values

                    lagmat = get(cc,'lag'); % get lag matrix

                    [row,col] = find((ccmat >= mincc) & (abs(lagmat) <= maxlag));

                    % If any templates correlate above min cc, combine them

                    if isempty(col)

                        checktemp = 0; % no correlated templates, so stop check

                    else % Merge first listed template with any matching, then re-check if any remain

                        % Make list of "unique" template matches (in case template N matches to more than 1)

                        ucol = unique(col);

                        % Find correlated template matches that don't include first listed template

                        if isempty(find((row((col~=ucol(1,1)),1) ~= ucol(1,1)),1))

                            checktemp = 0; % if none, no more checks needed

                        end

                        ftname = strcat('temp',num2str(data.(sta).tempdata(ucol(1,1),1)));

                        if isfield(data.(sta),ftname) == 1 % make sure template still exists

                            ucol(1,1);
                            matchtemp = find(row((col==ucol(1,1)),1) > col((col==ucol(1,1)),1)); % only consider bottom "half" of matrix

                            adjcc = adjusttrig(cc,'INDEX',ucol(1,1)); % Adjust trigger to first (most recent) template

                            % Make list of templates to stack ensuring lag
                            % time is acceptable & not duplicate

                            stackt = ucol(1,1); % start vector of template index #s

                            mergedata = data.(sta).(ftname).evdata; % start temp list of event data
                            seqstat = data.(sta).tempdata(ucol(1,1),3); % status for new template; equals 1 if any seq to stack are on
                            seqstatl2 = data.(sta).tempdata(ucol(1,1),4); % status for new template; equals 1 if any seq to stack are lv 2 on
                            
                            for m = 1:length(matchtemp)

                                tname = strcat('temp',num2str(data.(sta).tempdata(row(matchtemp(m),1),1)));
                                    
                                % Check absolute lag time
                                
                                if abs(lagmat(row(matchtemp(m),1),col(matchtemp(m),1))) <= maxlag
                                    
                                    stackt = [stackt;row(matchtemp(m),1)];

                                    % Combine event times list with first template (and sort)

                                    mergedata = sortrows(vertcat(mergedata,data.(sta).(tname).evdata));

                                    % remove event times if not enough events for sequence

                                    if data.(sta).tempdata(row(matchtemp(m)),2) < minev

                                        data.(sta) = rmfield(data.(sta),tname);

%                                         disp(['Template ' num2str(data.(sta).tempdata(row(matchtemp(m),1),1))...
%                                             ' is being consolidated.'])

                                    elseif data.(sta).tempdata(row(matchtemp(m)),3) == 1 % if sequence is on, turn off and record

                                        data.(sta).seqdect = [data.(sta).seqdect;data.(sta).tempdata(row(matchtemp(m),1),1),...
                                            (t + twin),data.(sta).tempdata(row(matchtemp(m),1),2),2];

                                        seqstat = 1; % status for new template; equals 1 if any seq to stack are on
                                        
                                        disp(['Template ' num2str(data.(sta).tempdata(row(matchtemp(m),1),1))...
                                            ' (' sta ') is currently in sequence but is being merged.'])
                                        
                                        % Write activity to log file
                                        
                                        fprintf(fID,'%s\n',['Template ' num2str(data.(sta).tempdata(row(matchtemp(m),1),1))...
                                            ' (' sta ') is currently in sequence but is being merged.']);

                                        % Remove from "on" list

                                        data.(sta).seqon(data.(sta).seqon==data.(sta).tempdata(row(matchtemp(m),1),1)) = [];

                                        % Check for lv 2 status
                                        
                                        if l2chk == 1 && data.(sta).tempdata(row(matchtemp(m)),4) == 1 % if template is lv 2
                                            
                                            data.(sta).seqonl2(data.(sta).seqonl2==data.(sta).tempdata(row(matchtemp(m),1),1)) = [];
                                            
                                            data.(sta).seqdectl2 = [data.(sta).seqdectl2;data.(sta).tempdata(row(matchtemp(m),1),1),...
                                                (t + twin),data.(sta).tempdata(row(matchtemp(m),1),2),2];
                                            
                                            seqstatl2 = 1; % status for new template; equals 1 if any seq to stack are lv 2 on
                                            
                                            disp(['Template ' num2str(data.(sta).tempdata(row(matchtemp(m),1),1))...
                                            ' (' sta ') is currently level 2 but is being merged.'])
       
                                            % Write activity to log file
                                        
                                            fprintf(fID,'%s\n',['Template ' num2str(data.(sta).tempdata(row(matchtemp(m),1),1))...
                                            ' (' sta ') is currently level 2 but is being merged.']);
                                        
                                        end
                                        
%                                     else
% 
%                                         disp(['Template ' num2str(data.(sta).tempdata(row(matchtemp(m),1),1))...
%                                             ' is being merged.'])

                                    end

                                end

                            end

                            % If at least one template match meets criteria, merge

                            if length(stackt) > 1

                                tID.(sta) = tID.(sta) + 1; % start new template; move to next ID #

                                ntname = strcat('temp',num2str(tID.(sta)));

                                % If disappearing (first) template doesn't have enough events for sequence, remove it

                                if data.(sta).tempdata(ucol(1,1),2) < minev

                                    data.(sta) = rmfield(data.(sta),ftname);

%                                     disp(['Template ' num2str(data.(sta).tempdata(ucol(1,1),1)) ' is being consolidated.'])

                                elseif data.(sta).tempdata(ucol(1,1),3) == 1 % if sequence is on, turn off and record

                                    data.(sta).seqdect = [data.(sta).seqdect;data.(sta).tempdata(ucol(1,1),1),(t + twin),...
                                        data.(sta).tempdata(ucol(1,1),2),2];

                                    disp(['Template ' num2str(data.(sta).tempdata(ucol(1,1),1))...
                                        ' (' sta ') is currently in sequence but is being merged.'])

                                    % Write activity to log file
                                        
                                    fprintf(fID,'%s\n',['Template ' num2str(data.(sta).tempdata(ucol(1,1),1))...
                                        ' (' sta ') is currently in sequence but is being merged.']);
                                    
                                    % Remove from "on" list

                                    data.(sta).seqon(data.(sta).seqon==data.(sta).tempdata(ucol(1,1),1)) = [];

                                    % Check for lv 2 status
                                    
                                    if l2chk == 1 && data.(sta).tempdata(ucol(1,1),4) == 1 % if template is lv 2
                                        
                                        data.(sta).seqonl2(data.(sta).seqonl2==data.(sta).tempdata(ucol(1,1),1)) = [];
                                        
                                        data.(sta).seqdectl2 = [data.(sta).seqdectl2;data.(sta).tempdata(ucol(1,1),1),(t + twin),...
                                            data.(sta).tempdata(ucol(1,1),2),2];
                                        
                                        disp(['Template ' num2str(data.(sta).tempdata(ucol(1,1),1))...
                                            ' (' sta ') is currently level 2 but is being merged.'])
                                        
                                        % Write activity to log file
                                        
                                        fprintf(fID,'%s\n',['Template ' num2str(data.(sta).tempdata(ucol(1,1),1))...
                                            ' (' sta ') is currently level 2 but is being merged.']);
                                        
                                    end
                                        
%                                 else
% 
%                                     disp(['Template ' num2str(data.(sta).tempdata(ucol(1,1),1))...
%                                         ' is being merged.'])

                                end

                                adjcc = crop(adjcc,0,evwin); % Crop to window of first (most recent) template

                                stackcobj = waveform(stack(adjcc,stackt)); % stack correlated templates

                                % Make entry for new template

                                data.(sta).(ntname).evdata = mergedata;
                                
                                if keepseqtemp == 1
                                
                                    data.(sta).(ntname).template = stackcobj(length(stackcobj));

                                end
                                
                                % Make list of non-correlated templates

                                nonmatch = nan((length(data.(sta).templates)-length(stackt)),1);
                                ind = 1;

                                for n = 1:length(data.(sta).templates)

                                    if isempty(find(stackt==n,1)) % if not on stack list, add

                                        nonmatch(ind,1) = n;
                                        ind = ind + 1; % nove to next row

                                    end

                                end

                                % Add new template to top and remove old templates (and associated info)

                                data.(sta).templates = vertcat(stackcobj(end),data.(sta).templates(nonmatch));

                                data.(sta).tempdata = vertcat([tID.(sta) size(data.(sta).(ntname).evdata,1) seqstat seqstatl2],...
                                    data.(sta).tempdata(nonmatch,:));

                                % Add entry to sequence timeline and sequence on lists if it's turned on

                                if seqstat == 1

                                    data.(sta).seqdect = [data.(sta).seqdect;data.(sta).tempdata(1,1),(t + twin),data.(sta).tempdata(1,2),3];

                                    data.(sta).seqon = [data.(sta).seqon; data.(sta).tempdata(1,1)];

                                    disp(['Merged template ' num2str(data.(sta).tempdata(1,1)) ' (' sta ') is being turned on.'])

                                    % Write activity to log file
                                        
                                    fprintf(fID,'%s\n',['Merged template ' num2str(data.(sta).tempdata(1,1))...
                                        ' (' sta ') is being turned on.']);
                                    
                                end

                                % Add entry to sequence timeline and sequence on lists if it's turned on lv 2
                                
                                if l2chk == 1 && seqstatl2 == 1

                                    data.(sta).seqdectl2 = [data.(sta).seqdectl2;data.(sta).tempdata(1,1),(t + twin),data.(sta).tempdata(1,2),3];

                                    data.(sta).seqonl2 = [data.(sta).seqonl2; data.(sta).tempdata(1,1)];

                                    disp(['Merged template ' num2str(data.(sta).tempdata(1,1)) ' (' sta ') is being turned on for level 2.'])

                                    % Write activity to log file
                                        
                                    fprintf(fID,'%s\n',['Merged template ' num2str(data.(sta).tempdata(1,1))...
                                        ' (' sta ') is being turned on for level 2.']);
                                    
                                end
                                
                            end

                        end

                        % clear temporary variables

                        clear adjcc cc stackobj ccmat row col lagmat m ucol stackt ftname tname mergedata ind nonmatch seqstat;

                    end

                end

            end
            
        else % keep track of gap times
            
            if data.(sta).gapon == 0 % if not already in a gap, initiate and get time
                
                data.(sta).gap = [data.(sta).gap;t,(t+twin)]; % use start and end times of this window for now
                
                data.(sta).gapon = 1;
                
            end
            
        end

        
        %% Check for Event Sequences
        
        % if there is template data (i.e. at least one event)
        
        if isfield(data.(sta),'tempdata')
            
            evseq = find(data.(sta).tempdata(:,2)>=minev); % look for event sequences
            
            if isempty(evseq) == 0 % if at least one template has enough events for sequence
%
%                 disp(' ')
%                 disp(['Current window end time: ' datestr(t+twin)])
                
                for row = 1:length(evseq)
                    
                    tname = strcat('temp',num2str(data.(sta).tempdata(evseq(row),1)));
                    
                    % Check if enough events within given time period (relative to current time)
                    
                    interevt = find(data.(sta).(tname).evdata(:,1) >= ((t + twin) - seqT));
                    interevtoff = find(data.(sta).(tname).evdata(:,1) >= ((t + twin) - seqToff - (twin * seqoffwin)));
                    
                    % if enough matches and not yet declared, declare sequence
                    
                    if length(interevt) >= minev && data.(sta).tempdata(evseq(row),3) == 0
                        
                        data.(sta).tempdata(evseq(row),3) = 1; % sequence on
                        
                        % display 1st event, most recent event, total # of events, on/off, station
                        
                        disp([' ' datestr(data.(sta).(tname).evdata(1,1)) '      ' datestr(max(data.(sta).(tname).evdata(:,1))) '            '...
                            num2str(data.(sta).tempdata(evseq(row),2)) '                 '...
                            num2str(data.(sta).tempdata(evseq(row),1)) '          on        ' sta])
                        
                        % Write to log file
                        
                        fprintf(fID,'%s\n',[' ' datestr(data.(sta).(tname).evdata(1,1)) '      ' datestr(max(data.(sta).(tname).evdata(:,1))) '            '...
                            num2str(data.(sta).tempdata(evseq(row),2)) '                 '...
                            num2str(data.(sta).tempdata(evseq(row),1)) '          on        ' sta]);
                        
                        % Make list of sequence on/off times for plotting
                        
                        if ~isfield(data.(sta),'seqdect')
                            
                            % [template_ID current_window_start_time total_events on/off] with on=1, off=0
                            data.(sta).seqdect = [data.(sta).tempdata(evseq(row),1),(t+twin),data.(sta).tempdata(evseq(row),2),1];
                            
                        else
                            
                            data.(sta).seqdect = [data.(sta).seqdect;data.(sta).tempdata(evseq(row),1),(t+twin),...
                                data.(sta).tempdata(evseq(row),2),1];
                            
                        end
                        
                        % Make list of sequences currently on
                        
                        data.(sta).seqon = [data.(sta).seqon; data.(sta).tempdata(evseq(row),1)];
                        
                        % if not enough matches and already declared, turn "off" sequence
                        
                    elseif length(interevtoff) < minevoff && data.(sta).tempdata(evseq(row),3) == 1
                        
                        data.(sta).tempdata(evseq(row),3) = 0; % sequence off
                        
                        % display 1st event, most recent event, total # of events, on/off, station
                        
                        disp([' ' datestr(data.(sta).(tname).evdata(1,1)) '      ' datestr(max(data.(sta).(tname).evdata(:,1))) '            '...
                            num2str(data.(sta).tempdata(evseq(row),2)) '                 '...
                            num2str(data.(sta).tempdata(evseq(row),1)) '          off       ' sta])
                        
                        % Write to log file
                        
                        fprintf(fID,'%s\n',[' ' datestr(data.(sta).(tname).evdata(1,1)) '      ' datestr(max(data.(sta).(tname).evdata(:,1))) '            '...
                            num2str(data.(sta).tempdata(evseq(row),2)) '                 '...
                            num2str(data.(sta).tempdata(evseq(row),1)) '          off       ' sta]);
                        
                        % Make list of sequence on/off times for plotting
                        
                        data.(sta).seqdect = [data.(sta).seqdect;data.(sta).tempdata(evseq(row),1),...
                            (t + twin),data.(sta).tempdata(evseq(row),2),0];
                        
                        % Remove from "on" list
                        
                        data.(sta).seqon(data.(sta).seqon==data.(sta).tempdata(evseq(row),1)) = [];
                        
                    end
                    
                    % Check for level 2 sequences on/off if wanted
                    
                    if l2chk == 1
                        
                        % Check if enough events within given time period (relative to current time)
                        
                        interevtl2 = find(data.(sta).(tname).evdata(:,1) >= ((t + twin) - seqTl2));
                        interevtoffl2 = find(data.(sta).(tname).evdata(:,1) >= ((t + twin) - seqToffl2));
                        
                        % if sequence on but not yet lv 2, turn on if requirements met
                        
                        if length(interevtl2) >= minevl2 && data.(sta).tempdata(evseq(row),3) == 1 && data.(sta).tempdata(evseq(row),4) == 0
                            
                            data.(sta).tempdata(evseq(row),4) = 1; % level 2 on
                            
                            % display 1st event, most recent event, total # of events, on/off, station
                            
                            disp([' ' datestr(data.(sta).(tname).evdata(1,1)) '      ' datestr(max(data.(sta).(tname).evdata(:,1))) '            '...
                                num2str(data.(sta).tempdata(evseq(row),2)) '                 '...
                                num2str(data.(sta).tempdata(evseq(row),1)) '          on 2      ' sta])
                            
                            % Write to log file
                            
                            fprintf(fID,'%s\n',[' ' datestr(data.(sta).(tname).evdata(1,1)) '      ' datestr(max(data.(sta).(tname).evdata(:,1))) '            '...
                                num2str(data.(sta).tempdata(evseq(row),2)) '                 '...
                                num2str(data.(sta).tempdata(evseq(row),1)) '          on 2      ' sta]);
                            
                            % Make list of sequence lv 2 on/off times for plotting
                            
                            if ~isfield(data.(sta),'seqdectl2')
                                
                                % [template_ID current_window_start_time total_events lv2_on/off] with on=1, off=0
                                data.(sta).seqdectl2 = [data.(sta).tempdata(evseq(row),1),(t+twin),data.(sta).tempdata(evseq(row),2),1];
                                
                            else
                                
                                data.(sta).seqdectl2 = [data.(sta).seqdectl2;data.(sta).tempdata(evseq(row),1),(t+twin),...
                                    data.(sta).tempdata(evseq(row),2),1];
                                
                            end
                            
                            % Make list of sequences currently on lv 2
                            
                            data.(sta).seqonl2 = [data.(sta).seqonl2; data.(sta).tempdata(evseq(row),1)];
                            
                        % if lv 2 requirements no longer met, turn off
                            
                        elseif length(interevtoffl2) < minevoffl2 && data.(sta).tempdata(evseq(row),4) == 1
                            
                            data.(sta).tempdata(evseq(row),4) = 0; % lv 2 off
                            
                            % display 1st event, most recent event, total # of events, on/off, station
                            
                            disp([' ' datestr(data.(sta).(tname).evdata(1,1)) '      ' datestr(max(data.(sta).(tname).evdata(:,1))) '            '...
                                num2str(data.(sta).tempdata(evseq(row),2)) '                 '...
                                num2str(data.(sta).tempdata(evseq(row),1)) '          off 2     ' sta])
                            
                            % Write to log file
                            
                            fprintf(fID,'%s\n',[' ' datestr(data.(sta).(tname).evdata(1,1)) '      ' datestr(max(data.(sta).(tname).evdata(:,1))) '            '...
                                num2str(data.(sta).tempdata(evseq(row),2)) '                 '...
                                num2str(data.(sta).tempdata(evseq(row),1)) '          off 2     ' sta]);
                            
                            % Make list of sequence level 2 on/off times for plotting
                            
                            data.(sta).seqdectl2 = [data.(sta).seqdectl2;data.(sta).tempdata(evseq(row),1),...
                                (t + twin),data.(sta).tempdata(evseq(row),2),0];
                            
                            % Remove from level 2 "on" list
                            
                            data.(sta).seqonl2(data.(sta).seqonl2==data.(sta).tempdata(evseq(row),1)) = [];
                            
                            % if seq is turned off but lv 2 is still on, turn it off
                            
                        elseif data.(sta).tempdata(evseq(row),3) == 0 && data.(sta).tempdata(evseq(row),4) == 1
                            
                            data.(sta).tempdata(evseq(row),4) = 0; % lv 2 off
                            
                            % display 1st event, most recent event, total # of events, on/off, station
                            
                            disp([' ' datestr(data.(sta).(tname).evdata(1,1)) '      ' datestr(max(data.(sta).(tname).evdata(:,1))) '            '...
                                num2str(data.(sta).tempdata(evseq(row),2)) '                 '...
                                num2str(data.(sta).tempdata(evseq(row),1)) '          off 2     ' sta])
                            
                            % Write to log file
                            
                            fprintf(fID,'%s\n',[' ' datestr(data.(sta).(tname).evdata(1,1)) '      ' datestr(max(data.(sta).(tname).evdata(:,1))) '            '...
                                num2str(data.(sta).tempdata(evseq(row),2)) '                 '...
                                num2str(data.(sta).tempdata(evseq(row),1)) '          off 2     ' sta]);
                            
                            % Make list of sequence level 2 on/off times for plotting
                            
                            data.(sta).seqdectl2 = [data.(sta).seqdectl2;data.(sta).tempdata(evseq(row),1),...
                                (t + twin),data.(sta).tempdata(evseq(row),2),0];
                            
                            % Remove from level 2 "on" list
                            
                            data.(sta).seqonl2(data.(sta).seqonl2==data.(sta).tempdata(evseq(row),1)) = [];
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end

    end
            
    
    %% Send Alert if Needed
    
    % Count how many stations have a sequence "on"
    
    onct = 0; % always (re)start at 0
    reqct = 0; % always (re)start at 0
    
    onctl2 = 0; % always (re)start at 0
    reqctl2 = 0; % always (re)start at 0
    
%     lv2staon = []; % string list of stations for alert message
    
    for ss = 1:size(params.sta.list,1)
        
        sta = strcat(params.sta.list(ss,:)); % get station name without whitespaces
        
        if ~isempty(data.(sta).seqon) % check station for sequence(s) on
            
            onct = onct + 1;
            
            % Check to see if it's on preferred list (if there's preferred station requirement)
            
            if minreq > 0
                
                for rr = 1:size(reqsta,1)
                    
                    if strcmpi(strcat(reqsta(rr,:)),sta) == 1 % if it is, count it
                        
                        reqct = reqct + 1;
                        
                        break; % leave loop
                        
                    end
                    
                end
                
            end
            
        end
        
        if l2chk == 1 && ~isempty(data.(sta).seqonl2) % check station for sequence(s) on level 2
            
            onctl2 = onctl2 + 1;
            
%             lv2staon = [lv2staon,sta,', ']; % make a list of "on" stations
            
            % Check to see if it's on preferred list (if there's preferred station requirement)
            
            if minreq > 0
                
                for rr = 1:size(reqsta,1)
                    
                    if strcmpi(strcat(reqsta(rr,:)),sta) == 1 % if it is, count it
                        
                        reqctl2 = reqctl2 + 1;
                        
                        break; % leave loop
                        
                    end
                    
                end
                
            end
            
        end

    end
    
    % If at least minsta stations have a sequence on, send alert (if new seq)
    
    if onct >= minsta && reqct >= minreq && curstat == 0
        
        disp(' ')
        disp(['** ' num2str(onct) ' stations (' num2str(reqct) ' required) have at least one sequence in progress!'])
        disp(['Current time: ' datestr(t+twin)])
        
        % Write to log file
        
        fprintf(fID,'%s\n',['** ' num2str(onct) ' stations (' num2str(reqct)...
            ' required) have at least one sequence in progress at ' datestr(t+twin) '!']);
        
        curstat = 1;
        
        alertt.lv1.on = [alertt.lv1.on;(t+twin),onct];
        
        alertt.lv1.curon = t + twin; % latest time alert is on
        
    elseif onct >= minsta && reqct >= minreq && curstat == 1 % if alert still valid, keep track of time
        
        alertt.lv1.curon = t + twin; % latest time alert is on
        
    elseif  (onct < minsta || reqct < minreq) && curstat == 1 % else, send off alert if dropped below min (and on)
        
        % if alert has been "on" for required period of time and is past waiting time to turn off
       
        if (alertt.lv1.curon + (twin * offwait)) < (t + twin) && ((t+twin) - alertt.lv1.on(end,1)) > alertont
        
            disp(' ')
            disp('** Stations with at least one sequence in progress no longer meet requirements.')
            disp(['Current time: ' datestr(t+twin)])
            
            % Write to log file
        
            fprintf(fID,'%s\n',['** Stations with at least one sequence in progress no longer meet the minimum at ' datestr(t+twin) '.']);
            
            curstat = 0;
            
            alertt.lv1.off = [alertt.lv1.off;(t+twin),onct];
            
            alertt.lv1.curon = []; % clear current on time
            
        end
     
    end
    
    % If at least minsta stations have a sequence on level 2, send alert (if new level up)
    
    if l2chk == 1
        
        if onctl2 >= minsta && reqctl2 >= minreq && curstatl2 == 0 && curstat == 1
            
            disp(' ')
            disp(['** ' num2str(onctl2) ' stations (' num2str(reqctl2) ' required) have at least one sequence at level 2!'])
            disp(['Current time: ' datestr(t+twin)])
            
            % Write to log file
            
            fprintf(fID,'%s\n',['** ' num2str(onctl2) ' stations (' num2str(reqctl2)...
                ' required) have at least one sequence at level 2 at ' datestr(t+twin) '!']);
            
            curstatl2 = 1;
            
            alertt.lv2.on = [alertt.lv2.on;(t+twin),onctl2];
            
%             % Make figure for alert
%         
%             figname = make_alert_fig(params,directory,'lv2_ON',(t+twin),seqTl2,mySource);
%         
%             % Send email/text notification
%             
%             sendmail(sendtoonl2,'Event Sequence at level 2!',...
%                 [num2str(onctl2) ' stations have at least one sequence at level 2! Time is ' datestr(t+twin) ' UTC. '...
%                 'Sequence parameters: min ' num2str(minevl2) ' events in ' num2str(seqTl2*24*60) ' minutes.' 10 10 ...
%                 'Stations on level 2: ' lv2staon(1:end-2)],figname);
                    
             alertt.lv2.curon = t + twin; % latest time alert is on

         elseif onctl2 >= minsta && reqctl2 >= minreq && curstatl2 == 1 && curstat == 1 % if alert still valid, keep track of time

             alertt.lv2.curon = t + twin; % latest time alert is on
        
         % else, send off alert if dropped below min (and on); also, force off if alert went off
            
         elseif  (onctl2 < minsta || reqctl2 < minreq || curstat == 0) && curstatl2 == 1
             
             % if alert has been "on" for required period of time or if level 1 turns off
             
             if (alertt.lv2.curon + (twin * offwaitl2)) < (t + twin) && ((t+twin) - alertt.lv2.on(end,1)) > alertontl2 || curstat == 0 
                 
                 disp(' ')
                 disp('** Stations with at least one sequence at level 2 no longer meet requirements.')
                 disp(['Current time: ' datestr(t+twin)])
                 
                 % Write to log file
        
                 fprintf(fID,'%s\n',['** Stations with at least one sequence at level 2 no longer meet requirements at ' datestr(t+twin) '.']);
                 
                 curstatl2 = 0;
                 
                 alertt.lv2.off = [alertt.lv2.off;(t+twin),onctl2];
                 
                 alertt.lv2.curon = []; % clear current on time
                 
             end
           
        end
        
    end
     
end

fclose(fID); % close log file


%% Plot Timeline & Sequence(s) on Waveform

% Make sure at least one station has sequences to plot

goodseq = 1;

for s = 1:size(params.sta.list,1)
    
    sta = strcat(params.sta.list(s,:)); % get station name without whitespaces
    
    if isfield(data.(sta),'seqdect') % if there are sequences
        
        goodseq = 1;
        
        break; % at least one stations has sequences
        
    end
    
end

if goodseq == 1
    
    % Plot waveform over analyzed time interval (with CSS data)
    
%     sampst = 1;
%     sampend = floor((tend - tstrt) * 24 * 3600 * cssdata.file1.(sta).(cha).samprate);
%     
%     figure;
%     
%     wft = (cssdata.file1.(sta).(cha).time(sampst:sampend) - cssdata.file1.(sta).(cha).begtime) / 3600; % in hours
%     
%     plot(wft,cssdata.file1.(sta).(cha).amp(sampst:sampend));
%     
%     hold all;
%             
%         scatter(ontimes,(ones(length(ont),1)*max(cssdata.file1.(sta).(cha).amp)/2),25,data.(sta).seqdect(ont,1),'filled'); % on times
%         scatter(offtimes,(ones(length(offt),1)*max(cssdata.file1.(sta).(cha).amp)/2),25,data.(sta).seqdect(offt,1)); % off times
    
    % Make timeline plot for each station
    
    if plstatl == 1
        
        % Read in known slow clap times
        
        if exist('scfile','var') == 1
            
            sc_evs = importdata(scfile); % in UTC
            
        end
        
        onsp = 0;
        
        figure;
        
        for s = 1:size(params.sta.list,1)
            
            sta = strcat(params.sta.list(s,:)); % get station name without whitespaces
            
            cha = params.cha.list(s,:);
            
            % Run plotting if station has sequences
            
            if isfield(data.(sta),'seqdect') == 1
            
            % Find number of unique templates in seqdect
                        
            seqs = unique(data.(sta).seqdect(:,1));
            
            onsp = onsp + 1;
            
            ons=find(data.(sta).seqdect(:,4)==1); % make list of on's
            offs=find(data.(sta).seqdect(:,4)==0); % make list of off's
            cons=find(data.(sta).seqdect(:,4)==3); % make list of consolidation on's
            coffs=find(data.(sta).seqdect(:,4)==2); % make list of consolidation off's
            
            ontimes=(data.(sta).seqdect(ons,2)-tstrt) * 24; % find on times for plotting
            offtimes=(data.(sta).seqdect(offs,2)-tstrt) * 24; % find off times for plotting
            contimes=(data.(sta).seqdect(cons,2)-tstrt) * 24; % find consolidation on times for plotting
            cofftimes=(data.(sta).seqdect(coffs,2)-tstrt) * 24; % find consolidation off times for plotting
            
            % Make vector with colors
            
            oncol = nan(length(ons),1);
            
            for a = 1:length(ons)
                
                oncol(a,1) = find(seqs == data.(sta).seqdect(ons(a),1));
                
            end
            
            offcol = nan(length(offs),1);
            
            for a = 1:length(offs)
                
                offcol(a,1) = find(seqs == data.(sta).seqdect(offs(a),1));
                
            end
            
            concol = nan(length(cons),1);
            
            for a = 1:length(cons)
                
                concol(a,1) = find(seqs == data.(sta).seqdect(cons(a),1));
                
            end
            
            coffcol = nan(length(coffs),1);
            
            for a = 1:length(coffs)
                
                coffcol(a,1) = find(seqs == data.(sta).seqdect(coffs(a),1));
                
            end
            
            % Plot sequence on/off times
            
            subplot(4,2,onsp);
            
            colormap(colorcube(length(seqs)));
                             
            hold all;
            
            scatter(ontimes,data.(sta).seqdect(ons,3),100,oncol,'filled','marker','^','markeredgecolor','k'); % on times
            scatter(offtimes,data.(sta).seqdect(offs,3),100,offcol,'filled','marker','o','markeredgecolor','k'); % off times
            scatter(contimes,data.(sta).seqdect(cons,3),25,concol,'filled','marker','d','markeredgecolor','k'); % consolidation on times
            scatter(cofftimes,data.(sta).seqdect(coffs,3),25,coffcol,'filled','marker','s','markeredgecolor','k'); % consolidation off times
            
            xlabel(['Hours since ',datestr(tstrt)]);
            ylabel('# of Events in Sequence/Template');
            title([{['Timeline of Sequences for ',sta,' ',cha]},{'\^ = on    o = off    <> = merger on    | | = merger off'}]);
            
            xlim([0 ((tend-tstrt)*24)]); % plot over entire time period examined
            
            c = colorbar('location','eastoutside');
            
            c.Label.String = 'Sequence/Template #';
                        
            set(gca,'fontsize',14);
            box on;
            grid on;
            
            hold all;
            
            % Plot alert on/off times
            
            if ~isempty(alertt.lv1.on)
                
                scatter(((alertt.lv1.on(:,1)-tstrt)*24),alertt.lv1.on(:,2),144,'r','filled','marker','p'); % alert on times
                
            end
            
            if ~isempty(alertt.lv1.off)
                
                scatter(((alertt.lv1.off(:,1)-tstrt)*24),alertt.lv1.off(:,2),144,'b','filled','marker','p'); % alert off times
                
            end
            
            % Plot lines of known sequence on/off times
            
            if exist('sc_evs','var') == 1
                
                for onev = 1:size(sc_evs,1)
                    
                    plot(([sc_evs(onev,1) sc_evs(onev,1)] - tstrt)*24,[4 80],'r'); % slow clap start
                    plot(([sc_evs(onev,2) sc_evs(onev,2)] - tstrt)*24,[4 80],'b'); % slow clap end
                    
                end
                
            end
            
            end
        
        end
        
    end
    
    % Plot sequence/alert timeline for all stations
    
    if pltl == 1
        
        % Read in known slow clap times if it hasn't been done already
        
        if exist('scfile','var') == 1 && ~exist('sc_evs','var')
            
            sc_evs = importdata(scfile); % in UTC
            
        end
         
        % Read in known explosion times if it hasn't been done already
        
        if exist('expfile','var') == 1 && ~exist('exp_evs','var')
            
            exp_evs = importdata(expfile); % in UTC
            
        end
        
        % Set up subplots (each up to a week)
        
        spls = ceil((tend - tstrt) / 7); % total # of subplots needed
            
        % Plot ons/offs for sequences and alerts on same plot
        
        h=figure;
                    
        onspl = 1; % start on first subplot
        
        if ~exist('splrows','var') || (tend - tstrt) <= 7
            
            splrows = 1; % default to one row/plot per fig if not defined
           
        elseif splrows > spls % if there aren't enough weeks to fill all subplots on fist figure
            
            splrows = spls; % reset to only have as many rows as needed
            
        end
            
        if (tend - tstrt) <= 7 % is time period is <= 1 week
            
           tstop = tstrt + 6; % make a "fake" tstop
           
        else
            
            tstop = tend; % stop at end of time period
            
        end
        
        for tt = tstrt:7:tstop % loop through time period in weeks

        for s = 1:size(params.sta.list,1)

            sta = strcat(params.sta.list(s,:)); % get station name without whitespaces
                                
            subplot(splrows,1,onspl);
                
            % Plot long data gaps as horizontal lines (if there are any)
            
            if ~isempty(data.(sta).gap)
                
                % Find gaps during current window
                
                gaps = find((data.(sta).gap(:,1) >= tt & data.(sta).gap(:,1) <= tt+7) |...
                    (data.(sta).gap(:,2) >= tt & data.(sta).gap(:,2) <= tt+7) |...
                    (data.(sta).gap(:,1) <= tt & data.(sta).gap(:,2) >= tt+7));
                                    
                hold all;
                    
                for gr = 1:length(gaps)
                    
                    stt = max(data.(sta).gap(gaps(gr),1),tt); % either the start time or the subplot start time
                    edt = min(data.(sta).gap(gaps(gr),2),(tt+7)); % either the start time or the subplot start time

                    plot([stt edt],[s s],'color',[0 0.75 0],'linewidth',1.5); % current station # (alert line is y=0)
                    
                end
                
            end
            
            % Plot station outages as horizontal lines (if there are any)
            
            % Read in outages from file
            
            if exist('staoutfile','var') == 1 && ~exist('sta_out','var')
                
                sta_out = importdata(staoutfile); % in UTC
                
            end
        
            if exist('sta_out','var')
                
                % Find outages during current window
                
                outs = find((cell2mat(sta_out(:,2)) >= tt & cell2mat(sta_out(:,2)) <= tt+7) |...
                    (cell2mat(sta_out(:,3)) >= tt & cell2mat(sta_out(:,3)) <= tt+7) |...
                    (cell2mat(sta_out(:,2)) <= tt & cell2mat(sta_out(:,3)) >= tt+7));
                
                % Look through list for current station
                
                souts = [];
                
                for so = 1:size(outs,1)
                    
                    if strcmpi(strcat(sta_out{outs(so),1}),sta) == 1
                        
                        souts = [souts;outs(so)]; % list of rows with relevant info
                        
                    end
                    
                end
                
                if ~isempty(souts) % if there are outages, plot them
                    
                    hold all;
                    
                    for sso = 1:length(souts)
                        
                        stt = max(sta_out{souts(sso),2},tt); % either the start time or the subplot start time
                        edt = min(sta_out{souts(sso),3},(tt+7)); % either the start time or the subplot start time
                        
                        plot([stt edt],[s s],'color',[0 0.5 0],'linewidth',1.5); % current station # (alert line is y=0)
                        
                    end
                    
                end
                
            end
            
            % Run plotting if station has sequences within current time period
            
            if isfield(data.(sta),'seqdect') == 1 && ~isempty(find(data.(sta).seqdect(:,2)>=tt & data.(sta).seqdect(:,2)<=(tt+7),1))
                
                % Make lists of relevant rows within current time period
                
                ons = find(data.(sta).seqdect(:,4)==1 & data.(sta).seqdect(:,2)>=tt & data.(sta).seqdect(:,2)<=(tt+7)); % on's
                offs = find(data.(sta).seqdect(:,4)==0 & data.(sta).seqdect(:,2)>=tt & data.(sta).seqdect(:,2)<=(tt+7)); % off's
                cons = find(data.(sta).seqdect(:,4)==3 & data.(sta).seqdect(:,2)>=tt & data.(sta).seqdect(:,2)<=(tt+7)); % merger on's
                coffs = find(data.(sta).seqdect(:,4)==2 & data.(sta).seqdect(:,2)>=tt & data.(sta).seqdect(:,2)<=(tt+7)); % merger off's
                
                % Find times for plotting
                
%                 ontimes=(data.(sta).seqdect(ons,2)-tstrt) * 24; % on
%                 offtimes=(data.(sta).seqdect(offs,2)-tstrt) * 24; % off
%                 contimes=(data.(sta).seqdect(cons,2)-tstrt) * 24; % merger on
%                 cofftimes=(data.(sta).seqdect(coffs,2)-tstrt) * 24; % merger off

                % Find times and plot
                
                hold all;             
                
                if ~isempty(ons)
                    
                    ontimes = data.(sta).seqdect(ons,2); % on
                    scatter(ontimes,ones(size(ontimes))*s,100,'k','filled','marker','^','markeredgecolor','k'); % on times
                    
                end
                
                if ~isempty(offs)
                    
                    offtimes = data.(sta).seqdect(offs,2); % off
                    scatter(offtimes,ones(size(offtimes))*s,100,'k','filled','marker','o','markeredgecolor','k'); % off times
                    
                end
                
                if ~isempty(cons)
                    
                    contimes = data.(sta).seqdect(cons,2); % merger on
                    scatter(contimes,ones(size(contimes))*s,25,'k','filled','marker','d','markeredgecolor','k'); % merger on times
                    
                end
                
                if ~isempty(coffs)
                    
                    cofftimes = data.(sta).seqdect(coffs,2); % merger off
                    scatter(cofftimes,ones(size(cofftimes))*s,25,'k','filled','marker','s','markeredgecolor','k'); % merger off times
                    
                end
                
                % Do the same for level 2 if it exists and has sequences within current time period
                
                if l2chk == 1 && isfield(data.(sta),'seqdectl2') == 1 &&...
                        ~isempty(find(data.(sta).seqdectl2(:,2)>=tt & data.(sta).seqdectl2(:,2)<=(tt+7),1))
                    
                    % Make lists of relevant rows within current time period
                    
                    ons = find(data.(sta).seqdectl2(:,4)==1 & data.(sta).seqdectl2(:,2)>=tt & data.(sta).seqdectl2(:,2)<=(tt+7)); % on's
                    offs = find(data.(sta).seqdectl2(:,4)==0 & data.(sta).seqdectl2(:,2)>=tt & data.(sta).seqdectl2(:,2)<=(tt+7)); % off's
                    cons = find(data.(sta).seqdectl2(:,4)==3 & data.(sta).seqdectl2(:,2)>=tt & data.(sta).seqdectl2(:,2)<=(tt+7)); % merger on's
                    coffs = find(data.(sta).seqdectl2(:,4)==2 & data.(sta).seqdectl2(:,2)>=tt & data.(sta).seqdectl2(:,2)<=(tt+7)); % merger off's
                    
                    % Find times for plotting
                    
%                     ontimes=(data.(sta).seqdectl2(ons,2)-tstrt) * 24; % on
%                     offtimes=(data.(sta).seqdectl2(offs,2)-tstrt) * 24; % off
%                     contimes=(data.(sta).seqdectl2(cons,2)-tstrt) * 24; % merger on
%                     cofftimes=(data.(sta).seqdectl2(coffs,2)-tstrt) * 24; % merger off

                    % Find times and plot
                    
                    hold all;             
                    
                    if ~isempty(ons)
                        
                        ontimes = data.(sta).seqdectl2(ons,2); % on
                        scatter(ontimes,ones(size(ontimes))*s,100,'y','filled','marker','^','markeredgecolor','k'); % on times
                        
                    end
                    
                    if ~isempty(offs)
                        
                        offtimes = data.(sta).seqdectl2(offs,2); % off
                        scatter(offtimes,ones(size(offtimes))*s,100,'y','filled','marker','o','markeredgecolor','k'); % off times
                        
                    end
                    
                    if ~isempty(cons)
                        
                        contimes = data.(sta).seqdectl2(cons,2); % merger on
                        scatter(contimes,ones(size(contimes))*s,25,'y','filled','marker','d','markeredgecolor','k'); % merger on times
                        
                    end
                    
                    if ~isempty(coffs)
                        
                        cofftimes = data.(sta).seqdectl2(coffs,2); % merger off
                        scatter(cofftimes,ones(size(cofftimes))*s,25,'y','filled','marker','s','markeredgecolor','k'); % merger off times
                        
                    end
                    
                end
                
                % Plot alert on/off times
                
                if ~isempty(alertt.lv1.off)
                    
                    % Find off alerts within current time period
                    
                    offs = find(alertt.lv1.off(:,1)>=tt & alertt.lv1.off(:,1)<=(tt+7));
                    
%                     scatter(((alertt.lv1.off(:,1)-tstrt)*24),zeros(size(alertt.lv1.off,1),1),196,'b','filled','marker','p'); % alert off times
                    
                    if ~isempty(offs)
     
                        scatter(alertt.lv1.off(offs,1),zeros(length(offs),1),196,'b','filled','marker','p','markeredgecolor','k'); % alert off times
                   
                    end

                end
                
                if ~isempty(alertt.lv1.on)
                    
                    % Find on alerts within current time period
                    
                    ons = find(alertt.lv1.on(:,1)>=tt & alertt.lv1.on(:,1)<=(tt+7));
                    
%                     scatter(((alertt.lv1.on(:,1)-tstrt)*24),zeros(size(alertt.lv1.on,1),1),196,'r','filled','marker','p'); % alert on times

                    if ~isempty(ons)
                        
                        scatter(alertt.lv1.on(ons,1),zeros(length(ons),1),196,'r','filled','marker','p','markeredgecolor','k'); % alert on times
                    
                    end
                    
                end
                
                if ~isempty(alertt.lv2.on)
                    
                    % Find alerts within current time period
                    
                    ons = find(alertt.lv2.on(:,1)>=tt & alertt.lv2.on(:,1)<=(tt+7)); % on's
                    offs = find(alertt.lv2.off(:,1)>=tt & alertt.lv2.off(:,1)<=(tt+7)); % off's
                    
%                     scatter(((alertt.lv2.on(:,1)-tstrt)*24),zeros(size(alertt.lv2.on,1),1),196,'m','filled','marker','p'); % alert lv 2 on times
%                     scatter(((alertt.lv2.off(:,1)-tstrt)*24),zeros(size(alertt.lv2.off,1),1),196,'c','filled','marker','p'); % alert lv 2 off times
                    
                    if ~isempty(offs)
                        
                        scatter(alertt.lv2.off(offs,1),zeros(length(offs),1),196,'c','filled','marker','p','markeredgecolor','k'); % alert lv 2 off times
                    
                    end

                    if ~isempty(ons)
                        
                        scatter(alertt.lv2.on(ons,1),zeros(length(ons),1),196,[1 0.75 0.75],'filled','marker','p',...
                            'markeredgecolor','k'); % alert lv 2 on times

                    end
                    
                end
                
            end
            
        end
       
         % Plot lines of known sequence on/off times
         
         hold all;
         
         if exist('sc_evs','var') == 1
             
             % Find sequences within current time period
             
             seqrows = find((sc_evs(:,1)>=tt & sc_evs(:,1)<=(tt+7)) | (sc_evs(:,2)>=tt & sc_evs(:,2)<=(tt+7))); % either start or stop in period
             
             if ~isempty(seqrows) % if there are some, plot them
                 
                 for onev = 1:length(seqrows)
                     
                     %plot(([sc_evs(seqrows(onev),1) sc_evs(seqrows(onev),1)] - tstrt)*24,[0 s],'r'); % slow clap start
                     %plot(([sc_evs(seqrows(onev),2) sc_evs(seqrows(onev),2)] - tstrt)*24,[0 s],'b'); % slow clap end
                     plot([sc_evs(seqrows(onev),1) sc_evs(seqrows(onev),1)],[-1 s+1],'r','linewidth',2); % slow clap start
                     plot([sc_evs(seqrows(onev),2) sc_evs(seqrows(onev),2)],[-1 s+1],'b','linewidth',2); % slow clap end
                     
                 end
                 
             end
             
         end
         
         % Plot lines of known explosion on/off times
         
         hold all;
         
         if exist('exp_evs','var') == 1
             
             % Find sequences within current time period
             
             seqrows = find((exp_evs(:,1)>=tt & exp_evs(:,1)<=(tt+7)) | (exp_evs(:,2)>=tt & exp_evs(:,2)<=(tt+7))); % either start or stop in period
             
             if ~isempty(seqrows) % if there are some, plot them
                 
                 for onev = 1:length(seqrows)
                     
                     %plot(([exp_evs(seqrows(onev),1) exp_evs(seqrows(onev),1)] - tstrt)*24,[0 s],'r'); % slow clap start
                     %plot(([exp_evs(seqrows(onev),2) exp_evs(seqrows(onev),2)] - tstrt)*24,[0 s],'b'); % slow clap end
                     plot([exp_evs(seqrows(onev),1) exp_evs(seqrows(onev),1)],[-1 s+1],':r','linewidth',2); % slow clap start
                     plot([exp_evs(seqrows(onev),2) exp_evs(seqrows(onev),2)],[-1 s+1],':b','linewidth',2); % slow clap end
                     
                 end
                 
             end
             
         end
         
         % Fix up plot labels, axes, etc.
         
         hold all;
         
         set(gca,'fontsize',14); % use larger font size
         box on;
         grid on;
         
         %xlabel(['Hours since ',datestr(tstrt)]);
         %xlim([0 ((tend-tstrt)*24)]); % plot over entire time period examined
         
         % Label x-axis with dates
         
         if tend < (tt+7) % if data ends before a week's time
         
             xlim([tt tend]); % plot over entire time period examined
         
         else
             
             xlim([tt (tt+7)]); % plot each week
             
         end
         
         datetick('x','mmm-dd HH:MM','keeplimits','keepticks'); % put date labels on x-axis
         
         if onspl == 1 % only put title on first subplot of figure
             
             title([{'Timeline of Sequences for All Stations'},{'\^ = on    o = off    <> = merger on    [] = merger off'},...
                 {'Black = level 1    Yellow = level 2'}]);
             
         end
         
         % Label y-axis with station names
         
         yticks(-1:(size(params.sta.list,1)+1)); % put an extra on both sides for box edges
         yticklabels(vertcat({''},{'Alert'},cellstr(params.sta.list),{''}));
         
         ylim([-1 (size(params.sta.list,1)+1)]);
         
         % If all subplots filled and there's more time to go, start new figure
         
         if onspl == splrows && tt < (tend - 7)
             
             figure;
             
             onspl = 1; % start over on first subplot
              
         elseif (tt + 7 + twin) > tstop % if there's not at least 1 time window of timeline remaining
             
             break; % stop plotting
             
         else % just move to next subplot
             
             onspl = onspl + 1;
             
         end
         
        end
        
    end
    
    
    % Plot each sequence separately
    
    if plseqs == 1
        
        for s = 1:1%size(params.sta.list,1)
            
            sta = strcat(params.sta.list(s,:)); % get station name without whitespaces
            
            cha = params.cha.list(s,:);
            
            % Find number of unique templates in seqdect
            
            seqs = unique(data.(sta).seqdect(:,1));
            
            figure;
            
            fignum = 1; % figure number
            
            prow = 5; % # of plots per row            
            pcol = 2; % # of plots per column
            pfig = prow * pcol; % # of subplots per figure
            
            onsp = 1;
            
            for snum = 1:length(seqs)
                                
                tname = strcat('temp',num2str(seqs(snum)));
                
                wfst = min(data.(sta).(tname).evdata(:,1)) - (seqT / 4); % seqT/4 before first event onset in seq
                
%                 wfst = min(data.(sta).seqdect(((data.(sta).seqdect(:,4)==1 | data.(sta).seqdect(:,4)==3) &...
%                     data.(sta).seqdect(:,1)==seqs(snum)),2)) - seqT; % seqT before first seq on
                
                wfend = max(data.(sta).seqdect(((data.(sta).seqdect(:,4)==0 | data.(sta).seqdect(:,4)==2) &...
                    data.(sta).seqdect(:,1)==seqs(snum)),2)) + (seqT / 4); % seqT/4 after last seq off
                
                % Get data from AVO Winston
                
                scnlList = scnlobject(sta,cha,params.net.list,'--');
                w = waveform(mySource,scnlList,wfst,wfend);
                
                if ~isempty(w) % only make plot if there's data
                    
                    % filter waveform
                    
                    f = filterobject('b', [lf hf], poles); % band-pass filter
                    filtw = filtfilt(f,w);
                    
                    wft = (get(filtw,'timevector') - wfst) * 24; % in hours
                    
                    % Plot waveform
                    
                    subplot(prow,pcol,onsp);
                    
                    plot(wft,get(filtw,'data'));
                    
                    % Plot on/off times as circles (on=open, off=filled, colored by ID #)
                    
                    % Find on/off times for plotting
                    
                    seqons = (data.(sta).seqdect((data.(sta).seqdect(:,1)==seqs(snum) &...
                        (data.(sta).seqdect(:,4)==1 | data.(sta).seqdect(:,4)==3)),2) - wfst) * 24; % in hours
                    seqoffs = (data.(sta).seqdect((data.(sta).seqdect(:,1)==seqs(snum)&...
                        (data.(sta).seqdect(:,4)==0 | data.(sta).seqdect(:,4)==2)),2) - wfst) * 24; % in hours
                    
                    amp = max(get(filtw,'data'));
                    
                    hold all;
                    
                    scatter(seqons,ones(size(seqons))*(3/4*amp),100,'r','filled','markeredgecolor','k'); % on times
                    scatter(seqoffs,ones(size(seqoffs))*(3/4*amp),100,'k'); % off times
                    
                    % Plot event onset times
                    
                    scatter(((data.(sta).(tname).evdata(:,1) - wfst) * 24),ones(size(data.(sta).(tname).evdata(:,1)))*(2/3*amp),36,...
                        'c','filled','markeredgecolor','k')
                    
                    xlabel(['Hours since ',datestr(wfst)]);
                    
                    title(['Template ',num2str(seqs(snum))]);
       
                    box on;
                    grid on;
                    
                end
                
                % If subplots filled and sequences still remain, start new figure
                
                if onsp == pfig && snum < length(seqs)
                    
                    figure;
                    
                    fignum = fignum + 1;
                    
                    onsp = 1;
                    
                else
                    
                    onsp = onsp + 1; % move to next subplot
                    
                end
                
            end
            
        end
        
    end
    
end


%% Function: Make Figure for Alert

function figname = make_alert_fig(params,directory,atype,atime,seqT,mySource)

figure('visible','off'); % make new figure window but don't show on screen

% Get waveforms

wfs = []; % start empty variable for waveforms

for s = 1:size(params.sta.list,1)
    
    sta = strcat(params.sta.list(s,:));
    
    scnlList = scnlobject(sta,params.cha.list(s,:),params.net.list,'--');
    w = waveform(mySource,scnlList,(atime-seqT),atime);
    
    if ~isempty(w) % if there's data, add to plot
        
        wfs = [wfs;w];
        
    end
    
end

% Make spectrogram plot

s = spectralobject(512,462,25,[20 120]);

specgram(s,wfs,'innerLabels',false);

% Save figure (and return filename)

figname = [directory,'alert_',atype,'_',datestr(atime,'yyyymmdd_HHMMSS'),'.png']; % make figure name (with directory)

print(figname,'-dpng'); % save figure

end


%% Function: Put input into useable vector lists

function listout = makelist(strlist,type,varargin)

%% Reformat list

if nargin > 2
    
    sortl = varargin{1};
    
else
    
    sortl = 'sorton';
   
end

if strcmpi(type,'num') == 1
    
    num1 = '!';
    num2 = '!';
    listout = '!';
    bound = 1;

    for a = 1:length(strlist)

        % if not slash or colon or list end, make array name

        if strcmpi('/',strlist(a)) == 0 && strcmpi(':',strlist(a)) == 0 && a ~= length(strlist)

            % see if temp number variables need to be reset

            if strcmpi('!',num1) == 1 || (strcmpi('!',num2) == 1 && strcmpi(':',strlist(a-1)) == 1)

                 vn = 1;

            else

                 vn = vn + 1;

            end

            if bound == 1

                num1(vn) = strlist(a);

            else

                num2(vn) = strlist(a);

            end

        % if colon, switch to find second number

        elseif strcmpi(':',strlist(a)) == 1

            bound = 2;

        % if slash or list end and not sequence end, add number to list

        elseif (strcmpi('/',strlist(a)) == 1 || a == length(strlist)) && bound == 1

            if a == length(strlist)

               if strcmpi(num1,'!') == 1

                    num1 = strlist(a);

                else

                    num1 = strcat(num1,strlist(a));

                end

            end

            if strcmpi(listout,'!') == 1

                listout = str2double(num1); % begin list

            else

                listout = vertcat(listout,str2double(num1)); % else add number to existing list

            end

            num1 = '!'; % reset temp name variable

        % if sequence end and slash or list end, make sequence and add to list

        elseif bound == 2 && (strcmpi('/',strlist(a)) == 1 || a == length(strlist))

            if a == length(strlist)

                if strcmpi(num2,'!') == 1

                    num2 = strlist(a);

                else

                    num2 = strcat(num2,strlist(a));

                end

            end

            templistout = transpose(str2double(num1):str2double(num2));

            if strcmpi(listout,'!') == 1

                listout = templistout; % begin list

            else

                listout = vertcat(listout,templistout); % else add number to existing list

            end

            % reset temp number variables and clear sequence counter

            num1 = '!'; 
            num2 = '!';
            bound = 1;

        end

    end
    
    % Sort list if wanted
    
    if strcmpi(sortl,'sortoff') == 0
        
        listout = sort(listout); % make sure final list is ordered
        
    end
    
elseif strcmpi(type,'str') == 1
    
    str1 = '!';
    listout = '!';
    
    for a = 1:length(strlist)

        % if not comma & not last char, make array name

        if strcmpi(',',strlist(a)) == 0 && a ~= length(strlist)

            if strcmpi('!',str1) == 1

                 n = 1;

            else

                 n = n + 1;

            end

            str1(n) = strlist(a);

            % if not comma & last char, make array name then add to list

        elseif strcmpi(',',strlist(a)) == 0 && a == length(strlist)

            if strcmpi('!',str1) == 1

                 n = 1;

            else

                 n = n + 1;

            end

            str1(n) = strlist(a);

            if strcmpi(listout,'!') == 1

                listout = str1; % begin list

            else

                listout = vertcat(listout,str1); % else add name to existing list

            end

            str1 = '!'; % reset temp name variable

        % if comma, add array name to list

        else

            if strcmpi(listout,'!') == 1

                listout = str1; % begin list

            else

                listout = vertcat(listout,str1); % else add name to existing list

            end

            str1 = '!'; % reset temp name variable

        end

    end
    
else
    
   disp(' ')
   disp('Error: The entered list type is invalid.')
    
end

end


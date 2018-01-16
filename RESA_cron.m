
% Real-time Detector for Repeating Event/Slow Clap Sequences (running as a cron job)
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
% Required functions:
%
% - Uses the ("old") GISMO Correlation Toolbox and Waveform Suite
%
%       By: Michael West and Celso Reyes (respectively)
%           Geophysical Institute, Alaska Volcano Observatory, U. Alaska Fairbanks
%
%       Website: https://geoscience-community-codes.github.io/GISMO/
%
% - Make sure the toolbox is installed. It is freely avaiable at the above website.
%
% - Requires updated version of @correlation/adjusttrig.m from GISMO suite to properly
%       align traces and updated version @correlation/private/xcorrrow.m
%
% - extractdatairis - only necessary if data will be retrieved via irisFetch

% --------------------

% By: Gabrielle Tepp, USGS AVO
% Created: 3/15/2017
% Last updated: 12/18/2017

%--------------------------------------------------------------------------%

warning off % don't display warnings


%% Read in parameters from file

pfid = fopen(paramfile);

while ~feof(pfid) % until end of file
    
    sline = fgetl(pfid); % get line
    
     % if it's not a blank line or comment (starting with % or #), evaluate it
     
    if ~isempty(sline) && strcmpi(sline(1),'%') == 0 && strcmpi(sline(1),'#') == 0
        
        while strcmpi(sline(end-2:end),'...') == 1 % if statement continued on next line, get it
            
            sline = [sline(1:end-3),fgetl(pfid)];
            
        end
        
        try
            
            eval(sline);
        
        catch
            
           disp(' ')
           disp('Bad line in parameter file. Can''t evaluate:')
           disp(sline)
            
        end
       
    end
    
end

% Put together sendto lists

sendtoon={};
sendtooff={};
sendtoonl2={};
sendtooffl2={};

for sr = 1:size(sendlist,1)
    
    if sendlist{sr,3} > 0
        sendtoon = [sendtoon,sendlist(sr,sendlist{sr,3})];
    end
    
    if sendlist{sr,4} > 0
        sendtooff = [sendtooff,sendlist(sr,sendlist{sr,4})];
    end
    
    if sendlist{sr,5} > 0
        sendtoonl2 = [sendtoonl2,sendlist(sr,sendlist{sr,5})];
    end
    
    if sendlist{sr,6} > 0
        sendtooffl2 = [sendtooffl2,sendlist(sr,sendlist{sr,6})];
    end
    
end


%% Get file information

% Put array, station, channel, and network info into useable lists

params.net.list = makelist(params.net.str,'str');
params.sta.list = makelist(params.sta.str,'str');
params.cha.list = makelist(params.cha.str,'str');

% Make sure there are station-channel pairs

if size(params.sta.list,1) ~= size(params.cha.list,1)
    
    disp(' ')
    disp('Number of stations and channels do not match. Will remove extra(s) from end of list.')
    
    if size(params.sta.list,1) > size(params.cha.list,1) % remove extra station(s)
        
        params.sta.list = params.sta.list(1:size(params.cha.list,1),:);
        
    else % remove extra channel(s)
        
        params.cha.list = params.cha.list(1:size(params.sta.list,1),:);
        
    end
    
    
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

% Make sure directory ends with "/"

if strcmpi(directory(length(directory)),'/') == 0
    
   directory = [directory '/']; 
    
end

% Load data, alert, and tID structures if saved file exists

if exist([directory,structfile],'file') == 2
    
    disp(' ')
    disp('Reading in templates and other data....')
    
    load([directory,structfile]);
    
else % if files don't exist, initiate alert and tID variables
    
    % Make list of available template ID #s
    
    for ss = 1:size(params.sta.list,1)
        
        sta = strcat(params.sta.list(ss,:));
        
        tID.(sta) = transpose(1:ceil(tempmax*1.5));
        
    end
    
    % Start alert status at 0
    
    alert.lv1.status = 0;
    
    if l2chk == 1
        
        alert.lv2.status = 0;
        
    end
    
end

ovlp = (l_lta + evwin) / 60; % set overlap of time windows; in min

if twin <= ovlp
    
   disp(' ')
   error('Error: Time window is shorter than overlap time. Can''t continue.') 
   
end

% Convert times to days

twin = twin / 60 / 24;
ovlp = ovlp / 60 / 24;
seqT = seqT / 60 / 24;
seqTl2 = seqTl2 / 60 / 24;
holdt = holdt / 60 / 24;
holdtl2 = holdtl2 / 60 / 24;
seqToff = seqToff /60 / 24;
seqToffl2 = seqToffl2 /60 / 24;

% Check that parameters are good

if xcend > evwin
    
    xcend = evwin;
    
end


%% Start Detector

% Set window start time

nowt = datevec(now_utc); % get current UTC time as a date vector

nowt(1,6) = 0; % set sec = 0 (i.e., floor to nearest minute)

t = datenum(nowt) - twin; % compute start of time window

tcheck = t;

% Initiate fields

if ~exist('data','var')
    
    for ss = 1:size(params.sta.list,1)
        
        sta = strcat(params.sta.list(ss,:));
        
        data.(sta).seqon = [];
        
        if l2chk == 1
            
            data.(sta).seqonl2 = []; % lv 2 on list
            
        end
        
    end
    
end

% Start new log file if it doesn't yet exist

if ~exist([directory logfile],'file')
    
    fID = fopen([directory logfile],'a'); % open log file for writing (append)
    
    fprintf(fID,'%s\n%s\n\n','      1st Event              Most Recent Event      Total # of Events    Template ID    On/Off    Station',...
        '----------------------    ----------------------    -----------------    -----------    ------    -------');
    
    fclose(fID);

end

    % Send check every 12 hours that code's still running
    
    if sendchk == 1
        
        [~,~,~,hh,mm,~] = datevec(t); % current time window hour & min
        
        if ((hh + mm/60 + ovlp/24) <= 5 && (hh + mm/60 + twin*24) > 5) ||...
                ((hh + mm/60 + ovlp/24) <= 17 && (hh + mm/60 + twin*24) > 17)
            
            sendmail(sendtochk,'Everything''s good with cron alarm!',['Still running at ' datestr(now) ' AKDT. Nothing to report.']);
            
            disp(['Check sent at ' datestr(t+twin) '.'])
            
        end
        
    end
    
    fID = fopen([directory logfile],'a'); % open log file for writing (append)
    
    for s = 1:size(params.sta.list,1)
        
        %% Get data
        
        sta = strcat(params.sta.list(s,:)); % get station name without any whitespaces
        
        if ~isempty(dbtype) && strcmpi(dbtype,'IRIS') == 1
         
            % Get data from IRIS
            
            tempw = irisFetch.Traces(params.net.list,sta,'--',params.cha.list(s,:),t,(t+twin));
            
            if isempty(tempw) == 1 % if no data retrieved
                
                data.(sta).wf = []; % make empty waveform
                
            else % put data into waveform object
                
                tempd = extractdatairis(tempw,tempw(1).sampleRate,t,(t+twin),NaN); % combine data "chunks"
                
                data.(sta).wf = waveform(sta,params.cha.list(s,:),tempw(1).sampleRate,t,tempd); % put into waveform object
                
            end
            
        else
            
            % Get data from AVO Winston
            
            scnlList = scnlobject(sta,params.cha.list(s,:),params.net.list,'--');
            data.(sta).wf = waveform(mySource,scnlList,t,(t+twin));
            
        end
    

        %% Run Detector if There's Data

        if ~isempty(data.(sta).wf)

            %% Filter

            data.(sta).wf = zero2nan(data.(sta).wf,10); % change zeros (min of 10 consecutive) to NaNs
            
            data.(sta).wf = demean(data.(sta).wf); % remove mean
            
            if ~isempty(find(isnan(get(data.(sta).wf,'data')),1)) % if there are NaNs, replace with 0 before filtering
                
                tempd = get(data.(sta).wf,'data');
                
                tempd(isnan(tempd)) = 0;
                
                data.(sta).wf = set(data.(sta).wf,'data',tempd);
                
                clear tempd; % clear temporary variable
                
            end
            
            % Make filter object

            f = filterobject('b', [lf hf], poles); % band-pass filter

            data.(sta).wf = filtfilt(f,data.(sta).wf);


            %% Do STA-LTA Filtering

            %         disp(' ')
            %         disp('Applying STA-LTA filter....')
            
            % Call STA-LTA function
            % [l_sta l_lta th_on th_off min_sep min_dur] (times in s)

            data.(sta).evcat = sta_lta(data.(sta).wf,'edp',[l_sta l_lta threson thresoff min_sep min_dur]);

            % if event(s) found, go to next step

            if isempty(data.(sta).evcat) == 0


                %% Check for Matching Waveforms

                for onev = 1:size(data.(sta).evcat,1)

                    wfst = data.(sta).evcat(onev,1) - (buff / 3600 / 24); % in days
                    wfend = wfst + (evwin / 3600 / 24); % in days
                    
                    % Make sure waveform fits in current data window
                    % Note that the STA/LTA won't pick up events within the first LTA window of
                    % the waveform, but that's enforced here anyway to avoid potential repeats
                    
                    if wfend <= get(data.(sta).wf,'end') && wfst >= get(data.(sta).wf,'start') &&...
                        wfst > (get(data.(sta).wf,'start')+(l_lta/3600/24))

                        evwf = extract(data.(sta).wf,'TIME',wfst,wfend);
                        
                        evwf = demean(evwf); % demean again just to make sure all event waveforms centered around zero

                        % if no templates currently exist, start list

                        if ~isfield(data.(sta),'templates')

                            data.(sta).templates = evwf;

                            data.(sta).tempdata = [1 1 0 0]; % [template_ID# total_#_matches sequence_declared lv2_declared]

                            data.(sta).temp1.evdata = data.(sta).evcat(onev,1); % keep event onset time
                            
                            tID.(sta)(1,:) = []; % remove # from available tID list

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

                                % Make corr object with just new waveform & matching template
                                
                                mcobj = correlation(waveform(cobj,[tnum length(lagmat)]),trigs([tnum length(lagmat)],1));
                                
                                % Adjust trigger to matching template

                                mcobj = xcorr(mcobj,[xcst xcend]); % need to repopulate lag matrix
                                
                                mcobj = adjusttrig(mcobj,'INDEX',1);

                                % Crop to matching template window

                                mcobj = crop(mcobj,0,evwin);

                                stackcobj = waveform(stack(mcobj)); % stack matching template and new event

                                nonmatches = find((1:length(data.(sta).templates) ~= tnum)); % list of non-matching template #s

                                % add event onset to template list

                                tname = strcat('temp',num2str(data.(sta).tempdata(tnum,1)));

                                % keep event onset time
                                
                                data.(sta).(tname).evdata = vertcat(data.(sta).(tname).evdata,data.(sta).evcat(onev,1));

                                %data.(sta).(tname).template = stackcobj(length(stackcobj)); % keep "new" stacked template

                                % update total match count

                                data.(sta).tempdata(tnum,2) = data.(sta).tempdata(tnum,2) + 1;

                                % update template and info

                                data.(sta).templates = vertcat(stackcobj(length(stackcobj)),data.(sta).templates(nonmatches));
                                data.(sta).tempdata = vertcat(data.(sta).tempdata(tnum,:),data.(sta).tempdata(nonmatches,:)); % re-order tempdata

                            % no match found and not a duplicate, add to list
                                 
                            elseif isempty(data.(sta).evcat(onev,1)==data.(sta).(['temp',num2str(data.(sta).tempdata(tnum,1))]).evdata(:,1)) == 0

                                if size(data.(sta).templates,1) < tempmax % if still below max allowed template #

                                    data.(sta).templates = vertcat(evwf,data.(sta).templates);

                                    data.(sta).tempdata = vertcat([tID.(sta)(1,1) 1 0 0],data.(sta).tempdata);

                                    tname = strcat('temp',num2str(tID.(sta)(1,1)));

                                    data.(sta).(tname).evdata = data.(sta).evcat(onev,1); % keep event onset time
                                    
                                    tID.(sta)(1,:) = []; % remove # from available tID list
                                    
                                else % Only keep (tempmax - 1) most recent templates + newest addition (and any in sequence)

                                    % Find last listed template without sequence "on"
                                    
                                    rmtnum = find((data.(sta).tempdata(:,3) == 0),1,'last');
                                    
                                    % Remove field entry of disappearing template
                                    
                                    tname = strcat('temp',num2str(data.(sta).tempdata(rmtnum,1)));
                                    data.(sta) = rmfield(data.(sta),tname);
                                    
                                    tID.(sta) = [tID.(sta);data.(sta).tempdata(rmtnum,1)]; % make tID available at end of list
                                                                 
                                    % Add new template to top of list and corresponding info
  
                                    keeptemp = find((1:length(data.(sta).templates)) ~= rmtnum); % make list of templates to keep
                                    
                                    data.(sta).templates = vertcat(evwf,data.(sta).templates(keeptemp,:));

                                    data.(sta).tempdata = vertcat([tID.(sta)(1,1) 1 0 0],data.(sta).tempdata(keeptemp,:));

                                    tname = strcat('temp',num2str(tID.(sta)(1,1)));

                                    data.(sta).(tname).evdata = data.(sta).evcat(onev,1); % keep event onset time
                                    
                                    tID.(sta)(1,:) = []; % remove # from available tID list

                                end

                            end

                        end

                    end

                end

            end


            %% Consolidate Highly Correlated Templates

            if isfield(data.(sta),'templates') == 1 % if there are templates

                checktemp = 1;

                while checktemp == 1

                    % x-correlate current template events

                    cc = xcorr(correlation(data.(sta).templates,get(data.(sta).templates,'start')),[xcst xcend]);

                    ccmat = get(cc,'corr'); % get correlation matrix

                    % Find any cc values above the minimum

                    ccmat = ccmat - eye(length(ccmat)); % subtract off diagonal to avoid auto-correlation values

                    lagmat = get(cc,'lag'); % get lag matrix

                    [row,col] = find((ccmat >= mincc) & (abs(lagmat) <= maxlag));

                    % If any templates correlate above min cc, combine them

                    if isempty(col)

                        checktemp = 0; % no correlated templates, so stop check

                    else % Consolidate first listed template with any matching, then re-check if any remain

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

                            % Make list of templates to stack ensuring lag time is acceptable

                            stackt = ucol(1,1); % start vector of template index #s

                            mergedata = data.(sta).(ftname).evdata; % start temp list of event data
                            seqstat = data.(sta).tempdata(ucol(1,1),3); % status for new template; equals 1 if any seq to stack are on
                            seqstatl2 = data.(sta).tempdata(ucol(1,1),4); % status for new template; equals 1 if any seq to stack are lv 2 on
                            
                            for m = 1:length(matchtemp)

                                % Check absolute lag time
                                
                                if abs(lagmat(row(matchtemp(m),1),col(matchtemp(m),1))) <= maxlag
                                    
                                    stackt = [stackt;row(matchtemp(m),1)];

                                    tname = strcat('temp',num2str(data.(sta).tempdata(row(matchtemp(m),1),1)));

                                    % Combine event times list with first template (and sort)

                                    mergedata = sortrows(vertcat(mergedata,data.(sta).(tname).evdata));

                                    % if sequence is on, turn off and record
                                    
                                    if data.(sta).tempdata(row(matchtemp(m)),3) == 1 

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
                                            
                                            seqstatl2 = 1; % status for new template; equals 1 if any seq to stack are lv 2 on
                                            
                                            disp(['Template ' num2str(data.(sta).tempdata(row(matchtemp(m),1),1))...
                                            ' (' sta ') is currently level 2 but is being merged.'])
                                        
                                            % Write activity to log file
                                        
                                            fprintf(fID,'%s\n',['Template ' num2str(data.(sta).tempdata(row(matchtemp(m),1),1))...
                                            ' (' sta ') is currently level 2 but is being merged.']);
                                        
                                        end

                                    end
                                    
                                    % remove event times
                                    
                                    data.(sta) = rmfield(data.(sta),tname);
                                    
                                    tID.(sta) = [tID.(sta);data.(sta).tempdata(row(matchtemp(m),1),1)]; % make tID available at end of list

                                end

                            end

                            % If at least one template match meets criteria, consolidate

                            if length(stackt) > 1

                                ntname = strcat('temp',num2str(tID.(sta)(1,1)));

                                % if sequence is on, turn off and record
                                
                                if data.(sta).tempdata(ucol(1,1),3) == 1 

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
                                         
                                        disp(['Template ' num2str(data.(sta).tempdata(ucol(1,1),1))...
                                            ' (' sta ') is currently level 2 but is being merged.'])
                                        
                                        % Write activity to log file
                                        
                                        fprintf(fID,'%s\n',['Template ' num2str(data.(sta).tempdata(ucol(1,1),1))...
                                            ' (' sta ') is currently level 2 but is being merged.']);
                                    
                                    end

                                end
                                
                                % remove event times of first matching template
                                
                                data.(sta) = rmfield(data.(sta),ftname);

                                tID.(sta) = [tID.(sta);data.(sta).tempdata(ucol(1,1),1)]; % make tID available at end of list

                                % Make entry for new template

                                data.(sta).(ntname).evdata = mergedata;
%                                 data.(sta).(ntname).template = stackcobj(length(stackcobj));

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
                                
                                adjcc = crop(adjcc,0,evwin); % Crop to window of first (most recent) template

                                stackcobj = waveform(stack(adjcc,stackt)); % stack correlated templates
                                
                                data.(sta).templates = vertcat(stackcobj(length(stackcobj)),data.(sta).templates(nonmatch));

                                data.(sta).tempdata = vertcat([tID.(sta)(1,1) size(data.(sta).(ntname).evdata,1) seqstat seqstatl2],...
                                    data.(sta).tempdata(nonmatch,:));
                                
                                tID.(sta)(1,:) = []; % remove # from available list

                                % Add entry to sequence timeline and sequence on lists if it's turned on

                                if seqstat == 1

                                    data.(sta).seqon = [data.(sta).seqon; data.(sta).tempdata(1,1)];

                                    disp(['Merged template ' num2str(data.(sta).tempdata(1,1)) ' (' sta ') is being turned on.'])
                                    
                                    % Write activity to log file
                                        
                                    fprintf(fID,'%s\n',['Merged template ' num2str(data.(sta).tempdata(1,1))...
                                        ' (' sta ') is being turned on.']);

                                end

                                % Add entry to sequence timeline and sequence on lists if it's turned on lv 2
                                
                                if l2chk == 1 && seqstatl2 == 1

                                    data.(sta).seqonl2 = [data.(sta).seqonl2; data.(sta).tempdata(1,1)];

                                    disp(['Merged template ' num2str(data.(sta).tempdata(1,1)) ' (' sta ') is being turned on for level 2.'])

                                    % Write activity to log file
                                        
                                    fprintf(fID,'%s\n',['Merged template ' num2str(data.(sta).tempdata(1,1))...
                                        ' (' sta ') is being turned on for level 2.']);
                                    
                                end
                                
                            end

                        end

                        % clear temporary variables

                        clear adjcc cc stackobj ccmat row col lagmat m ucol stackt ftname tname consoltimes ind nonmatch seqstat;

                    end

                end

            end
            
        else
            
            disp(['No data found for station ' sta])
            
        end
        
        
        %% Check for Event Sequences
        
        % if there is template data (i.e. at least one event)
        
        if isfield(data.(sta),'tempdata')
            
            evseq = find(data.(sta).tempdata(:,2)>=minev); % look for event sequences
            
            if isempty(evseq) == 0 % if at least one template has enough events for sequence

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
    
    lv1staon = []; % reset lv 1 on list
    lv2staon = []; % reset lv 2 on list
    
    for ss = 1:size(params.sta.list,1)
        
        sta = strcat(params.sta.list(ss,:)); % get station name without whitespaces
        
        if ~isempty(data.(sta).seqon)
            
            onct = onct + 1; % increase count by one
            
            lv1staon = [lv1staon,sta,', ']; % make a list of "on" stations
            
            % Check to see if it's on preferred list (if there is one)
            
            if ~isempty(reqsta)
                
                for rr = 1:size(reqsta,1)
                    
                    if strcmpi(reqsta(rr,:),sta) == 1 % if it is, count it
                        
                        reqct = reqct + 1;
                        
                        break; % leave loop
                        
                    end
                    
                end
                
            end
            
        end
        
        if l2chk == 1 && ~isempty(data.(sta).seqonl2) % check station for sequence(s) on level 2
            
            onctl2 = onctl2 + 1;
            
            lv2staon = [lv2staon,sta,', ']; % make a list of "on" stations
            
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
    
    if onct >= minsta && reqct >= minreq && alert.lv1.status == 0
        
        disp(' ')
        disp(['** ' num2str(onct) ' stations have at least one sequence in progress!'])
        disp(['Current time: ' datestr(t+twin)])
        
        % Write to log file
        
        fprintf(fID,'%s\n',['** ' num2str(onct) ' stations (' num2str(reqct)...
            ' required) have at least one sequence in progress at ' datestr(t+twin) '!']);
        
        alert.lv1.status = 1; % change alert status to on
        
        alert.lv1.ontime = t+twin; % keep current alert on time
        
        alert.lv1.curon = t + twin; % latest time alert is on
        
        % Make spectrogram for alert
        
        if strcmpi(dbtype,'IRIS') == 1
            
            figname = make_alert_fig(params,directory,'ON',(t+twin),seqT,dbtype);
            
        else
            
            figname = make_alert_fig(params,directory,'ON',(t+twin),seqT,dbtype,mySource);
            
        end
        
        % Send email/text notification
        
        try
            
            sendmail(sendtoon,'Event Sequence in Progress!',...
                [num2str(onct) ' stations have at least one sequence in progress! Time is ' datestr(t+twin) ' UTC. '...
                10 10 'Sequence parameters: min ' num2str(minev) ' events in ' num2str(seqT*24*60) ' minutes.' 10 10 ...
                'Stations on: ' lv1staon(1:end-2)],figname);
            
            disp('Notifications successfully sent.')
            
        catch %if it fails to send, try each recipient separately
            
            for l = 1:size(sendtoon,2)
                
                try
                    
                    sendmail(sendtoon{1,l},'Event Sequence in Progress!',...
                        [num2str(onct) ' stations have at least one sequence in progress! Time is ' datestr(t+twin) ' UTC. '...
                        10 10 'Sequence parameters: min ' num2str(minev) ' events in ' num2str(seqT*24*60) ' minutes.' 10 10 ...
                        'Stations on: ' lv1staon(1:end-2)],figname);
                    
                catch
                    
                    disp(['Error: Cannot send alert to ' sendtoon{1,l} '.'])
                    
                end
                
                
            end
            
        end
        
    elseif onct >= minsta && reqct >= minreq && alert.lv1.status == 1 % if alert still valid, keep track of time
        
        alert.lv1.curon = t + twin; % latest time alert is on
        
    elseif  (onct < minsta || reqct < minreq) && alert.lv1.status == 1 % else, send off alert if dropped below min (and on)
        
        % if alert has been "on" for required period of time and is past waiting time to turn off or if level 1 turns off
        
        if (alert.lv1.curon + (twin * offwait)) < (t + twin) && ((t+twin) - alert.lv1.ontime(end)) > holdt
            
            disp(' ')
            disp('** Stations with at least one sequence in progress no longer meet the minimum.')
            disp(['Current time: ' datestr(t+twin)])
            
            % Write to log file
            
            fprintf(fID,'%s\n',['** Stations with at least one sequence in progress no longer meet the minimum at ' datestr(t+twin) '.']);
            
            alert.lv1.status = 0; % change alert status to off
            
            alert.lv1.curon = []; % clear latest alert on time
            
            params.alert.lv1.ontime = []; % clear alert on time
            
            % Make figure for alert
            
            if strcmpi(dbtype,'IRIS') == 1
                
                figname = make_alert_fig(params,directory,'OFF',(t+twin),seqT,dbtype);
                
            else
                
                figname = make_alert_fig(params,directory,'OFF',(t+twin),seqT,dbtype,mySource);
                
            end
            
            % Send email/text notification
            
            try
                
                if strcmpi(figname,'none') == 1
                    
                    sendmail(sendtooff,'Sequence Has Ended',...
                        ['Stations with at least one sequence in progress no longer meet the minimum. Time is ' datestr(t+twin) ' UTC.'...
                        10 10 'Sequence parameters: min ' num2str(minev) ' events in ' num2str(seqT*24*60) ' minutes.' 10 10 ...
                        'Stations on: ' lv1staon(1:end-2) 10 10 'No data found.']);
                    
                else
                    
                    sendmail(sendtooff,'Sequence Has Ended',...
                        ['Stations with at least one sequence in progress no longer meet the minimum. Time is ' datestr(t+twin) ' UTC.'...
                        10 10 'Sequence parameters: min ' num2str(minev) ' events in ' num2str(seqT*24*60) ' minutes.' 10 10 ...
                        'Stations on: ' lv1staon(1:end-2)],figname);
                    
                end
                
            catch %if it fails to send, try each recipient separately
                
                for l = 1:size(sendtooff,2)
                    
                    try
                        
                        if strcmpi(figname,'none') == 1
                            
                            sendmail(sendtooff{1,l},'Sequence Has Ended',...
                                ['Stations with at least one sequence in progress no longer meet the minimum. Time is ' datestr(t+twin) ' UTC.'...
                                10 10 'Sequence parameters: min ' num2str(minev) ' events in ' num2str(seqT*24*60) ' minutes.' 10 10 ...
                                'Stations on: ' lv1staon(1:end-2) 10 10 'No data found.']);
                            
                        else
                            
                            sendmail(sendtooff{1,l},'Sequence Has Ended',...
                                ['Stations with at least one sequence in progress no longer meet the minimum. Time is ' datestr(t+twin) ' UTC.'...
                                10 10 'Sequence parameters: min ' num2str(minev) ' events in ' num2str(seqT*24*60) ' minutes.' 10 10 ...
                                'Stations on: ' lv1staon(1:end-2)],figname);
                            
                        end
                        
                    catch
                        
                        disp(['Error: Cannot send alert to ' sendtooff{1,l} '.'])
                        
                    end
                    
                    
                end
                
            end
            
        end
        
    end
    
    % If at least minsta stations have a sequence on level 2, send alert (if new level up)
    
    if l2chk == 1
        
        if onctl2 >= minsta && reqctl2 >= minreq && alert.lv2.status == 0 && alert.lv1.status == 1
            
            disp(' ')
            disp(['** ' num2str(onctl2) ' stations (' num2str(reqctl2) ' required) have at least one sequence at level 2!'])
            disp(['Current time: ' datestr(t+twin)])
            
            % Write to log file
            
            fprintf(fID,'%s\n',['** ' num2str(onctl2) ' stations (' num2str(reqctl2)...
                ' required) have at least one sequence at level 2 at ' datestr(t+twin) '!']);
            
            alert.lv2.status = 1;
            
            alert.lv2.ontime = t+twin; % keep current alert on time
            
            alert.lv2.curon = t + twin; % latest time alert is on
            
            % Make figure for alert
            
            if strcmpi(dbtype,'IRIS') == 1
                
                figname = make_alert_fig(params,directory,'lv2_ON',(t+twin),seqTl2,dbtype);
                
            else
                
                figname = make_alert_fig(params,directory,'lv2_ON',(t+twin),seqTl2,dbtype,mySource);
                
            end
            
            % Send email/text notification
            
            try
                
                sendmail(sendtoonl2,'Event Sequence at level 2!',...
                    [num2str(onctl2) ' stations have at least one sequence at level 2! Time is ' datestr(t+twin) ' UTC. '...
                    'Sequence parameters: min ' num2str(minevl2) ' events in ' num2str(seqTl2*24*60) ' minutes.' 10 10 ...
                    'Stations on level 2: ' lv2staon(1:end-2)],figname);
                sendmail(sendtoadd,'Event Sequence at level 2!',...
                    [num2str(onctl2) ' stations have at least one sequence at level 2! Time is ' datestr(t+twin) ' UTC. '...
                    'Sequence parameters: min ' num2str(minevl2) ' events in ' num2str(seqTl2*24*60) ' minutes.' 10 10 ...
                    'Stations on level 2: ' lv2staon(1:end-2)],figname);
            catch %if it fails to send, try each recipient separately
                
                for l = 1:size(sendtoonl2,2)
                    
                    try
                        
                        sendmail(sendtoonl2{1,l},'Event Sequence at level 2!',...
                            [num2str(onctl2) ' stations have at least one sequence at level 2! Time is ' datestr(t+twin) ' UTC. '...
                            'Sequence parameters: min ' num2str(minevl2) ' events in ' num2str(seqTl2*24*60) ' minutes.' 10 10 ...
                            'Stations on level 2: ' lv2staon(1:end-2)],figname);
                        
                    catch
                        
                        disp(['Error: Cannot send alert to ' sendtoonl2{1,l} '.'])
                        
                    end
                    
                    
                end
                
            end
            
        elseif onct >= minsta && reqctl2 >= minreq && alert.lv2.status == 1 && alert.lv1.status == 1 % if alert still valid, keep track of time
            
            alert.lv2.curon = t + twin; % latest time alert is on
            
            % else, send off alert if dropped below min (and on); also, force off if lv 1 alert went off
            
        elseif  (onctl2 < minsta || reqctl2 < minreq || alert.lv1.status == 0) && alert.lv2.status == 1
            
            % if alert has been "on" for required period of time and is past waiting time to turn off or if level 1 turns off
            
            if (alert.lv2.curon + (twin * offwaitl2)) < (t + twin) && ((t+twin) - alert.lv2.ontime(end)) > holdtl2 || alert.lv2.status == 0
                
                disp(' ')
                disp('** Stations with at least one sequence at level 2 no longer meet requirements.')
                disp(['Current time: ' datestr(t+twin)])
                
                % Write to log file
                
                fprintf(fID,'%s\n',['** Stations with at least one sequence at level 2 no longer meet requirements at ' datestr(t+twin) '.']);
                
                alert.lv2.status = 0;
                
                alert.lv2.ontime = []; % clear alert on time
                
                alert.lv2.curon = []; % clear latest alert on time
                
                % Make figure for alert
                
                if strcmpi(dbtype,'IRIS') == 1
                    
                    figname = make_alert_fig(params,directory,'lv2_OFF',(t+twin),seqTl2,dbtype);
                    
                else
                
                    figname = make_alert_fig(params,directory,'lv2_OFF',(t+twin),seqTl2,dbtype,mySource);
                
                end
                
                % Send email/text notification
                
                try
                    
                    if strcmpi(figname,'none') == 1
                        
                        sendmail(sendtooffl2,'Sequence No Longer at Level 2',...
                            ['Stations with at least one sequence at level 2 no longer meet the minimum. Time is ' datestr(t+twin) ' UTC.'...
                            10 10 'Sequence parameters: min ' num2str(minevl2) ' events in ' num2str(seqTl2*24*60) ' minutes.' 10 10 ...
                            'Stations on level 2: ' lv2staon(1:end-2) 10 10 'No data found.']);
                        
                    else
                        
                        sendmail(sendtooffl2,'Sequence No Longer at Level 2',...
                            ['Stations with at least one sequence at level 2 no longer meet the minimum. Time is ' datestr(t+twin) ' UTC.'...
                            10 10 'Sequence parameters: min ' num2str(minevl2) ' events in ' num2str(seqTl2*24*60) ' minutes.' 10 10 ...
                            'Stations on level 2: ' lv2staon(1:end-2)],figname);
                        
                    end
                    
                catch %if it fails to send, try each recipient separately
                    
                    for l = 1:size(sendtooffl2,2)
                        
                        try
                            
                            if strcmpi(figname,'none') == 1
                                
                                sendmail(sendtooffl2{1,l},'Sequence No Longer at Level 2',...
                                    ['Stations with at least one sequence at level 2 no longer meet the minimum. Time is ' datestr(t+twin) ' UTC.'...
                                    10 10 'Sequence parameters: min ' num2str(minevl2) ' events in ' num2str(seqTl2*24*60) ' minutes.' 10 10 ...
                                    'Stations on level 2: ' lv2staon(1:end-2) 10 10 'No data found.']);
                                
                            else
                                
                                sendmail(sendtooffl2{1,l},'Sequence No Longer at Level 2',...
                                    ['Stations with at least one sequence at level 2 no longer meet the minimum. Time is ' datestr(t+twin) ' UTC.'...
                                    10 10 'Sequence parameters: min ' num2str(minevl2) ' events in ' num2str(seqTl2*24*60) ' minutes.' 10 10 ...
                                    'Stations on level 2: ' lv2staon(1:end-2)],figname);
                                
                            end
                            
                        catch
                            
                            disp(['Error: Cannot send alert to ' sendtooffl2{1,l} '.'])
                            
                        end
                        
                        
                    end
                    
                end
                
            end
            
        end
        
    end
   
    
    %% Save Structures for Next Load
    
    % Clear time window waveforms to save space
    
    for ss = 1:size(params.sta.list,1)

        sta = strcat(params.sta.list(s,:)); % get station name without any whitespaces
        
        data.(sta).wf = [];
    
    end
    
    % Save structures
    
    save([directory,structfile],'data','alert','tID'); % will over-write any existing file
    
    fclose(fID); % close log file
    
%end


%% Function: Make Figure for Alert

function figname = make_alert_fig(params,directory,atype,atime,seqT,dbtype,varargin)

if nargin > 6 && strcmpi(dbtype,'IRIS') == 0 % if Winston data, mySource needs to be included as argument
    
    mySource = varargin{1};
    
end

figure('visible','off'); % make new figure window but don't show on screen

% Get waveforms

wfs = []; % start empty variable for waveforms

for s = 1:size(params.sta.list,1)
    
    sta = strcat(params.sta.list(s,:));
    
    if ~isempty(dbtype) && strcmpi(dbtype,'IRIS') == 1
        
        % Get data from IRIS
        
        tempw = irisFetch.Traces(params.net.list,sta,'--',params.cha.list(s,:),(atime-seqT),atime);
        
        if isempty(tempw) == 1 % if no data retrieved
            
            w = []; % make empty waveform
            
        else % put data into waveform object
            
            tempd = extractdatairis(tempw,tempw(1).sampleRate,(atime-seqT),atime,NaN); % combine data "chunks"
            
            w = waveform(sta,params.cha.list(s,:),tempw(1).sampleRate,(atime-seqT),tempd); % put into waveform object
            
        end
        
    else
        
        % Get data from AVO Winston
        
        scnlList = scnlobject(sta,params.cha.list(s,:),params.net.list,'--');
        w = waveform(mySource,scnlList,(atime-seqT),atime);
        
    end
    
    if ~isempty(w) % if there's data, add to plot
        
        wfs = [wfs;w];
        
    end
    
end

% Make spectrogram plot

if ~isempty(wfs)

s = spectralobject(512,462,25,[20 120]);

specgram(s,wfs,'fontsize',12);

set(gcf,'Units','pixels','OuterPosition',[50 50 800 1100]); % [left bottom width height]

% Save figure (and return filename)

figname = [directory,'alert_',atype,'_',datestr(atime,'yyyymmdd_HHMMSS'),'.png']; % make figure name (with directory)

print(figname,'-dpng'); % save figure

else
    
   figname = 'none'; 
    
end

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


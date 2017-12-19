function data = opencss(directory,varargin)

% opencss.m - Gets data from CSS files using info in wfdisc files.
%
% --------------------
%
% USAGE:
%
%   To use optional inputs: data = opencss(directory,'varname',value,...)
%
% INPUT:
%
%   'directory': directory containing wfdisc files (will prompt for if not specified)
%
%   varargin: series of parameters for waveforms (will use default values if not specified)
%
%       - 'openall': 1 to open and return all CSS files, 0 to open specified files only (default)
%
%       - 'files': column vector list of wfdisc file names if not opening all (will prompt for if not specified)
%
%       - 'station': column vector list of stations (default all)
% 
%       - 'stastr': comma separated list of stations (default all, all names same length)
% 
%       - 'channel': column vector list of stations (default all)
% 
%       - 'chastr': comma separated list of channels (default all, all names same length)
% 
%       * Note: only include one each of station/stastr and channel/chastr
% 
% OUTPUT:
% 
%   'data': structure containing data obtained from files
%          Format: data.file#.station.channel, contains structure output by readcss
% 
% --------------------
% 
% NOTES:
% 
%   Uses externally available readcss.m function
% 
%   Requires function makelist.m

% --------------------

% By: Gabrielle Tepp, USGS AVO
% Created: 12/27/16
% Last updated: 8/7/17

%--------------------------------------------------------------------------%

%% Get Inputs and Set Defaults if Needed

% Get optional inputs, if any

if nargin > 1
    
    % Assume varname, value pairs
    
    for v = 1:2:size(varargin,2)
    
        if strcmpi(varargin{v},'openall') == 1
            
            openall = varargin{v+1};
            
        elseif strcmpi(varargin{v},'files') == 1
            
            filelist = varargin{v+1};
            
        elseif strcmpi(varargin{v},'station') == 1
            
            stalist = varargin{v+1};
            
        elseif strcmpi(varargin{v},'stastr') == 1
            
            stastr = varargin{v+1};
           
        elseif strcmpi(varargin{v},'channel') == 1
            
            chalist = varargin{v+1};
            
        elseif strcmpi(varargin{v},'chastr') == 1
            
            chastr = varargin{v+1};
            
        end
        
    end
    
end

% Set defaults if needed

if ~exist('openall','var')
    
    openall = 0;
    
end

if ~exist('stalist','var') && ~exist('stastr','var')
    
    stalist = '*';
    
end

if ~exist('chalist','var') && ~exist('chastr','var')
    
    chalist = '*';
    
end

% Make useable lists of stations and channels (if lists not provided)

if exist('stastr','var')
    
    stalist = makelist(stastr,'str');
    
end

if exist('chastr','var')
    
    chalist = makelist(chastr,'str');
    
end

% Get list of files if not specified

if openall == 0 && ~exist('filelist','var')
    
    disp(' ')
    files = input('List the wfdisc file names (comma separated, no spaces): ','s');
    
    % Make useable list of files
    
    filelist = makelist(files,'str');

end

if openall == 1
    
    % Make list of all wfdisc files in specified directory
    
    filedata = dir([directory,'*.wfdisc']);
    
    filelist = '';
    
    for f = 1:size(filedata,1)
        
        filelist = [filelist;filedata(f).name];
        
    end
    
end


%% Get Data from Files

for onf = 1:size(filelist,1)
    
    disp(' ')
    disp(['Working on file ' num2str(onf) ' of ' num2str(size(filelist,1)) '....'])
    
    % Read Wfdisc File to Get Needed Info to Open Files

    fID = fopen([directory,filelist(onf,:)],'r'); % open file

    % This format doesn't work 100% for some reason, but it gets the needed info
    % Maybe formatting changed at some point?
    lformat = '%6c %8c %17.5f %9u %8u %8u %17.5f %8u %11.7f %16.6f %16.6f %6c %1c %2c %1c %64c %32c %10u %9u %18c\';
    
    % Loop through file until end is reached
    
    fileinfo = [];
    
    while ~feof(fID)
       
        line = textscan(fID,lformat);
        
        fileinfo = vertcat(fileinfo,line);
        
    end
    
    % If station/channel lists are wildcards, replace
    
    if strcmpi(stalist,'*') == 1
        
       stalist = char(deblank(unique(fileinfo(:,1))));
        
    end
    
    if strcmpi(chalist,'*') == 1
        
       chalist = char(deblank(unique(fileinfo(:,2))));
        
    end
    
    % Open Data File(s)

    % Make name for structure
    
    fname = ['file' num2str(onf)];
    
    % Loop through all specified stations and channels
    
    for c = 1:size(chalist,1)
       
        cha = chalist(c,:);
        
        for s = 1:size(stalist,1)
                        
            sta = stalist(s,:);
            
            row = []; % clear row variable
            
            % Find row with matching station and channel
            
            for onr = 1:size(fileinfo,1)
                
                row = find(strcmpi(deblank(fileinfo{onr,2}),cha) & strcmpi(deblank(fileinfo{onr,1}),sta));
                
                % When match found, give notice, get data, and leave loop
                
                if ~isempty(row)
                    
                    disp(['Data found for ' sta ' ' cha '.'])
                    
                    data.(fname).(sta).(cha) = readcss([directory,filelist(onf,:)],sta,cha,fileinfo{onr,3},fileinfo{onr,7});
                   
                    break;
                    
                elseif onr == size(fileinfo,1) && isempty(row) == 1
                    
                    disp(['No data found for ' sta ' ' cha '.'])
                    
                end
    
            end
            
        end
        
    end

    fclose(fID);
    
end    


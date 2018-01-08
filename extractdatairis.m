% a function to extract data from output of irisFetch

% By: Matt Haney, Alaska Volcano Observatory

function dtaf = extractdatairis(dta,sr1,tc1,tc2,nd1)
% 
% dta: the irisFetch structure
% sr1: sample rate
% tc1: start of time window in calendar format used by irisFetch
% tc2: end of time window in calendar format used by irisFetch
% nd1: the no data value or no data flag (e.g., -(2^31), NaN, etc)

% expected data length
%dtaf = zeros(sr1*round((datenum(tc2)-datenum(tc1))*24*60*60)+1,1);
%Note: the one above returns # of samples + 1 for non-gappy data
dtaf = zeros(sr1*round((datenum(tc2)-datenum(tc1))*24*60*60),1);

% the previous time (updated in loop), set to beginning of window initially
dprev = datenum(tc1);
% sample countr number
countr = 1;
% number of data fragments (bb)
[aa bb] = size(dta);
for ii=1:bb
    
    % how many samples in gap since previous?
    ngp = round((dta(ii).startTime-dprev)*24*60*60*sr1 - 1);
    
    % if a gap, fill it in
    if (ngp > 0)
        dtaf(countr:(countr+ngp-1)) = nd1;
        countr = countr + ngp;
    else
    end
    
    % now past gap, put in some real data
    stwin = countr;
    fnwin = stwin + dta(ii).sampleCount-1;
    dtaf(stwin:fnwin) = dta(ii).data;
    
    % update the new previous time and advance the counter by 1
    dprev = dta(ii).endTime;
    countr = fnwin+1;
end

% how many samples in gap at end?
% could also assume we did the above correct and add in the necessary 
% number of missing samples to achieve the expected data length

% if no data
if (bb == 0)
    
    %dtaf = nd1*ones(sr1*round((datenum(tc2)-datenum(tc1))*24*60*60)+1,1);
    dtaf = nd1*ones(sr1*round((datenum(tc2)-datenum(tc1))*24*60*60),1);
    
% if there was data
else

% dlast = datenum(tc2);
% ngp = round((dlast-dta(ii).endTime)*24*60*60*sr1 - 1);
% if (ngp > 0)
%     dtaf(countr:(countr+ngp-1)) = nd1;
%     countr = countr + ngp-1;
% else
% end

dtaf(countr:end) = nd1;

end


% tlenw = t2-t1;
% 
% % was some data returned?
% if (sum(size(tr2)) == 2)
%     % grab the text output info
%     qq = char(tr2.toString);
%     % what is the real start time %
%     tr2_t1 = qq((length(qq)-54):(length(qq)-32));
%     % what is the real end time %
%     tr2_t2 = qq((length(qq)-22):(length(qq)-0));
%     % get data and make it double
%     d2=tr2.buffer;
%     %tstvari = tr2.buffer;
%     d2 = double(d2');
%     % find value for no data and save it
%     nd2 = tr2.NO_DATA;
%     % the sampling rate
%     sr2 = tr2.getSamplingRate;
%     % initialize
%     d2n = zeros(1,tlenw*sr2);
%     
%     % how do the real start and end times compare to what was asked for?
%     %
%     % convert the real start time string to [YYYY MM DD HH MM SS]
%     t1v2 = [ str2num(tr2_t1(1:4)) str2num(tr2_t1(6:7)) str2num(tr2_t1(9:10)) ...
%         str2num(tr2_t1(12:13)) str2num(tr2_t1(15:16)) str2num(tr2_t1(18:23)) ];
%     % how does this compare to what we asked for?
%     t1d2 = t1 - cal2sec(t1v2);
%     % convert the real end time string to [YYYY MM DD HH MM SS]
%     t2v2 = [ str2num(tr2_t2(1:4)) str2num(tr2_t2(6:7)) str2num(tr2_t2(9:10)) ...
%         str2num(tr2_t2(12:13)) str2num(tr2_t2(15:16)) str2num(tr2_t2(18:23)) ];
%     % how does this compare to what we asked for?
%     t2d2 = cal2sec(t2v2)-t2;
%     %
%     % assume in the following the sampling rate is 100 Hz
%     % cut out the correct time or put no data into the time series
%     %
%     % if we got more than what we asked for, cut out what we asked for
%     if (t1d2 >= 0 & t2d2 >= 0)
%         d2n = d2(int32(t1d2*sr2+1):int32(length(d2)-t2d2*sr2));
%     % if we got less than what we asked for at the start, 
%     % put no data values there
%     elseif (t1d2 < 0 & t2d2 >= 0)
%         d2n = [ nd2*ones(1,abs(int32(t1d2*sr2))) d2(1:int32(length(d2)-t2d2*sr2)) ];
%     % if we got less than what we asked for at the end, 
%     % put no data values there
%     elseif (t1d2 >= 0 & t2d2 < 0) 
%         d2n = [ d2(int32(t1d2*sr2+1):length(d2)) nd2*ones(1,abs(int32(t2d2*sr2))) ];
%     % if we got less than what we asked for at both the start and end, 
%     % put no data values in both places
%     else
%         d2n = [ nd2*ones(1,abs(int32(t1d2*sr2))) d2(1:length(d2)) nd2*ones(1,abs(int32(t2d2*sr2))) ];
%     end
%     
% else
%     
%     % leave d2n set to zero if nothing returned, also treat the other
%     % variables
%     sr2 = sre;
%     nd2 = -2^31;
%     d2n = nd2*ones(1,tlenw*sr2);
%     %d2n = zeros(1,tlenw*sr2);
%     %nd2 = -2^31;
%     
% end
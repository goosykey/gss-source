%function chanalyze (filename, optionsflag)

% FLAGINPUT
%   1  basic pie only
%   2  message curve
%   4  char curve
%   8  media curve
%   16 time chatting
%   32 dissappearing

filename = 'cm.txt';
optionsflag = 63;

options = de2bi(optionsflag,6);

%% common data
chiper = 0.1;
factor = 7;

%% file ops

f = fopen (filename,'r','n','UTF-8');

formatSpec = '%s';

A = textscan(f,formatSpec,'delimiter','\n');

A = A{1,1};

%% remove security info

if strfind(A{1,1}, 'шифрованием')
    A(1,:) = [];
end

%% check for faulty EOLs

i = 1;
while i < numel(A)
    if isempty(A{i})
        badstring = true;
    elseif ~isstrprop(A{i}(1),'digit')
        badstring = true;
    elseif A{i}(3) ~= '.'
        badstring = true;
    elseif ~isstrprop(A{i}(4),'digit')
        badstring = true;
    else
        badstring = false;
    end
    
    if ~badstring
        lastgood = i;
        i = i+1;
    else
        A{lastgood} = [A{lastgood}, ' ', A{i}];
        A(i) = [];
    end
end



%% into struct

namearray = cell(0,1);

for i = 1:numel(A)
    %id
    msg(i).id = i;
    
    %datetime
    str = A{i};
    d = datetime(str(1:17),'InputFormat','dd.MM.yyyy'', ''HH:mm');
    msg(i).timestamp = datestr(d);
    msg(i).datetime = d-5/24;
    msg(i).realdatetime = d;
    
    %name
    j = 21;
    namehere = '';
    while A{i}(j) ~= ':'
        namehere = [namehere, A{i}(j)];
        j = j+1;
    end
    msg(i).name = namehere;
    
    %nameid
    nameexists = false;
    nameidhere = numel(namearray)+1;
    for k = 1:numel(namearray)
        if strcmp(namehere, namearray{k})
            nameexists = true;
            nameidhere = k;
        end
    end
    if ~nameexists
        namearray{end+1} = namehere;
    end
    msg(i).nameid = nameidhere;
    
    %text
    msg(i).text = A{i}(j+2:end);
    
    %ifmedia
    if strcmp(msg(i).text,'<Без медиафайлов>')
        msg(i).media = true;
        msg(i).text = '';
    else
        msg(i).media = false;
    end
    
    %charlength
    msg(i).chars = numel(msg(i).text);
    
end

namecount = numel(namearray);

if options(1)
    %% preplotting
    nameids = [msg.nameid]';
    medias = [msg.media]';
    charss = [msg.chars]';
    
    msgcount = zeros(namecount,1);
    chrcount = zeros(namecount,1);
    mdacount = zeros(namecount,1);
    
    for i = 1:namecount
        msgcount(i) = sum((nameids==i));
        chrcount(i) = sum(charss(nameids==i));
        mdacount(i) = sum(medias(nameids==i));
    end
    
    
    %% plotting - general
    
    figure; hold on;
    explode = ones(namecount,1);
    
    subplot(3,1,1);
    pie(msgcount,explode,cellstr(num2str(msgcount)));
    title(sprintf('Total messages : %g',sum(msgcount)));
    
    subplot(3,1,2);
    pie(chrcount,explode,cellstr(num2str(chrcount)));
    title(sprintf('Total characters : %g',sum(chrcount)));
    
    subplot(3,1,3);
    pie(mdacount,explode,cellstr(num2str(mdacount)));
    title(sprintf('Total media : %g',sum(mdacount)));
    
    legend(namearray);
    
end


if any(options(2:4))
%% preplot - intensity curves
    dt0 = dateshift(msg(1).datetime,'start','day');
    dtend = dateshift(msg(end).datetime,'end','day');
    days = (dt0:dtend)';
    ids = [msg.id]';
    
    msgday = zeros(numel(days),namecount);
    chrday = zeros(numel(days),namecount);
    mdaday = zeros(numel(days),namecount);
    
    for i = 1:numel(days)
        dd = days(i);
        idshere = ids(dateshift([msg.datetime],'start','day')' == dd);
        msghere = msg(idshere);
        nameidshere = [msghere.nameid]';
        for j = 1:namecount
            msgday(i,j) = sum(nameidshere==j);
            chrday(i,j) = sum([msghere(nameidshere==j).chars]);
            mdaday(i,j) = sum([msghere(nameidshere==j).media]);
        end
    end
    
    msgday(:,end+1) = sum(msgday,2);
    chrday(:,end+1) = sum(chrday,2);
    mdaday(:,end+1) = sum(mdaday,2);
    
end

if options(2)
%% msg intensity
    
    figure;
    
    [daysq, gagli2] = chipize (days, msgday, chiper, factor);
    
    % gagli1 = gaavg(msgday,factor);
    % daysq = days(1):chiper:days(end);
    % gagli2 = [];
    % for i = 1:namecount+1
    %     garg = pchip(1:numel(days),gagli1(:,i),1:chiper:numel(days))';
    %     gagli2 = [gagli2,garg];
    % end
    plot(daysq,gagli2(:,1:end-1),'--'); hold on; grid;
    plot(daysq,gagli2(:,end),'-r','linewidth',3);
    plot(days,msgday(:,end),':r');
    legend([namearray,'TOTAL','RAW']);
    title('Messages sent');
end

if options(3)
    %% chr intensity
    
    figure;
    
    [daysq, gagli2] = chipize (days, chrday, chiper, factor);
    
    % gagli1 = gaavg(chrday,factor);
    % daysq = days(1):chiper:days(end);
    % gagli2 = [];
    % for i = 1:namecount+1
    %     garg = pchip(1:numel(days),gagli1(:,i),1:chiper:numel(days))';
    %     gagli2 = [gagli2,garg];
    % end
    plot(daysq,gagli2(:,1:end-1),'--'); hold on; grid;
    plot(daysq,gagli2(:,end),'-r','linewidth',3);
    plot(days,chrday(:,end),':r');
    legend([namearray,'TOTAL','RAW']);
    title('Characters sent');
end

if options(4)
    %% media intensity
    
    figure;
    
    [daysq, gagli2] = chipize (days, mdaday, chiper, factor);
    
    % gagli1 = gaavg(mdaday,factor);
    % daysq = days(1):chiper:days(end);
    % gagli2 = [];
    % for i = 1:namecount+1
    %     garg = pchip(1:numel(days),gagli1(:,i),1:chiper:numel(days))';
    %     gagli2 = [gagli2,garg];
    % end
    plot(daysq,gagli2(:,1:end-1),'--'); hold on; grid;
    plot(daysq,gagli2(:,end),'-r','linewidth',3);
    plot(days,mdaday(:,end),':r');
    legend([namearray,'TOTAL','RAW']);
    title('Media sent');
end


%% time in chat

if options(5)
    
    WOCHAT = zeros(numel(days),1);
    chattime = 5; %minutes
    
    for i = 1:numel(days)
        wochat = ones(1,1440);
        dd = days(i);
        idshere = ids(dateshift([msg.datetime],'start','day')' == dd);
        msghere = msg(idshere);
        mi = hour([msghere.datetime])*60 + minute([msghere.datetime]);
        for j = 1:numel(msghere)
            if mi(j) <= chattime
                mi(j) = chattime+1;
            elseif mi(j) > 1440-chattime
                mi(j) = 1440-chattime;
            end
            wochat((mi(j)-chattime):(mi(j)+chattime)) = 0;
        end
        
        WOCHAT(i) = sum(wochat,2)/60;
    end
    
    figure;
    gagli1 = gaavg(WOCHAT,factor);
    daysq = days(1):chiper:days(end);
    gagli2 = pchip(1:numel(days),gagli1,1:chiper:numel(days))';
    plot(daysq,24-gagli2,'-m','linewidth',2); hold on; grid;
    title('Time chatting');
    
end

if options(6)
%% DISSAPPEARING

dt0 = dateshift(msg(1).realdatetime,'start','day');
dtend = dateshift(msg(end).realdatetime,'end','day');
realdays = (dt0:dtend)';

responsetime = zeros(numel(realdays),namecount);
omitday = zeros(numel(realdays),1);

for i = 1:numel(realdays)
    dd = realdays(i);
    idshere = ids(dateshift([msg.realdatetime],'start','day')' == dd);
    msghere = msg(idshere);
    nameidshere = namecount - [msghere(1:end-1).nameid] + 1;
    mi = hour([msghere.realdatetime])*60 + minute([msghere.realdatetime]);
    diffmi = diff(mi);
    
    if isempty(diffmi)
        responsetime(i,:) = responsetime(i-1,:);
        omitday(i) = 1;
        continue
    end
    
    [maxwait,maxi] = max(diffmi); % adjust to not include sleepy time
    if maxwait > 210 && mi(maxi) < 360
        diffmi(maxi) = 0;
    end
    
    for j = 1:namecount
        responsetime(i,j) = mean(diffmi(nameidshere==j));
        if isnan(responsetime(i,j))
            responsetime(i,j) = responsetime(i-1,j);
            omitday(i) = 1;
        end
    end
    
    if i == 226
        testmsghere = msghere;
        testmi = mi;
        testdiffmi = diffmi;
        testmaxwait = maxwait;
        testmaxi = maxi;
    end
    
end

figure;
[daysq, gagli2] = chipize (realdays, responsetime, chiper, factor);
omitday1 = interp1(days,omitday,daysq);
gagli2(omitday1 > 0.5,:) = nan;

plot(daysq,gagli2,'linewidth',2); hold on; grid;
title('Mean time taken for response');
legend(namearray);

end

%end


%% functions


function colmod = gaavg (cols, factor)

n = numel(cols(:,1));
colmod = cols;

for i = 1:n
    
    rowi = cols(i,:)*0;
    
    if i < (factor+1)/2
        f = 2*i-1;
    elseif i > n-(factor-1)/2
        f = 2*(n-i)+1;
    else
        f = factor;
    end
    
    mtotal = 0;
    for j = (i - (f-1)/2):(i + (f-1)/2)
        m = (f+1)/2 - abs(i-j); %[i,j]
        rowj = cols(j,:);
        rowi = rowi + m*rowj;
        mtotal = mtotal + m;
    end
    
    rowi = rowi./mtotal;
    colmod(i,:) = rowi;
    
end


end

function [daysq, gagli2] = chipize (days, data, chiper, factor)

gagli1 = gaavg(data,factor);
daysq = days(1):chiper:days(end);
gagli2 = [];

for i = 1:numel(data(1,:))
    garg = pchip(1:numel(days),gagli1(:,i),1:chiper:numel(days))';
    gagli2 = [gagli2,garg];
end

end


%% UNUSED

function colmod = gaga1 (cols,maxgrowth)

F = 2*maxgrowth/pi; 
A = maxgrowth/2;
B = maxgrowth/6; %max deltadelta
colmod = cols;
delta = colmod(2,:) - colmod(1,:);
delta(abs(delta)>A) = F*atan(delta(abs(delta)>A)/A);
colmod(2,:) = colmod(1,:) + delta;

for i = 3:numel(cols(:,1))
    deltanew = colmod(i,:) - colmod(i-1,:);
    deltanew(abs(deltanew)>A) = F*atan(deltanew(abs(deltanew)>A)/A);
    deltadelta = deltanew-delta;
    deltadelta(abs(deltadelta) > B) = B*sign(deltadelta(abs(deltadelta) > B));
    deltanew = delta + deltadelta;
    colmod(i,:) = colmod(i-1,:) + deltanew;
    delta = deltanew;
end

colmod(colmod<0) = 0;

end







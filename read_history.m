%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read History 

% Needs the following variables
%   files:       List of files with the history

% Nabil 29/9/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clear files;


files{1} = 'output/hist_0128x0128x0153.dat';
% files{2} = '../post_coll_16x32_24pts_h18/statistics/hist_0512x0256x0153.1.dat';
% files{3} = '../post_coll_16x32_24pts_h18/statistics/hist_0512x0256x0153.2.dat';
% files{4} = '../post_coll_16x32_24pts_h18/statistics/hist_0512x0256x0153.3.dat';
% files{5} = '../post_coll_16x32_24pts_h18/statistics/hist_0512x0256x0153.4.dat';
% files{6} = '../post_coll_16x32_24pts_h18/statistics/hist_0512x0256x0153.5.dat';
% files{7} = '../post_coll_16x32_24pts_h18/statistics/hist_0512x0256x0153.6.dat';
% files{8} = '../post_coll_16x32_24pts_h18/statistics/hist_0512x0256x0153.7.dat';
% files{9} = '../post_coll_16x32_24pts_h18/statistics/hist_0512x0256x0153.8.dat';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = 0;
mpgx_index_correction = 0;
for j = 1:size(files,2)
  fname  = files{j};

  fid    = fopen(fname,'r');
           fread(fid,1,'int32');
  t_     = fread(fid,1,'float64');
  Re     = fread(fid,1,'float64');
  alp    = fread(fid,1,'float64');
  bet    = fread(fid,1,'float64');
  Lx     = 2*pi/alp;
  Lz     = 2*pi/bet;
  Ly     = 2;
  mpgx_  = fread(fid,1,'float64');
  nband  = fread(fid,1,'int32');
           fread(fid,2,'int32');
  Ngal   = fread(fid,[4 nband+2],'int32');
           fread(fid,1,'int32');
  
                 fread(fid,1,'int32');
  flag_ctpress = fread(fid,1,'int32');
  
  mpgx_index_correction(j) = i+1;   % Eliminate the first element of mpgx of every file

  while(flag_ctpress==0)
    i = i+1;
    iter(i)      = fread(fid,1,'int32');
    t   (i)      = fread(fid,1,'float64');
    dtv (i)      = fread(fid,1,'float64');
    dtc (i)      = fread(fid,1,'float64');
    dt  (i)      = fread(fid,1,'float64');
    err (i)      = fread(fid,1,'float64');
    Qx  (i)      = fread(fid,1,'float64');
    QxT (i)      = fread(fid,1,'float64');
    mpgx(i)      = fread(fid,1,'float64'); 
    dgx (i)      = fread(fid,1,'float64');
    Umax(i)      = fread(fid,1,'float64');
    Uslp(i)      = fread(fid,1,'float64');
                   fread(fid,1,'int32');
                   fread(fid,1,'int32');
    flag_ctpress = fread(fid,1,'int32');
  end
  
  while(flag_ctpress==1)
    i = i+1;
    iter (i)     = fread(fid,1,'int32');
    t    (i)     = fread(fid,1,'float64');
    dtv  (i)     = fread(fid,1,'float64');
    dtc  (i)     = fread(fid,1,'float64');
    dt   (i)     = fread(fid,1,'float64');
    err  (i)     = fread(fid,1,'float64');
    Qx   (i)     = fread(fid,1,'float64');
    mpgx (i)     = fread(fid,1,'float64');
    Umax (i)     = fread(fid,1,'float64');
    Uslp(i)      = fread(fid,1,'float64');
                   fread(fid,1,'int32');
                   fread(fid,1,'int32');
    flag_ctpress = fread(fid,1,'int32');
  end
end

t1                           = t;
t1(mpgx_index_correction)    = [];
mpgx1                        = mpgx;
mpgx1(mpgx_index_correction) = [];

startskip = 0;
disp(['skipping first ', num2str(startskip) ,'........'])
startskip = startskip + 1;

%% Plotting
  % Pretty colours
    %colours;
    colors = [[0     ,0.4470,0.7410];...
          [0.8500,0.3250,0.0980];...
          [0.9290,0.6940,0.1250];...
          [0.4940,0.1840,0.5560];...
          [0.4660,0.6740,0.1880];...
          [0.3010,0.7450,0.9330];...
          [0.6350,0.0780,0.1840]];
  
  % Setting up the figure
    figure(1), clf, hold on;
    xlabel('$t$','Interpreter','LaTex')
    set(gca,'fontsize',16);
    set(gca,'XScale','lin');
    box on;
  % Plots
    % Error 
    subplot(4,1,1), hold on;
    plot(t(startskip:end) ,err(startskip:end)                          ,'LineWidth',2,'Color',colors(5,:));
    plot(t(startskip:end) ,err(startskip:end)*0 + mean(err(startskip:end))    ,'--','LineWidth',2,'Color',colors(1,:));
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    ylabel('err','Interpreter','LaTex')
    % Mass flow
    subplot(4,1,2), hold on;
    plot(t(startskip:end) ,Qx(startskip:end)                           ,'LineWidth',2,'Color',colors(5,:));
    plot(t(startskip:end) ,Qx(startskip:end)*0 + mean(Qx(startskip:end))      ,'--','LineWidth',2,'Color',colors(1,:));
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    ylabel('Q','Interpreter','LaTex')
    % Mean pressure gradient
%    subplot(4,1,3), hold on;
%    plot(t1,mpgx1                        ,'LineWidth',2,'Color',colors(5,:));
%    plot(t1,mpgx1*0 + mean(mpgx1(:)),'--','LineWidth',2,'Color',colors(1,:));
%    plot(t1(startskip:end),mpgx1(startskip:end)*0 + mean(mpgx1(startskip:end)),'--','LineWidth',2,'Color',colors(2,:));
%    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
%    ylabel('$dp/dx$','Interpreter','LaTex')
    % Maximum velocity
    subplot(4,1,4), hold on;
    plot(t(startskip:end) ,Umax(startskip:end)                         ,'LineWidth',2,'Color',colors(5,:));
    plot(t(startskip:end) ,Umax(startskip:end)*0 + mean(Umax(startskip:end))  ,'--','LineWidth',2,'Color',colors(1,:));
    ylabel('$U_{max}$','Interpreter','LaTex')
    % Slip velocity
    subplot(4,1,3), hold on;
    plot(t(startskip:end) ,Uslp(startskip:end)                         ,'LineWidth',2,'Color',colors(5,:));
    plot(t(startskip:end) ,Uslp(startskip:end)*0 + mean(Uslp(startskip:end))  ,'--','LineWidth',2,'Color',colors(1,:));
    ylabel('$U_{slip}$','Interpreter','LaTex')


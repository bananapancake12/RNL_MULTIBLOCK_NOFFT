%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;

time = 05901;
%xx = 0128;
xx = 0128;
zz = 0128;
%zz = 0128;
%xx = 0128;
yy = 0153;
% xx = 0064;
% zz = 0064;
 %yy = 0033;

% yy = 0065;
%file = 'u2_0128x0256x0154_t00000.dat';

file = sprintf('output/u1_%4.4dx%4.4dx%4.4d_t%5.5d.dat',xx,zz,yy,time);
file

%% Reading
fid = fopen(file,'r');
% a = fread(fid,[1,1],'int');
h       = fread(fid,1,'*uint32');
t       = fread(fid,1,'real*8');
Re      = fread(fid,1,'real*8');
alp     = fread(fid,1,'real*8');
bet     = fread(fid,1,'real*8');
mpgx    = fread(fid,1,'real*8');
nband   = fread(fid,1,'int');
iter    = fread(fid,1,'int');
dummit  = fread(fid,90,'int');
N       = fread(fid,4*5,'int');
N       = reshape(N,4,5);
    fread(fid,1,'real*8'); %skip % fseek(fid,8,'cof');
y       = fread(fid,N(4,nband+1)-N(4,0+1)+2,'real*8');
dthetai = fread(fid,1,'real*8');
    fread(fid,1,'real*8'); %skip
dthdy   = fread(fid,N(4,nband+1)-N(4,0+1)+2,'real*8');


     %fread(fid,1,'int'); %skip
for index=1:N(4,1+1)-N(4,0+1)+1
    j1(index)        = fread(fid,1,'int');
    one1(index)      = fread(fid,1,'int');
    nx1(index)       = fread(fid,1,'int');
    nz1(index)       = fread(fid,1,'int');
    yj1(index)       = fread(fid,1,'real*8');
    u11(index,:)     = fread(fid,nx1(index)*nz1(index),'real*8');
                      fread(fid,1,'real*8');
end
for index=1:N(4,2+1)-N(4,1+1)
    j2(index)        = fread(fid,1,'int');
    one2(index)      = fread(fid,1,'int');
    nx2(index)       = fread(fid,1,'int');
    nz2(index)       = fread(fid,1,'int');
    yj2(index)       = fread(fid,1,'real*8');
    u12(index,:)     = fread(fid,nx2(index)*nz2(index),'real*8');
                      fread(fid,1,'real*8');
end
for index=1:N(4,3+1)-N(4,2+1)+1
    j3(index)        = fread(fid,1,'int');
    one3(index)      = fread(fid,1,'int');
    nx3(index)       = fread(fid,1,'int');
    nz3(index)       = fread(fid,1,'int');
    yj3(index)       = fread(fid,1,'real*8');
    u13(index,:)     = fread(fid,nx3(index)*nz3(index),'real*8');
                      fread(fid,1,'real*8');
end
fclose(fid);

u1 = [u11(:,1);u12(:,1);u13(:,1)];

%% Display
disp(['t=',num2str(t),' iter=',num2str(iter)]);
disp(['Re=',num2str(Re),' alp=',num2str(alp),' bet=',num2str(bet),...
    ' mpgx=',num2str(mpgx)]);
disp(['nband=',num2str(nband)]);
disp(N);
disp(['y =[',num2str(y(1)),',',num2str(y(end)),']']);
disp(['Umax=',num2str(max(u1))]);

%% Plot
colors = [[0     ,0.4470,0.7410];...
          [0.8500,0.3250,0.0980];...
          [0.9290,0.6940,0.1250];...
          [0.4940,0.1840,0.5560];...
          [0.4660,0.6740,0.1880];...
          [0.3010,0.7450,0.9330];...
          [0.6350,0.0780,0.1840]];

figure(5);

%clf
%hold on
%axis equal;
%plot(u1,-y,'linewidth',2,'color',colors(5,:));
hold on
plot(u1,y,'linewidth',2,'color',colors(1,:));
title('u1')
file = sprintf('output/u3_%4.4dx%4.4dx%4.4d_t%5.5d.dat',xx,zz,yy,time);
ylim([-1 1])

%% Reading
fid = fopen(file,'r');
% a = fread(fid,[1,1],'int');
h       = fread(fid,1,'*uint32');
t       = fread(fid,1,'real*8');
Re      = fread(fid,1,'real*8');
alp     = fread(fid,1,'real*8');
bet     = fread(fid,1,'real*8');
mpgx    = fread(fid,1,'real*8');
nband   = fread(fid,1,'int');
iter    = fread(fid,1,'int');
dummit  = fread(fid,90,'int');
N       = fread(fid,4*5,'int');
N       = reshape(N,4,5);
    fread(fid,1,'real*8'); %skip % fseek(fid,8,'cof');
y       = fread(fid,N(4,nband+1)-N(4,0+1)+2,'real*8');
dthetai = fread(fid,1,'real*8');
    fread(fid,1,'real*8'); %skip
dthdy   = fread(fid,N(4,nband+1)-N(4,0+1)+2,'real*8');


%     fread(fid,1,'int'); %skip
for index=1:N(4,1+1)-N(4,0+1)+1
    j1(index)        = fread(fid,1,'int');
    one1(index)      = fread(fid,1,'int');
    nx1(index)       = fread(fid,1,'int');
    nz1(index)       = fread(fid,1,'int');
    yj1(index)       = fread(fid,1,'real*8');
    u31(index,:)     = fread(fid,nx1(index)*nz1(index),'real*8');
                      fread(fid,1,'real*8');
end
for index=1:N(4,2+1)-N(4,1+1)
    j2(index)        = fread(fid,1,'int');
    one2(index)      = fread(fid,1,'int');
    nx2(index)       = fread(fid,1,'int');
    nz2(index)       = fread(fid,1,'int');
    yj2(index)       = fread(fid,1,'real*8');
    u32(index,:)     = fread(fid,nx2(index)*nz2(index),'real*8');
                      fread(fid,1,'real*8');
end
for index=1:N(4,3+1)-N(4,2+1)+1
    j3(index)        = fread(fid,1,'int');
    one3(index)      = fread(fid,1,'int');
    nx3(index)       = fread(fid,1,'int');
    nz3(index)       = fread(fid,1,'int');
    yj3(index)       = fread(fid,1,'real*8');
    u33(index,:)     = fread(fid,nx3(index)*nz3(index),'real*8');
                      fread(fid,1,'real*8');
end
fclose(fid);

u3 = [u31(:,1);u32(:,1);u33(:,1)];

%% Display
disp(['t=',num2str(t),' iter=',num2str(iter)]);
disp(['Re=',num2str(Re),' alp=',num2str(alp),' bet=',num2str(bet),...
    ' mpgx=',num2str(mpgx)]);
disp(['nband=',num2str(nband)]);
disp(N);
disp(['y =[',num2str(y(1)),',',num2str(y(end)),']']);
disp(['Umax=',num2str(max(u3))]);

%% Plot
colors = [[0     ,0.4470,0.7410];...
          [0.8500,0.3250,0.0980];...
          [0.9290,0.6940,0.1250];...
          [0.4940,0.1840,0.5560];...
          [0.4660,0.6740,0.1880];...
          [0.3010,0.7450,0.9330];...
          [0.6350,0.0780,0.1840]];

figure(6);
%clf
%hold on
%axis equal;
plot(u3,y,'linewidth',2,'color',colors(1,:));
title('u3')
ylim([-1 1])

file = sprintf('output/u2_%4.4dx%4.4dx%4.4d_t%5.5d.dat',xx,zz,yy,time);


%% Reading
fid = fopen(file,'r');
% a = fread(fid,[1,1],'int');
h       = fread(fid,1,'*uint32');
t       = fread(fid,1,'real*8');
Re      = fread(fid,1,'real*8');
alp     = fread(fid,1,'real*8');
bet     = fread(fid,1,'real*8');
mpgx    = fread(fid,1,'real*8');
nband   = fread(fid,1,'int');
iter    = fread(fid,1,'int');
dummit  = fread(fid,90,'int');
N       = fread(fid,4*5,'int');
N       = reshape(N,4,5);
    fread(fid,1,'real*8'); %skip % fseek(fid,8,'cof');
y       = fread(fid,N(3,nband+1)-N(3,0+1)+2,'real*8');
dthetai = fread(fid,1,'real*8');
    fread(fid,1,'real*8'); %skip
dthdy   = fread(fid,N(3,nband+1)-N(3,0+1)+2,'real*8');


%     fread(fid,1,'int'); %skip
for index=1:N(3,1+1)-N(3,0+1)+1
    j1(index)        = fread(fid,1,'int');
    one1(index)      = fread(fid,1,'int');
    nx1(index)       = fread(fid,1,'int');
    nz1(index)       = fread(fid,1,'int');
    yj1(index)       = fread(fid,1,'real*8');
    u21(index,:)     = fread(fid,nx1(index)*nz1(index),'real*8');
                      fread(fid,1,'real*8');
end
for index=1:N(3,2+1)-N(3,1+1)
    j2(index)        = fread(fid,1,'int');
    one2(index)      = fread(fid,1,'int');
    nx2(index)       = fread(fid,1,'int');
    nz2(index)       = fread(fid,1,'int');
    yj2(index)       = fread(fid,1,'real*8');
    u22(index,:)     = fread(fid,nx2(index)*nz2(index),'real*8');
                      fread(fid,1,'real*8');
end
for index=1:N(3,3+1)-N(3,2+1)+1
    j3(index)        = fread(fid,1,'int');
    one3(index)      = fread(fid,1,'int');
    nx3(index)       = fread(fid,1,'int');
    nz3(index)       = fread(fid,1,'int');
    yj3(index)       = fread(fid,1,'real*8');
    u23(index,:)     = fread(fid,nx3(index)*nz3(index),'real*8');
                      fread(fid,1,'real*8');
end
fclose(fid);

u2 = [u21(:,1);u22(:,1);u23(:,1)];

%% Display
disp(['t=',num2str(t),' iter=',num2str(iter)]);
disp(['Re=',num2str(Re),' alp=',num2str(alp),' bet=',num2str(bet),...
    ' mpgx=',num2str(mpgx)]);
disp(['nband=',num2str(nband)]);
disp(N);
disp(['y =[',num2str(y(1)),',',num2str(y(end)),']']);
disp(['Umax=',num2str(max(u3))]);

%% Plot
colors = [[0     ,0.4470,0.7410];...
          [0.8500,0.3250,0.0980];...
          [0.9290,0.6940,0.1250];...
          [0.4940,0.1840,0.5560];...
          [0.4660,0.6740,0.1880];...
          [0.3010,0.7450,0.9330];...
          [0.6350,0.0780,0.1840]];

figure(8);
%clf
%hold on
%axis equal;
plot(u2,y,'linewidth',2,'color',colors(1,:));
title('u2')
ylim([-1 1])

file = sprintf('output/p_%4.4dx%4.4dx%4.4d_t%5.5d.dat',xx,zz,yy,time);


%% Reading
fid = fopen(file,'r');
% a = fread(fid,[1,1],'int');
h       = fread(fid,1,'*uint32');
t       = fread(fid,1,'real*8');
Re      = fread(fid,1,'real*8');
alp     = fread(fid,1,'real*8');
bet     = fread(fid,1,'real*8');
mpgx    = fread(fid,1,'real*8');
nband   = fread(fid,1,'int');
iter    = fread(fid,1,'int');
dummit  = fread(fid,90,'int');
N       = fread(fid,4*5,'int');
N       = reshape(N,4,5);
    fread(fid,1,'real*8'); %skip % fseek(fid,8,'cof');
y       = fread(fid,N(4,nband+1)-N(4,0+1)+2,'real*8');
dthetai = fread(fid,1,'real*8');
    fread(fid,1,'real*8'); %skip
dthdy   = fread(fid,N(4,nband+1)-N(4,0+1)+2,'real*8');


%     fread(fid,1,'int'); %skip
for index=2:N(4,1+1)-N(4,0+1)
    j1(index)        = fread(fid,1,'int');
    one1(index)      = fread(fid,1,'int');
    nx1(index)       = fread(fid,1,'int');
    nz1(index)       = fread(fid,1,'int');
    yj1(index)       = fread(fid,1,'real*8');
    p1(index,:)     = fread(fid,nx1(index)*nz1(index),'real*8');
                      fread(fid,1,'real*8');
end
for index=1:N(4,2+1)-N(4,1+1)+2
    j2(index)        = fread(fid,1,'int');
    one2(index)      = fread(fid,1,'int');
    nx2(index)       = fread(fid,1,'int');
    nz2(index)       = fread(fid,1,'int');
    yj2(index)       = fread(fid,1,'real*8');
    p2(index,:)     = fread(fid,nx2(index)*nz2(index),'real*8');
                      fread(fid,1,'real*8');
end
for index=1:N(4,3+1)-N(4,2+1)-1
    j3(index)        = fread(fid,1,'int');
    one3(index)      = fread(fid,1,'int');
    nx3(index)       = fread(fid,1,'int');
    nz3(index)       = fread(fid,1,'int');
    yj3(index)       = fread(fid,1,'real*8');
    p3(index,:)     = fread(fid,nx3(index)*nz3(index),'real*8');
                      fread(fid,1,'real*8');
end
fclose(fid);

p = [p1(:,1);p2(:,1);p3(:,1)];

%% Display
disp(['t=',num2str(t),' iter=',num2str(iter)]);
disp(['Re=',num2str(Re),' alp=',num2str(alp),' bet=',num2str(bet),...
    ' mpgx=',num2str(mpgx)]);
disp(['nband=',num2str(nband)]);
disp(N);
disp(['y =[',num2str(y(1)),',',num2str(y(end)),']']);
disp(['Umax=',num2str(max(p))]);

%% Plot
colors = [[0     ,0.4470,0.7410];...
          [0.8500,0.3250,0.0980];...
          [0.9290,0.6940,0.1250];...
          [0.4940,0.1840,0.5560];...
          [0.4660,0.6740,0.1880];...
          [0.3010,0.7450,0.9330];...
          [0.6350,0.0780,0.1840]];

figure(7);
%clf
%hold on
%axis equal;
plot(p(2:end),y(2:end-1),'linewidth',2,'color',colors(1,:));
title('p')
ylim([-1 1])

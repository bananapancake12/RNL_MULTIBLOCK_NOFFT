clear all; clc;

timematrix=[4886];
slice = 16;
camin=0;
camax=1.8;

alpha=1;
beta=2;

contourlvl=1

%slice = 14;
%slice = 9;

%camin=-0.0125;
%camax=0.0125;



r1=255;
g1=0;15;
b1=0;15;
r2=255;
g2=255;
b2=195;
r3=0;15;
g3=0;15;
b3=255;


r3=r2;
g3=g2;
b3=b2;

r2=(r3+r1)/2;
g2=(g3+g1)/2;
b2=(b3+b1)/2;


map=0;
map(1,1:3)=[r1,g1,b1];

for j = 1:127
    map(j,1:3)=[r1+(r2-r1)*j/128,g1+(g2-g1)*j/128,b1+(b2-b1)*j/128]/255;
end
map(128,1:3)=[r2,g2,b2]/255;
for j = 1:127
    map(128+j,1:3)=[r2+(r3-r2)*j/128,g2+(g3-g2)*j/128,b2+(b3-b2)*j/128]/255;
end

map=flipud(map);

for timeit=1:length(timematrix)

 
time=timematrix(timeit)  ;
    
%time = 00063;

%xx = 0128;
%xx = 0256;
%zz = 0256;
zz = 0768;
xx = 0128;
yy = 0153;
% xx = 0064;
% zz = 0064;
 %yy = 0033;

% yy = 0065;
%file = 'u2_0128x0256x0154_t00000.dat';

file = sprintf('output/u1_phys_%4.4dx%4.4dx%4.4d_t%5.5d.dat',xx,zz,yy,time);
%file

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
for index=1:N(4,2+1)-N(4,1+1)+1
    j2(index)        = fread(fid,1,'int');
    one2(index)      = fread(fid,1,'int');
    nx2(index)       = fread(fid,1,'int');
    nz2(index)       = fread(fid,1,'int');
    yj2(index)       = fread(fid,1,'real*8');
    u12(index,:)     = fread(fid,nx2(index)*nz2(index),'real*8');
                      fread(fid,1,'real*8');
end
for index=1:N(4,3+1)-N(4,2+1)+2
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



usurf=u11(contourlvl,:);
usurf=reshape(usurf,xx*3/2+2,zz*3/2);
usurf=usurf(1:xx*3/2,:);
xpl=[0:xx*3/2-1]/(xx*3/2)*2*3.14/alpha;
zpl=[0:zz*3/2-1]/(zz*3/2)*2*3.14/beta;
figure(1)
clf
pcolor(zpl,xpl,usurf)
colorbar
shading flat
axis equal
view(0,90)
%caxis([camin camax])
min(min(usurf))
max(max(usurf))

%colormap(map)

pause(0.1)



usurface=(u11(contourlvl+1,:)*(-1-yj1(1))+u11(contourlvl,:)*(yj1(2)+1))/(yj1(2)-yj1(1));

usurface=reshape(usurface,xx*3/2+2,zz*3/2);
usurface=usurface(1:xx*3/2,:);

ughost=u11(contourlvl,:);
ughost=reshape(ughost,xx*3/2+2,zz*3/2);
ughost=ughost(1:xx*3/2,:);

uplus1=u11(contourlvl+1,:);
uplus1=reshape(uplus1,xx*3/2+2,zz*3/2);
uplus1=uplus1(1:xx*3/2,:);

uplus2=u11(contourlvl+2,:);
uplus2=reshape(uplus2,xx*3/2+2,zz*3/2);
uplus2=uplus2(1:xx*3/2,:);

figure(2)
clf
plot(usurface(1:xx,1),'b-')
hold on
plot(ughost(1:xx,1),'r-')
plot(uplus1(1:xx,1),'c-')
plot(uplus2(1:xx,1),'m-')
figure(3)
clf
plot(usurface(1,:),'b-')
hold on
plot(ughost(1,:),'r-')
plot(uplus1(1,:),'c-')
plot(uplus2(1,:),'m-')


end
% for j=1:25
% usurf=u11(j,:);
% usurf=reshape(usurf,194,192);
% usurf=usurf(1:192,:);
% figure(1)
% clf
% surf(usurf)
% %caxis([0 0.15])
% colorbar 
% shading interp
% view(0,90)
% pause(0.25)
% end

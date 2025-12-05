clear all; clc;

%timematrix=[11157];
timematrix=[01857];

xx = 1152;%0128;
zz = 0576;
yy = 0153;

slice = 16;
camin=0;
camax=1.8;

alpha=1;
beta=2;

contourlvl=7

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

 
time=timematrix(1)  ;
    

file = sprintf('output/Qcrit_%4.4dx%4.4dx%4.4d_t%5.5d.dat',xx,zz,yy,time);
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
y       = fread(fid,N(3,nband+1)-N(3,0+1)+2,'real*8');
dthetai = fread(fid,1,'real*8');
    fread(fid,1,'real*8'); %skip
dthdy   = fread(fid,N(3,nband+1)-N(3,0+1)+2,'real*8');


     %fread(fid,1,'int'); %skip
for index=1:N(3,1+1)-N(3,0+1)+1
    j1(index)        = fread(fid,1,'int');
    one1(index)      = fread(fid,1,'int');
    nx1(index)       = fread(fid,1,'int');
    nz1(index)       = fread(fid,1,'int');
    yj1(index)       = fread(fid,1,'real*8');
    u11(index,:)     = fread(fid,nx1(index)*nz1(index),'real*8');
                      fread(fid,1,'real*8');
end
for index=1:N(3,2+1)-N(3,1+1)+1-1
    j2(index)        = fread(fid,1,'int');
    one2(index)      = fread(fid,1,'int');
    nx2(index)       = fread(fid,1,'int');
    nz2(index)       = fread(fid,1,'int');
    yj2(index)       = fread(fid,1,'real*8');
    u12(index,:)     = fread(fid,nx2(index)*nz2(index),'real*8');
                      fread(fid,1,'real*8');
end
for index=1:N(3,3+1)-N(3,2+1)+2-2+1
    j3(index)        = fread(fid,1,'int');
    one3(index)      = fread(fid,1,'int');
    nx3(index)       = fread(fid,1,'int');
    nz3(index)       = fread(fid,1,'int');
    yj3(index)       = fread(fid,1,'real*8');
    u13(index,:)     = fread(fid,nx3(index)*nz3(index),'real*8');
                      fread(fid,1,'real*8');
end
fclose(fid);

%u1 = [u11(:,:);u12(:,:);u13(:,:)];


Q11 = reshape(u11,length(j1),nx1(1),nz1(1));
Q12 = reshape(u12,length(j2),nx2(1),nz2(1));
Q13 = reshape(u13,length(j3),nx3(1),nz3(1));

for i = 1:nx1(1)-2
    for k = 1:nz1(1)
        for j = 1:length(j1)
            Q(i,k,j1(j)+1) = Q11(j,i,k);
        end
        for j = 1:length(j3)
            Q(i,k,j3(j)+1) = Q13(j,i,k);
        end
    end
end

for i = 1:nx2(1)-2
    for k = 1:nz2(1)
        for j = 1:length(j2)
            Q(i,k,j2(j)+1) = Q12(j,i,k);
        end
    end
end

for i = 1:nx2(1)-2
    for k = 1:nz2(1)
        for j = 1:yy/2+1
            Qplt(i,k,j) = Q(i,k,j);
        end
    end
end
return
[x,z,y] = meshgrid((0:nx2(1)-2-1)*2*pi/alp/(nx2(1)-2),(0:nz2(1)-1)*2*pi/bet/nz2(1),(0:yy/2)/(yy/2));

figure(2)
clf

hold on
%p = patch(isosurface(x,z,y,Q,5))
p = patch(isosurface(x,z,y,Qplt,5))
%p = patch(isosurface(x,z,y,Q,9))
isonormals(x,z,y,Qplt,p)
%isocolors(x,y,z,Q,p)
isocolors(x,z,y,y,p)
p.FaceColor = 'interp';
p.EdgeColor = 'none';
shading flat
view(60,20)
set(gca,'xcolor','k')
set(gca,'ycolor','k')
set(gca,'zcolor','k')
set(gcf,'Color','w')
set(gca,'Color','w')
daspect([1 1 1])
axis tight
camlight(0,0)
lighting GOURAUD
zlim([0 1])

xlabel('x')
ylabel('z')
zlabel('y')


nx = 1152;
nz = 576;


    ntilex = 96;
    ntilez = 48;
    Lfracx = 3;
    Lfracz = 3;

    planeBC = 2*ones(nx,nz); %2 - Free-shear

    texturepointsx = nx/ntilex;
    texturepointsz = nz/ntilez;
    
    postpointsx = (nx/ntilex/Lfracx)+1;
    postpointsz = (nz/ntilez/Lfracz)+1;
    
    for i=1:postpointsx
      for k=1:postpointsz
        planeBC(i,k) = 1; %1 - No-slip
      end
    end 
  
    
    for ip=1:ntilex
      for kp=1:ntilez
	for i=1:postpointsx
	  for k=1:postpointsz
	    planeBC(i+(ip-1)*texturepointsx,k+(kp-1)*texturepointsz) = planeBC(i,k);
      end        
    end
      end
    end
    

  
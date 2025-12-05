clear

time = 07417;
xx = 0128;
zz = 0128;
yy = 0153;


file = sprintf('output/u3_%4.4dx%4.4dx%4.4d_t%5.5d.dat',xx,zz,yy,time);

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

u31 = reshape(u31,nx1(1),nz1(1),N(4,1+1)-N(4,0+1)+1);
u32 = reshape(u32,nx2(1),nz2(1),N(4,2+1)-N(4,1+1));
u33 = reshape(u33,nx3(1),nz3(1),N(4,3+1)-N(4,2+1)+1);

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



u21 = reshape(u21,nx1(1),nz1(1),N(3,1+1)-N(3,0+1)+1);
u22 = reshape(u22,nx2(1),nz2(1),N(3,2+1)-N(3,1+1));
u23 = reshape(u23,nx3(1),nz3(1),N(3,3+1)-N(3,2+1)+1);



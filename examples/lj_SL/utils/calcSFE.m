clear all;

s1="step1_work.dat";
s2="step2_work.dat";
s3="step3_work.dat";
s4="step4_work.dat";

fileout=s1;

formatSpec = '%f %f %f';
sizeA = [3 inf];

fileID1 = fopen(fileout,'r');
A1=fscanf(fileID1,formatSpec,sizeA);
x=A1(1,:);
y=A1(2,:);
z=A1(3,:);



Q1 = trapz(x,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step ='step2';

fileout=s2;

fileID2 = fopen(fileout,'r');
A2=fscanf(fileID3,formatSpec,sizeA);

x=A2(1,:);
y=A2(2,:);
z=A2(3,:);

Q2 = trapz(x,y);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step ='step3';

fileout=s3;

fileID3 = fopen(fileout,'r');
A3=fscanf(fileID3,formatSpec,sizeA);

x=A3(1,:);
y=A3(2,:);
z=A3(3,:);

Q3 = trapz(x,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

step ='step4';

fileout=s4;

fileID4 = fopen(fileout,'r');
A4=fscanf(fileID4,formatSpec,sizeA);
x=A4(1,:);
y=A4(2,:);
z=A4(3,:);


Q4 = trapz(x,y);

workMat = [Q1 Q2 Q3 Q4]
work = Q1 + Q2 + Q3 + Q4




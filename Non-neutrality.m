grid on;

% area under the base station
r = 250;
t = 0 : pi/1000 : pi*2;
x = r*cos(t);
y = r*sin(t);
figure(1);
plot (x,y);
hold on;
grid on;
axis on;
axis ([-300 300 -300 300])
xlabel('users under BS, red denoting Normal users and yellow for Priority users' );


% layout of the base station
Bs_X = 0;
Bs_Y = 0;
plot(Bs_X, Bs_Y, '^g');


% genration of users in base station
a=30;
b=240;
n=1000; 
OpUserR = a + (b-a).*rand(n,1);
c=0;
d=2*pi;
theta = c + (d-c).*rand(n,1);
userX= (OpUserR.*cos(theta));
userY= (OpUserR.*sin(theta));
normal = 0.5 *n;
priority = n - normal;
userX1= zeros(1,(normal));
userY1= zeros(1,(normal));
userX2= zeros(1,(priority));
userY2= zeros(1,(priority));
for i = 1:(normal)
    userX1(i) = userX(i);
    userY1(i) = userY(i);
end
k = 1;
for i = (normal+1):n
    userX2(k) = userX(i);
    userY2(k) = userY(i);
    k = k+1;
end

figure(1);
hold on;
plot(userX1, userY1,'+r')
plot(userX2, userY2,'+y')


%parameter calculation for users in the network
d=zeros(1,n);
pathLoss_M=zeros(1,n);
pathLoss2_M=zeros(1,n);
shadowing_M=zeros(1,n);
hbs=25;
hms=1.5;
fc=2.6;
std_dev_M= 8;
for i= 1:n
    d(i)= sqrt((userX(i)^2) +(userY(i)^2));
    pathLoss_M(i)= (44.9-6.55*log10(hbs))*log10(d(i))+34.46+5.83*log10(hbs)+23*log10(fc/5);
    pathLoss2_M(i)= 40*log10(d(i)) + 0.65 - 16.2*log10(hbs) - 16.2*log10(hms) + 3.8*log10(fc/5); 
    shadowing_M= std_dev_M.* randn(n, 1);
end;




% parameter calculation for the first sector
i=1;
Aa1 = zeros(1,n) ;
Mg = 17; % max gain
Am = 20; %
W = 70;
gain = zeros(3, n);

while(i<=n)
if (0<=userX(i) && 0<userY(i))
  Aa1(i) = abs(atand(userY(i)/userX(i)))  ;

                 
elseif (0>userX(i) && 0<userY(i))
  Aa1(i) = 180 - abs(atand(userY(i)/userX(i)))  ;

elseif (0>userX(i) && 0>userY(i))
  Aa1(i) = 180 + abs(atand(userY(i)/userX(i)))  ;

elseif (0<userX(i) && 0>userY(i))
  Aa1(i) = 360 - abs(atand(userY(i)/userX(i)))  ;

end;
  if Aa1(i)>180
      Aa1(i)= 360- Aa1(i);
  end;
  if Aa1(i)< -180
      Aa1(i)= 360+Aa1(i);
  end;
      
   

 x=12*((Aa1(i)/W)^2);
gain(1,i)= Mg - min(x, Am);

i =i + 1 ;
end;



% parameter calculation for the second sector
i=1;
Aa2 = zeros(1,n) ;

while(i<=n)
if (0<userX(i) && 0<userY(i)) 
  Aa2(i) = 120 - abs(atand(userY(i)/userX(i)))  ;

elseif (0>userX(i) && 0<userY(i) && theta(i)>2*pi/3)
  Aa2(i) = 60 - abs(atand(userY(i)/userX(i)));
  
elseif (0>userX(i) && 0<userY(i) && theta(i)<2*pi/3)
  Aa2(i) = abs(atand(userY(i)/userX(i)))-60;

elseif (0>userX(i) && 0>userY(i))
  Aa2(i) = 60 + abs(atand(userY(i)/userX(i)))  ;

else 
  Aa2(i) = 240 - abs(atand(userY(i)/userX(i)))  ;

end;
 if Aa2(i)>180
      Aa2(i)= 360- Aa2(i);
  end;
  if Aa2(i)< -180
      Aa2(i)= 360+Aa2(i);
  end;
  
 x2=12*((Aa2(i)/W)^2);
gain(2,i)= Mg - min(x2, Am);

i =i + 1 ;
end;

% parameter calculation for the third sector
i=1;
Aa3 = zeros(1,n) ;

while(i<=n)
if (0<userX(i) && 0<userY(i))
  Aa3(i) = 240 - abs(atand(userY(i)/userX(i)))  ;


elseif (0>userX(i) && 0<userY(i))
  Aa3(i) = 60 + abs(atand(userY(i)/userX(i))) ;

elseif (0>userX(i) && 0>userY(i))
  theta1 =abs(atand(userY(i)/userX(i)) ) ;
  if theta1<60
      Aa3(i)= 60- theta1;
  elseif theta1>=60
      Aa3(i)=theta1-60;
  end;
  
 elseif (0<userX(i) && 0>userY(i))
  Aa3(i) = 120 - abs(atand(userY(i)/userX(i)))  ;

 end;
 
  if Aa3(i)>180
      Aa3(i)= 360- Aa3(i);
  end;
  if Aa3(i)< -180
      Aa3(i)= 360+Aa3(i);
  end;
      
 x3=12*((Aa3(i)/W)^2);
gain(3,i)= Mg - min(x3, Am);

i =i + 1 ;
end;


RS= zeros(3,n);
RSw= zeros(3,n);% Received signal in Watt
BS_Tx= 43;
gainMs= 0;
SNRw= zeros(3,n);
SNRdB= zeros(3,n);
bandwidth= 5000000;
noise_dB= -174 + 10*log10(bandwidth);
noise_w= 10^(noise_dB/10); % Noise power in Watt

alpha = 0.65;
dataRatePerH = zeros(1,n);
dataRate = zeros(1,n);
uSNRw= zeros(3,n);
uSNRdB = zeros(1,n);
RBbandwidth = 180000;
Unoise_dB= -174 + 10*log10(RBbandwidth);
Unoise_w= 10^(Unoise_dB/10);
fileSize = 16777216; % 2MB in bits

for i= 1:3
    for j= 1:n
        RS(i,j)= BS_Tx + gain(i,j) + gainMs - pathLoss_M(j) - shadowing_M(j);
        RSw(i,j)= 10^(RS(i,j)/10);
        SNRw(i,j)= RSw(i,j)/noise_w;
        SNRdB(i,j)= 10 *log10(SNRw(i,j));
        
        
        
        uSNRw(i,j)= RSw(i,j)/Unoise_w;
        uSNRdB(i,j)= 10 *log10(uSNRw(i,j));
        if uSNRdB(i,j)>21
            uSNRdB(i,j) =21;
        end
        dataRatePerH = alpha*log2(1+10.^(uSNRdB/10));
        dataRate = dataRatePerH*RBbandwidth; 
        
       
    end;
end;


RdataRate = max(dataRate,[],1); % data rate(tx rate ) of each user.
serviceTime = zeros(1,n); % Service time for each user.

i=1;
for i=i:n
serviceTime(i) = fileSize/(RdataRate(i));
end;

SNRdB1= zeros(3,normal);
SNRdB2= zeros(3,priority);

for i= 1:3
    for j= 1:normal
        SNRdB1(i,j)= SNRdB(i,j);
    end;
end;

k=1;
for i= 1:3
    for j= (normal+1):n
        SNRdB2(i,k)= SNRdB(i,j);
        k=k+1;
    end;
    k=1;
end;


[maxSNRdB, indexSNRdB] =  max(SNRdB,[],1);
cell1 = find(indexSNRdB==1); %index of All users connected to cell1
cell2 = find(indexSNRdB==2); %index of All users connected to cell2
cell3 = find(indexSNRdB==3); %index of All users connected to cell3

[maxSNRdB1, indexSNRdB1] =  max(SNRdB1,[],1);
[maxSNRdB2, indexSNRdB2] =  max(SNRdB2,[],1);
cell11 = find(indexSNRdB1==1); %index of Normal users connected to cell1
cell12 = find(indexSNRdB1==2); %index of Normal users connected to cell2
cell13 = find(indexSNRdB1==3); %index of Normal users connected to cell3

cell21 = find(indexSNRdB2==1); %index of Priority users connected to cell1
cell22 = find(indexSNRdB2==2); %index of Priority users connected to cell1
cell23 = find(indexSNRdB2==3); %index of Priority users connected to cell1

figure(2);
r = 250;
t = 0 : pi/1000 : pi*2;
x = r*cos(t);
y = r*sin(t);
plot (x,y);
hold on;

figure(2);
plot(userX(cell1), userY(cell1),'+r', userX(cell2), userY(cell2),'+g',userX(cell3), userY(cell3),'+b');
xlabel('connetion to anntenna based on SNRdB for three Sectors' );


% Radio Resource Allocation Segment

Lamda = (0.5:1:12.5);
z = length(Lamda);
blockingProb = zeros(1,z);
throughput = zeros(1,z);
bNum = zeros(1,z); %number of blocked users
sCall = zeros(1,z);%number of served users
timeNow = -1*ones(1,z);
averageDelay = zeros(1,z);
delay = zeros(1,z);
numChannelsPerSector = 25;
perNormal = 0.5; % Percentage of total bandwidth allocated to Normal users
perPriority = 1- perNormal; % Percentage of total bandwidth allocated to Priority users
numChannelsNormal =ceil( numChannelsPerSector * perNormal);
numChannelsPriority = numChannelsPerSector - numChannelsNormal;
l=1; %  Iterates up to the length of lambda (z)


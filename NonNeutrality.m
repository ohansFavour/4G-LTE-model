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
fileSize = 16000000; % 2MB in bits

for i= 1:3
    for j= 1:normal
        RS(i,j)= BS_Tx + gain(i,j) + gainMs - pathLoss_M(j) - shadowing_M(j);
        RSw(i,j)= 10^(RS(i,j)/10);
        SNRw(i,j)= RSw(i,j)/noise_w;
        SNRdB(i,j)= 10 *log10(SNRw(i,j));
        
        uSNRw(i,j)= RSw(i,j)/Unoise_w;
        uSNRdB(i,j)= 10 *log10(uSNRw(i,j));
        
        if uSNRdB(i,j)>21
            uSNRdB(i,j) =21;
        end
        
        dataRatePerH(i,j) = alpha*log2(1+10.^(uSNRdB(i,j)/10));
        dataRate(i,j) = 900000; % 900Kbps
        
    end;
    
    
       for j= normal+1 : n
        RS(i,j)= BS_Tx + gain(i,j) + gainMs - pathLoss_M(j) - shadowing_M(j);
        RSw(i,j)= 10^(RS(i,j)/10);
        SNRw(i,j)= RSw(i,j)/noise_w;
        SNRdB(i,j)= 10 *log10(SNRw(i,j));
        
        uSNRw(i,j)= RSw(i,j)/Unoise_w;
        uSNRdB(i,j)= 10 *log10(uSNRw(i,j));
        
        if uSNRdB(i,j)>21
            uSNRdB(i,j) =21;
        end
        
        dataRatePerH(i,j) = alpha*log2(1+10.^(uSNRdB(i,j)/10));
        dataRate(i,j) = 2500000; % 2.5Mbps
        
       end;
    
    
end;



RdataRate = max(dataRate,[],1); % data rate(tx rate ) of each user.
RBbandwidthNormal = RdataRate(1,1)/dataRatePerH(1,1); % Adjusted to fit 900kb/s data rate
RBbandwidthPriority = RdataRate(1,normal+1)/dataRatePerH(1,normal+1); % Adjusted to fit 2.5MB/s data rate
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

Lamda = 0.5:1.5:12.5;
z = length(Lamda);
blockingProb = zeros(1,z);
throughput = zeros(1,z);
Energy = zeros(1,z);
ECR = zeros(1,z);
bNum = zeros(1,z); %number of blocked users
sCall = zeros(1,z);%number of served users
timeNow = -1*ones(1,z);
averageDelay = zeros(1,z);
delay = zeros(1,z);
numChannelsPerSector = 32;
perNormal = 0.5; % Percentage of total bandwidth allocated to Normal users
perPriority = 1- perNormal; % Percentage of total bandwidth allocated to Priority users
numChannelsNormal = ceil( numChannelsPerSector * perNormal); % Number of channels available for normal user
numChannelsPriority = numChannelsPerSector - numChannelsNormal; % Number of channels available for priority user
l=1; %  Iterates up to the length of lamda (z)
p_max = 20; %in watt
p_not = 130; %in watt
deltaP = 4.7;
Ntrx = 1;


for lamda = 0.5:1.5:12.5;  % Adjust one on top too
    
 tNumCalls = 10000;  
 lamdaNormal = 0.5 * lamda;
 lamdaPriority = lamda - lamdaNormal;

 XNormal =  -1/lamdaNormal*log(1-rand(1,tNumCalls/2)); % Interarrival for Normal
 arrivalTimeNormal = cumsum(XNormal);  % vector of arrival times for Normal 

 XPriority =  -1/lamdaPriority*log(1-rand(1,tNumCalls/2)); % Interarrival for Priority
 arrivalTimePriority = cumsum(XPriority);  % vector of arrival times for Priority

 combineArrivalTime = sort(union(arrivalTimeNormal,arrivalTimePriority)); % in Ascending form
 
 departureTime = inf*ones(1,tNumCalls);
 depTime = inf*ones(1,n); % for each user in the system
 
 tNum = 1;  % Current number of arrived calls
 channelStatus1 = zeros(1,numChannelsPerSector);
 channelStatus2 = zeros(1,numChannelsPerSector);
 channelStatus3 = zeros(1,numChannelsPerSector);
 pp=0;
 
 cellSel = zeros(1,n);  % Cell selected by user
 uState = zeros(1,n);   % Current state of each user
 uStateNormal = zeros(1,normal); % state of Normal Users
 uStatePriority = zeros(1,priority); % state of Priority Users
 uChannel = zeros(1,n); % Current channel connected to
 userAtime = zeros(1,n);% Current user arrival time
 userDtime = zeros(1,n);% current user departure time
 minaggTime =0;
 
 
  while(tNum<=tNumCalls-1)
     
    idle = 0;
    busy = 1;
    
    % Check for the class corresponding to the current Arrival time
    normalArrival = any(arrivalTimeNormal(:) == combineArrivalTime(tNum) );
    priorityArrival = any(arrivalTimePriority(:) == combineArrivalTime(tNum) );
    
    
    if normalArrival == 1    % if user is from the Normal class
        allFree = find(uState == 0); % vector of indexes of all free users
        freeUsers = find(allFree <= normal);
        uCall = allFree(freeUsers(randi(numel(freeUsers))));   % index of selected user to make call  
    end
    
    if priorityArrival == 1    % if user is from the Priority class
        allFree = find(uState == 0); % vector of indexes of all free users
        freeUsers = find(allFree > normal);
        uCall = allFree(freeUsers(randi(numel(freeUsers))));  % index of selected user to make call  
    end
     
    userAtime(uCall) = combineArrivalTime(tNum); % Arrival Time of selected user
    
    id1 = find(cell1==uCall);
    id2 = find(cell2==uCall);
    id3 = find(cell3==uCall); 
    
    TF1 = isempty(id1);
    TF2 = isempty(id2);
    TF3 = isempty(id3); 
    
   
     if TF1==0  % if user to call is connected to cell1
        if (1<=uCall && uCall<= normal) % if user is a Normal user
            selCh = find(channelStatus1==0); 
            freeCh = find(selCh <= numChannelsNormal);
            cellSel(uCall)=11; % normal cell1
            
        else  % user is a priority user
            selCh = find(channelStatus1==0); 
            freeCh = find(selCh > numChannelsNormal);
            cellSel(uCall)=21; % priority cell1 
        end
        
    end;
    
    if TF2==0  % if user to call is connected to cell2
        if (1<=uCall && uCall<= normal) % if user is a Normal user
            selCh = find(channelStatus2==0); 
            freeCh = find(selCh <= numChannelsNormal);
            cellSel(uCall)=12; % normal cell2
            
        else  % user is a priority user
            selCh = find(channelStatus2==0); 
            freeCh = find(selCh >numChannelsNormal);
            cellSel(uCall)=22; % priority cell2 
        end
    end;
    
    if TF3==0  % if user to call is connected to cell3
         if (1<=uCall && uCall<= normal) % if user is a Normal user
            selCh = find(channelStatus3==0); 
            freeCh = find(selCh <= numChannelsNormal);
            cellSel(uCall)=13; % normal cell3
            
        else  % user is a priority user
            selCh = find(channelStatus3==0); 
            freeCh = find(selCh > numChannelsNormal);
            cellSel(uCall)=23; % priority cell3 
        end
    end;
     
    h = isempty(freeCh);
    if RdataRate(uCall) >= 155.5740 %corresponding datarate for 1.8dB user snrdB
        m = 0;
    end;
    
   
    if h==0 && m==0   % free channel and high enough snrDB
        
            uState(uCall) = busy;
            uChannel(uCall)= selCh(freeCh(1));
        
        
            if (cellSel(uCall) == 11 || cellSel(uCall) == 21)
                
                 %Energy Modelling
              busyRatio1 = length(find(channelStatus1 == 1))/length(channelStatus1);
              busyRatio2 = length(find(channelStatus2 == 1))/length(channelStatus2);
              busyRatio3 = length(find(channelStatus3 == 1))/length(channelStatus3);
           
              p_out1 = busyRatio1 * p_max;
              p_out2 = busyRatio2 * p_max;
              p_out3 = busyRatio3 * p_max;
              power = Ntrx*((p_not + deltaP*p_out1)+ (p_not + deltaP*p_out2)+ (p_not + deltaP*p_out3)) ; 
              
              if(indexaggTime == 1)
              Energy(l) = Energy(l) + power * (combineArrivalTime(tNum) - combineArrivalTime(tNum-1) );
              else
                  Energy(l) = Energy(l) + power * (combineArrivalTime(tNum) - minaggTime );
              end  
              
              
              channelStatus1(selCh(freeCh(1))) = busy;
              
            elseif (cellSel(uCall) == 12 || cellSel(uCall) == 22)
                
                  %Energy Modelling
              busyRatio1 = length(find(channelStatus1 == 1))/length(channelStatus1);
              busyRatio2 = length(find(channelStatus2 == 1))/length(channelStatus2);
              busyRatio3 = length(find(channelStatus3 == 1))/length(channelStatus3);
         
              p_out1 = busyRatio1 * p_max;
              p_out2 = busyRatio2 * p_max;
              p_out3 = busyRatio3 * p_max;
              power = Ntrx*((p_not + deltaP*p_out1)+ (p_not + deltaP*p_out2)+ (p_not + deltaP*p_out3)) ;  
             
              if(indexaggTime == 1)
              Energy(l) = Energy(l) + power * (combineArrivalTime(tNum) - combineArrivalTime(tNum-1) );
              else
                  Energy(l) = Energy(l) + power * (combineArrivalTime(tNum) - minaggTime );
              end  
              
                
              channelStatus2(selCh(freeCh(1))) = busy;
              
            elseif (cellSel(uCall) == 13 || cellSel(uCall) == 23) 
                
                  %Energy Modelling
              busyRatio1 = length(find(channelStatus1 == 1))/length(channelStatus1);
              busyRatio2 = length(find(channelStatus2 == 1))/length(channelStatus2);
              busyRatio3 = length(find(channelStatus3 == 1))/length(channelStatus3);
           
              p_out1 = busyRatio1 * p_max;
              p_out2 = busyRatio2 * p_max;
              p_out3 = busyRatio3 * p_max;
              power = Ntrx*((p_not + deltaP*p_out1)+ (p_not + deltaP*p_out2)+ (p_not + deltaP*p_out3)) ;  
              
              if(indexaggTime == 1)
              Energy(l) = Energy(l) + power * (combineArrivalTime(tNum) - combineArrivalTime(tNum-1) );
              else
                  Energy(l) = Energy(l) + power * (combineArrivalTime(tNum) - minaggTime );
              end  
              
                
            channelStatus3(selCh(freeCh(1))) = busy;
            end;
    
            departureTime(tNum) = combineArrivalTime(tNum) + serviceTime(uCall);
            depTime(uCall) = departureTime(tNum);uState(uCall) = busy;
            uChannel(uCall)= selCh(freeCh(1));
        
        
            if (cellSel(uCall) == 11 || cellSel(uCall) == 21)
            channelStatus1(selCh(freeCh(1))) = busy;
            elseif (cellSel(uCall) == 12 || cellSel(uCall) == 22)
            channelStatus2(selCh(freeCh(1))) = busy;
            elseif (cellSel(uCall) == 13 || cellSel(uCall) == 23) 
            channelStatus3(selCh(freeCh(1))) = busy;
            end;
    
            departureTime(tNum) = combineArrivalTime(tNum) + serviceTime(uCall);
            depTime(uCall) = departureTime(tNum);
        
    
    else bNum(l) = bNum(l)+1; 
    end;
    
 
    
    
     SdepartureTime = min(departureTime(departureTime>timeNow(l)));
     SarrivalTime = combineArrivalTime(tNum+1);
     aggTime = [SarrivalTime;SdepartureTime];
     [minaggTime, indexaggTime] =  min(aggTime,[],1);
     
     
      if indexaggTime == 1 % Arrival
          
          %Energy Modelling
          
              busyRatio1 = length(find(channelStatus1 == 1))/length(channelStatus1);
              busyRatio2 = length(find(channelStatus2 == 1))/length(channelStatus2);
              busyRatio3 = length(find(channelStatus3 == 1))/length(channelStatus3);
           
              p_out1 = busyRatio1 * p_max;
              p_out2 = busyRatio2 * p_max;
              p_out3 = busyRatio3 * p_max;
              power = Ntrx*((p_not + deltaP*p_out1)+ (p_not + deltaP*p_out2)+ (p_not + deltaP*p_out3)) ; 
              Energy(l) = Energy(l) + power * (SarrivalTime - timeNow(l) );
              
      end
  
     
     

     
      if indexaggTime == 2 % Departure Time
         TL1 = find(depTime<=SarrivalTime);
         TL2 = find(uState==busy);
         idUser = intersect(TL1,TL2);
         q = length(idUser); 
         kk = zeros(1,q);
         i=1;
         ii=1;
       
         for ii = ii:q
             kk(ii) = cellSel(idUser(ii));
         end; 
         
         initTime = timeNow(l); 
       
         while i<=q
             uState(idUser(i))=idle;% reseting a user's status after a successful call
             userDtime(idUser(i)) = depTime(idUser(i)); 
             delay(l) = delay(l) + userDtime(idUser(i)) - userAtime(idUser(i));
             sCall(l) = sCall(l)+1;
               
            if kk(i) ==11 || kk(i) == 21
              
                %Energy Modelling
              busyRatio1 = length(find(channelStatus1 == 1))/length(channelStatus1);
              busyRatio2 = length(find(channelStatus2 == 1))/length(channelStatus2);
              busyRatio3 = length(find(channelStatus3 == 1))/length(channelStatus3);
            
              p_out1 = busyRatio1 * p_max;
              p_out2 = busyRatio2 * p_max;
              p_out3 = busyRatio3 * p_max;
              power = Ntrx*((p_not + deltaP*p_out1)+ (p_not + deltaP*p_out2)+ (p_not + deltaP*p_out3)) ; 
              Energy(l) = Energy(l) + power * (depTime(idUser(i)) - initTime );
              initTime = depTime(idUser(i));
              
              channelStatus1(uChannel(idUser(i)))=idle;
         
            elseif kk(i)==12 || kk(i) == 22
                
               %Energy Modelling
              busyRatio1 = length(find(channelStatus1 == 1))/length(channelStatus1);
              busyRatio2 = length(find(channelStatus2 == 1))/length(channelStatus2);
              busyRatio3 = length(find(channelStatus3 == 1))/length(channelStatus3);
           
              p_out1 = busyRatio1 * p_max;
              p_out2 = busyRatio2 * p_max;
              p_out3 = busyRatio3 * p_max;
              power = Ntrx*((p_not + deltaP*p_out1)+ (p_not + deltaP*p_out2)+ (p_not + deltaP*p_out3)) ;  
              Energy(l) = Energy(l) + power * (depTime(idUser(i)) - initTime );
              initTime = depTime(idUser(i));
              
              channelStatus2(uChannel(idUser(i)))=idle; 
         
            elseif kk(i)==13 || kk(i) == 23
                
                   %Energy Modelling
              busyRatio1 = length(find(channelStatus1 == 1))/length(channelStatus1);
              busyRatio2 = length(find(channelStatus2 == 1))/length(channelStatus2);
              busyRatio3 = length(find(channelStatus3 == 1))/length(channelStatus3);
             
              p_out1 = busyRatio1 * p_max;
              p_out2 = busyRatio2 * p_max;
              p_out3 = busyRatio3 * p_max;
              power = Ntrx*((p_not + deltaP*p_out1)+ (p_not + deltaP*p_out2)+ (p_not + deltaP*p_out3)) ;  
              Energy(l) = Energy(l) + power * (depTime(idUser(i)) - initTime );
              initTime = depTime(idUser(i));
              
              channelStatus3(uChannel(idUser(i)))=idle;
            end;
            depTime(idUser(i))= inf;
            pp=pp+1;
            i=i+1;
         end;
      
     
     end;
    
   timeNow(l) = minaggTime;
   tNum = tNum + 1;
    
    
  end;
 
 
 blockingProb(l) = (bNum(l)/(tNum-1))*100;

 throughput(l) = (sCall(l)*fileSize)/timeNow(l);

 averageDelay(l) = delay(l)/sCall(l) ;
 
 ECR(l) = Energy(l)/(sCall(l)*fileSize);
 
 l = l+1;
 
end;

figure(3);
hold on;
grid on;
plot(Lamda, blockingProb,'^-b');
xlabel('arrival rate' );
ylabel('blocking pobability (%)' );


figure(4);
hold on;
grid on;
plot(Lamda, throughput, '^-b');
xlabel('arrival rate' );
ylabel('throughput (bit/s)' );


figure(5);
hold on;
grid on;
plot(Lamda, averageDelay, '^-b');
xlabel('arrival rate' );
ylabel('average delay (s)' );

figure(6);
hold on;
grid on;
plot(Lamda, ECR, '^-b');
xlabel('arrival rate' );
ylabel('Energy Consumption Ratio (Joules/bit)' );





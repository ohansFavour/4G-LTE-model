numCalls = 5000; % total number of calls
erlang= 0.5:1.5:12.5;  % Arrival Rate
serviceTimeNormal = 17.7778;
serviceTimePriority = 17.7778;


NNormal = 45; % total number of channels available to the normal users
normalErlang = 0.5; % Ratio of the arrival rate
erlangNormal = erlang*normalErlang; % Normal subscriber's Arrival rate
PbNormal = zeros(1,length(erlangNormal)); % Blocking Probability for the normal users
ErlgNormal = zeros(1,length(erlangNormal));
throughputNormal = zeros(1,length(erlangNormal));% Throughput for the normal users
yNormal = 1;


for Erlang = (erlang)*normalErlang
num = ((Erlang*serviceTimeNormal)^NNormal)/ factorial(NNormal); % Numerator of the formula
syms k ;
dem = symsum((((Erlang*serviceTimeNormal)^k)/factorial(k)), k, 0, NNormal);

PbNormal(yNormal) = num/dem;
ErlgNormal(yNormal)= Erlang;
throughputNormal(yNormal) = Erlang*8*2*(1000^2)* (1-PbNormal(yNormal));
yNormal = yNormal+1;
end

NbNormal= (PbNormal *(numCalls/2)); % Total number of blocked calls for Normal








NPriority = 45;   % total number of channels available to the Priority users
erlangPriority = (Erlang)*0.5;   % Ratio of the arrival rate
PbPriority = zeros(1,length(erlangPriority));
ErlgPriority = zeros(1,length(erlangPriority));
throughputPriority = zeros(1,length(erlangPriority));
yPriority = 1;
priorityErlang = 1- normalErlang;

for Erlang = (erlang)*priorityErlang
num = ((Erlang*serviceTimePriority)^NPriority)/ factorial(NPriority);
syms k ;
dem = symsum((((Erlang*serviceTimePriority)^k)/factorial(k)), k, 0, NPriority);

PbPriority(yPriority) = num/dem;
ErlgPriority(yPriority)= Erlang;
throughputPriority(yPriority) = Erlang*8*2*(1000^2)* (1-PbPriority(yPriority));
yPriority = yPriority+1;
end

NbPriority = PbPriority*(numCalls/2);  % blocked priority calls

PbTotal = ((NbNormal + NbPriority)/numCalls)*100;
ErlgTotal =  ErlgNormal + ErlgPriority ;
throughputTotal = throughputNormal + throughputPriority;



figure(3)
hold on;
grid on;
plot(erlang, PbTotal,'-or');

figure(4)
hold on;
grid on;
plot(erlang, throughputTotal, '-or');





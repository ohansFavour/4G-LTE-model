numCalls = 1000000;
erlang= 0.5:1.5:12.5;
serviceTimeNormal = 17.7778;
serviceTimePriority = 6.4;


NNormal = 45;
normalErlang = 0.5;
erlangNormal = Erlang*normalErlang;
PbNormal = zeros(1,length(erlangNormal));
ErlgNormal = zeros(1,length(erlangNormal));
throughputNormal = zeros(1,length(erlangNormal));
yNormal = 1;


for Erlang = (erlang)*normalErlang
num = ((Erlang*serviceTimeNormal)^NNormal)/ factorial(NNormal);
syms k ;
dem = symsum((((Erlang*serviceTimeNormal)^k)/factorial(k)), k, 0, NNormal);

PbNormal(yNormal) = num/dem;
ErlgNormal(yNormal)= Erlang;
throughputNormal(yNormal) = Erlang*8*2*(1000^2)* (1-PbNormal(yNormal));
yNormal = yNormal+1;
end

NbNormal= PbNormal *(numCalls/2);

NPriority = 45;
erlangPriority = (Erlang)*0.5;
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

NbPriority = PbPriority*(numCalls/2);

PbTotal = (NbNormal + NbPriority)/numCalls;
ErlgTotal =  ErlgNormal + ErlgPriority ;
throughputTotal = throughputNormal + throughputPriority;



figure(3)
hold on;
grid on;
plot(ErlgTotal, PbTotal*100,'-*g');

figure(4)
hold on;
grid on;
plot(ErlgTotal, throughputTotal, '-*g');


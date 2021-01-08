clear;
%% input and target
input = linspace(-2,2,10);
target =1+sin(input*(pi/8));

% mapping
[targetMinMax,mapping] = mapminmax(target,0,1);

%% create network (one hidden layer with 3 or 15 nodes)
%net.layers{1}.transferFcn = 'hardlim';
net = newfit(input, targetMinMax, [3], {exp(-(input).^2)});
net.trainParam.epochs = 3;
net.trainParam.lr = 0.15;

view(net)

%% training
net = init(net);                            % init
[net,tr] = train(net, input, targetMinMax ); % train
output = sim(net, input);                   % predict

%% view prediction
plot(input, mapminmax('reverse', output, mapping), 'r', 'linewidth',2), hold on
plot(input, target, 'o')
%plot(input, sin(input), 'g')
hold off
legend({'predicted' 'target' 'sin()'})
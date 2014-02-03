function [ fmin, ...        % function value of xmin
    FES, ...                % number of function evaluations done
    xmin, ...               % minimum search point of last iteration
    bestever ...            % overall best search point
    ] = SamACO( ...
    dimension, ...          % scalar shows the dimension count of arguments for fitfunc, F
    F, ...                  % handle of objective/fitness function
    Lower, ...              % scalar shows lower bound of the arguments
    Upper, ...              % scalar shows upper bound of the arguments
    inopts )                % struct with options, see defopts below
% SamACO.m, version 1.00, last change: 2014/2/3
% SamACO implements a Sampling Ant Colony Optimization for nonlinear
% function minimization.
%
% Input arguments: 
%   dimension is a scalar, indicates the dimension count of the arguments
%   for objective/fitness function.
%
%   F can be a function handle like @myfunc. F takes X as argument. X is a
%   n*m arguments matrix. n is the number of dimension of the
%   arguments. the F function can compute different m argument column
%   vectors at the same time.
%
%   Lower is a scalar, indicates the lower bound of the arguments
%
%   Upper is a scalar, indicates the upper bound of the arguments
%
%   inopts (an optional argument) is a struct holding additional input
%   options for SamACO. Valid field names and a short documentation can be
%   discovered by looking at the default options. Empty or missing fields
%   in inopts invoke the default value, i.e. inopts needs not to have all
%   valid field names.
% Output arguments: 
%   fmin, the minimum function value.
%
%   FES, number of function evaluations done.
%   
%   xmin, minimum search point of last iteration.
%
%   bestever, overall best search point.
%
% 输入：
%   dimension：作为标量，表明目标函数的维度数
%
%   F：目标函数的句柄，比如@myfun. F函数以X作为参数，X是n*m的参数矩阵。n是参数
%   的维度数，函数F可以同时计算m组参数列向量。
%
%   Lower：作为标量，是参数的下界
%
%   Upper：作为标量，是参数的上界。对于不同维度的参数的上下界是一样的
%
%   inopts：作为结构体，含有SamACO的各项参数，具体元素详见defopts。inopts中没有
%   的参数将使用默认的数值代替


%% Options
defopts.n           = dimension;        %维度
defopts.m           = 40;               %蚂蚁的数量20, 30, 40
defopts.v           = 40;               %局部搜索的次数
defopts.theta       = 1;                %随机搜索的次数，替换最差的theta个解
defopts.T0          = 0.1;              %信息素初始化的量
defopts.rou         = 0.5;              %信息素挥发常数
defopts.T_min       = defopts.T0;       %信息素最小值
defopts.T_max       = 1.0;              %信息素最大值
defopts.alpha       = 0.3;              %信息素加强系数
defopts.phi         = 1;                %加强phi个最优解的信息素
defopts.Vr          = 0.7;              %步长缩小常量
defopts.Ve          = 1 / defopts.Vr;   %步长增大常量
defopts.q0          = 0.1;              %蚂蚁选择某条路径的概率
defopts.times_max   = 600;              %最大迭代次数600
defopts.eps         = 1e-8;

if nargin < 5 || isempty(inopts)
    inopts = []; 
    opts = defopts;
else
    opts = getoptions(inopts, defopts);
end
%% Global Variables
r = zeros(opts.n, 1);    %记录每个维度的局部搜索步长
g = zeros(opts.n, 1);    %记录每个维度随机生成点的个数
x = zeros(opts.n, opts.m + opts.v);            %记录所有蚂蚁的坐标位置
x_best_so_far = zeros(opts.n, 1);    %记录目前最优蚂蚁的位置
x_best_val = 0;                 %记录目前最优蚂蚁的值
x_temp = zeros(opts.n, 1);           %记录蚂蚁位置的临时变量
tau = zeros(opts.n, opts.m + opts.v);          %记录信息素的分布
tau_best_so_far = zeros(opts.n, 1);  %记录最优点信息素在各个维度上的分布
Fx = zeros(1, opts.m + opts.v);           %记录所有蚂蚁的目标函数值
q = 0;              %记录随机数
l = zeros(opts.n, opts.m);    %构造每个蚂蚁的解的时候的记录标尺

%% Static Variables
FES = 0;  %计算目标函数的次数
bestever = [];     %每次迭代时的最优解

%% Initialization

r = (Upper - Lower) / (2 * opts.m) * ones(opts.n, 1);
g = (opts.m + opts.v) * ones(opts.n, 1);
tau = opts.T0 * ones(opts.n, opts.m + opts.v);
x = Lower + (Upper - Lower) * rand(opts.n, opts.m + opts.v);

Fx = F(x');
FES = FES + (opts.m + opts.v) * size(x, 2);

% 按照目标函数值对所有蚂蚁进行排序
temp = sortrows([Fx; x; tau]')';
x = temp([2 : (opts.n + 1)], :);
x_best_val = temp(1, 1);
x_best_so_far = x(:, 1);    %记录最优解的点
tau = temp([(opts.n + 2) : (2 * opts.n + 1)], :);
%排序完毕

tau_best_so_far = opts.T0 * ones(opts.n, 1);

%% Generation of the candidate varible values
%Perform a dynamic exploitation process to the best-so-far solution
times = 0;
while (1)
    times = times + 1;
    g = zeros(opts.n, 1);
    flag = 0;
    for k = 1 : opts.v
        q = rand(opts.n, 1);
        
        i_1 = find(q < 1 / 3);
        x_temp(i_1, 1) = min([x_best_so_far(i_1, 1) + r(i_1, 1) ...
            .* rand(size(i_1)), Upper * ones(size(i_1))], [], 2);
        
        i_2 = find(q  >= 2 / 3);
        x_temp(i_2, 1) = max([x_best_so_far(i_2, 1) - r(i_2, 1) ...
            .* rand(size(i_2)), Lower * ones(size(i_2))], [], 2);
        
        i_3 = setdiff(1:opts.n, [i_1; i_2])';
        x_temp(i_3) = x_best_so_far(i_3, 1);
        
        i_4 = union(i_1', i_2')';
        x(i_4 + opts.n * (opts.m + g(i_4, 1) - 1)) = x_temp(i_4, 1);
        tau(i_4 + opts.n * (opts.m + g(i_4, 1) - 1)) = opts.T0;

        g(i_4, 1) = g(i_4, 1) + 1;
        Fx_temp = F(x_temp');
        FES = FES + 1 * size(x_temp, 2);
        if Fx_temp < x_best_val
            x_best_so_far = x_temp;
            x_best_val = Fx_temp;

            flag = 1;
        end
    end
    if flag ~= 0
        r = r * opts.Ve;
    else
        r = r * opts.Vr;
    end
    x(1:opts.n, opts.m - opts.theta + 1:opts.m) = Lower + (Upper - Lower) ...
        * rand(opts.n, opts.theta);
    tau(1:opts.n, opts.m - opts.theta + 1:opts.m) = opts.T0;

    %% Ant's solution construction
    for i = 1 : opts.n
        for k = 1 : opts.m
            q = rand;
            if q < opts.q0
                [~, l(i, k)] = max(tau(i, :));
            else
                %计算信息素的和
                %tau_sum = 0;
                tau_sum = sum(tau(i, 1:opts.m + g(i, 1)));
                tau_sum = tau_sum + tau_best_so_far(i, 1);
                %计算各点信息素的百分比
                p_best = tau_best_so_far(i, 1) / tau_sum;
                temp = p_best;
                p_other = zeros(opts.m + g(i, 1), 1);
                for j = 1 : opts.m + g(i, 1)
                    temp = temp + tau(i, j);
                    p_other(j, 1) = temp / tau_sum;
                end

                %轮盘赌算法
                q_temp = rand;
                if q_temp < p_best
                    l(i, k) = 0;
                else
                    for j = 1 : opts.m + g(i, 1)
                        if q_temp < p_other(j, 1)
                            l(i, k) = j;
                            break;
                        end
                    end
                end
            end
            if l(i, k) == 0
                x_temp1(i, k) = x_best_so_far(i, 1);
                tau_temp(i, k) = tau_best_so_far(i, 1);
            else
                x_temp1(i, k) = x(i, l(i, k));
                tau_temp(i, k) = tau(i, l(i, k));
            end
        end
    end
    
    x(1:opts.n, 1:opts.m) = x_temp1(1:opts.n, 1:opts.m);
    tau(1:opts.n, 1:opts.m) = tau_temp(1:opts.n, 1:opts.m);
    Fx(1, 1:opts.m) = F(x(:, 1:opts.m)');
    FES = FES + opts.m;
    %% Pheromone update
    % 按照目标函数值对所有蚂蚁进行排序
    temp = sortrows([Fx; x; tau]')';
    x = temp([2 : (opts.n + 1)], :);
    x_best_val = temp(1, 1);
    x_best_so_far = x(:, 1);    %更新最优解
    tau = temp([(opts.n + 2) : (2 * opts.n + 1)], :);
    % 排序完毕

    tau(1:opts.n, 1:opts.m) = (1 - opts.rou) * tau(1:opts.n, 1:opts.m) ...
        + opts.rou * opts.T_min;            %信息素挥发
    tau(1:opts.n, 1:opts.phi) = (1 - opts.alpha) * tau(1:opts.n, 1:opts.phi) ...
        + opts.alpha * opts.T_max;    %信息素加强
    
    tau_best_so_far = tau(:, 1);
    bestever = [bestever; x_best_val, FES];
    
    if (FES >= dimension * 10000 || (times > 1 && abs(bestever(times, 1) ...
            - bestever(times - 1, 1)) < opts.eps))
        break;
    end
end
xmin = x_best_so_far;
fmin = x_best_val;
end


function net = nn_net(NH, NIN, GB)
%NN_NET  Inverse-model network NIN-NH(1)-NH(2)-6 (tanh, tanh, linear) as a dlnetwork.
%
%   net = NN_NET()        6-64-64-6
%   net = NN_NET(NH)      hidden widths NH, e.g. [32 32] or [16 16]
%   net = NN_NET(NH, 5)   5 inputs: x = [p/R; phi/90; theta/180 - 1]   (nn_train.m INP = 'ang')
%   net = NN_NET(Nh)      scalar NH: ONE hidden layer, 6-Nh-6 (reference/Force_model/Real_time_inverse model/
%                         forward_nn.pdf, Sec. 3):  I_hat = W2 tanh(W1 x + b1) + b2,  x = [r; pbar],
%                         Nh = 21 M 3 (630 for M = 10), 13 Nh + 6 parameters
%
%   net = NN_NET(Nh, 6, true)   biases b1, b2 also Glorot-uniform, U(-L, L), L = sqrt(6/(fan_in + fan_out))
%
%   Input  x = [r_hat(phi,theta); p/R]   6 x 1   (R set by the caller, nn_train.m; actuator frame)
%   Output I_hat_NN                      6 x 1   unit current, pole order P1..P6
%                                                (same order as rho_solve's I_opt)
%
%       h1 = tanh(W1 x  + b1)     W1 NH1 x 6,    b1 NH1 x 1
%       h2 = tanh(W2 h1 + b2)     W2 NH2 x NH1,  b2 NH2 x 1
%       I  =      W3 h2 + b3      W3 6 x NH2,    b3 6 x 1
%
%   The biases are kept on purpose: without them the network is odd,
%   I(-x) = -I(x), and since the force is quadratic in I, +r_hat and -r_hat
%   would give the same force at p = 0.
%   Weights start Glorot-uniform, biases at zero (the fullyConnectedLayer
%   defaults). No input normalisation: x is already O(1) by construction.
%
%   Requires Deep Learning Toolbox.

    if nargin < 1 || isempty(NH), NH = [64 64]; end
    if nargin < 2 || isempty(NIN), NIN = 6; end
    if nargin < 3, GB = false; end
    NOUT = 6;
    if isscalar(NH)                                            % one hidden layer: W1, b1 (fc1), W2, b2 (I)
        bi = {};  bo = {};
        if GB
            bi = {'BiasInitializer', @(sz) gb_(sz, NIN, NH)};  bo = {'BiasInitializer', @(sz) gb_(sz, NH, NOUT)};
        end
        layers = [
            featureInputLayer(NIN, 'Normalization','none', 'Name','x')
            fullyConnectedLayer(NH, 'Name','fc1', bi{:})
            tanhLayer('Name','tanh1')
            fullyConnectedLayer(NOUT, 'Name','I', bo{:})
        ];
    else
        layers = [
            featureInputLayer(NIN, 'Normalization','none', 'Name','x')
            fullyConnectedLayer(NH(1), 'Name','fc1')
            tanhLayer('Name','tanh1')
            fullyConnectedLayer(NH(2), 'Name','fc2')
            tanhLayer('Name','tanh2')
            fullyConnectedLayer(NOUT, 'Name','I')
        ];
    end
    net = dlnetwork(layers);
end

function b = gb_(sz, ni, no)
% Glorot-uniform bias: U(-L, L), L = sqrt(6 / (ni + no))
    b = sqrt(6 / (ni + no)) * (2*rand(sz) - 1);
end

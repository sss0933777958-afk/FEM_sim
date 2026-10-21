function net = nn_net(NH)
%NN_NET  Inverse-model network 6-NH(1)-NH(2)-6 (tanh, tanh, linear) as a dlnetwork.
%
%   net = NN_NET()        6-64-64-6
%   net = NN_NET(NH)      hidden widths NH, e.g. [32 32] or [16 16]
%
%   Input  x = [r_hat(phi,theta); p/R]   6 x 1   (R = 150 um, actuator frame)
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
    NIN = 6;  NOUT = 6;
    layers = [
        featureInputLayer(NIN, 'Normalization','none', 'Name','x')
        fullyConnectedLayer(NH(1), 'Name','fc1')
        tanhLayer('Name','tanh1')
        fullyConnectedLayer(NH(2), 'Name','fc2')
        tanhLayer('Name','tanh2')
        fullyConnectedLayer(NOUT, 'Name','I')
    ];
    net = dlnetwork(layers);
end

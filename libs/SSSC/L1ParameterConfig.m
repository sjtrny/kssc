%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
% Author: Xi PENG@milab.org, Sichuan University.
% pangsaai@gmail.com
% Date: 6, Apr. 2013
% Description:  This code is developed for parameter configurations. 

% IF your used any part of this code, PLEASE approximately cited our works;

% Reference:
% [1] Xi Peng, Lei Zhang, Zhang Yi,
%     Scalable Sparse Subspace Clustering,
%     The 26th IEEE Conference on Computer Vision and Pattern Recognition (CVPR), Portland, Oregon, USA, June, 2013.
% [2] Xi Peng, Lei Zhang, Zhang Yi, 
%     Constructing L2-Graph for Subspace Learning and Segmentation, 
%     arXiv:1209.0841.

% All rights reserved.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
% EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL HOLDER AND CONTRIBUTORS BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/

function par = L1ParameterConfig(par)
% --------- data configuration
par.nDim               =   114; % the eigenfaces dimension
par.landmarkNO         =   228; % the number of in-sample data
% --- optimization algorithms configuration
par.lambda             =   [1e-7];
par.tolerance          =   [0.001];
par.maxIteration = 5000;
par.isNonnegative = false;
STOPPING_GROUND_TRUTH = -1;
STOPPING_DUALITY_GAP = 1;
STOPPING_SPARSE_SUPPORT = 2;
STOPPING_OBJECTIVE_VALUE = 3;
STOPPING_SUBGRADIENT = 4;
par.stoppingCriterion = STOPPING_OBJECTIVE_VALUE;



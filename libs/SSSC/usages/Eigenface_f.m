%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
% Author: Xi PENG@milab.org, Sichuan University.
% pangsaai@gmail.com
% Date: 6, Apr. 2013
% Description:  This code is developed for reducing the dimensionality of
% data based on the code of CRC [3].

% IF your used any part of this code, PLEASE approximately cited the following works;

% Reference:
% [1] Xi Peng, Lei Zhang, Zhang Yi,
%     Scalable Sparse Subspace Clustering,
%     The 26th IEEE Conference on Computer Vision and Pattern Recognition (CVPR), Portland, Oregon, USA, June, 2013.
% [2] Xi Peng, Lei Zhang, Zhang Yi, 
%     Constructing L2-Graph for Subspace Learning and Segmentation, 
%     arXiv:1209.0841.
% [3] L. Zhang, M. Yang, and X. Feng. 
%     Sparse representation or collaborative representation: Which helps face recognition? 
%     In IEEE International Conference on Computer Vision,2011.


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

function [disc_set,disc_value,Mean_Image]=Eigenface_f(Train_SET,Eigen_NUM)
%------------------------------------------------------------------------
% Eigenface computing function

[NN,Train_NUM]=size(Train_SET);

if NN<=Train_NUM % not small sample size case
    
   Mean_Image=mean(Train_SET,2);  
   Train_SET=Train_SET-Mean_Image*ones(1,Train_NUM);
   R=Train_SET*Train_SET'/(Train_NUM-1);
   
   [V,S]=Find_K_Max_Eigen(R,Eigen_NUM);
   disc_value=S;
   disc_set=V;

else % for small sample size case
    
   Mean_Image=mean(Train_SET,2);  
   Train_SET=Train_SET-Mean_Image*ones(1,Train_NUM);

  R=Train_SET'*Train_SET/(Train_NUM-1);
  
  [V,S]=Find_K_Max_Eigen(R,Eigen_NUM);
  disc_value=S;
  disc_set=zeros(NN,Eigen_NUM);
  clear R S;
  Train_SET=Train_SET/sqrt(Train_NUM-1);
  
  for k=1:Eigen_NUM
    a = Train_SET*V(:,k);
    b = (1/sqrt(disc_value(k)));
    disc_set(:,k)=b*a;
  end

end

function [Eigen_Vector,Eigen_Value]=Find_K_Max_Eigen(Matrix,Eigen_NUM)

[NN,NN]=size(Matrix);
[V,S]=eig(Matrix); %Note this is equivalent to; [V,S]=eig(St,SL); also equivalent to [V,S]=eig(Sn,St); %

S=diag(S);
[S,index]=sort(S);

Eigen_Vector=zeros(NN,Eigen_NUM);
Eigen_Value=zeros(1,Eigen_NUM);

p=NN;
for t=1:Eigen_NUM
    Eigen_Vector(:,t)=V(:,index(p));
    Eigen_Value(t)=S(p);
    p=p-1;
end
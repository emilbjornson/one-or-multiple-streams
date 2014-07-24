function [sumRate,W,Wnorm] = functionBlockDiagonalization(H,N,M,K,P,E)
%Compute block-diagonalization precoding (or zero-forcing for after a fixed
%receive combining has been performed). This is used in the article:
%
%Emil Björnson, Marios Kountouris, Mats Bengtsson, Björn Ottersten,
%“Receive Combining vs. Multi-Stream Multiplexing in Downlink Systems with
%Multi-Antenna Users,” IEEE Transactions on Signal Processing, vol. 61, no.
%13, pp. 3431-3446, July 2013.
%
%Download article: http://arxiv.org/pdf/1207.2776.pdf
%
%This is version 1.0 (Last edited: 2014-07-24)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
%
%
%The implementation is based on the original implementation of block-
%diagonalization proposed for "Sum Capacity Block Diagonalization" in:
%
%Q. Spencer, A. Swindlehurst, and M. Haardt, “Zero-forcing methods
%for downlink spatial multiplexing in multiuser MIMO channels,” IEEE
%Trans. Signal Process., vol. 52, no. 2, pp. 461–471, 2004.
%
%INPUT:
%H       = MK x N channel matrix
%N       = Number of transmit antennas
%M       = Number of receive antenans per user
%K       = Number of users
%P       = Vector with different total transmit powers
%E       = Optional parameter that determine the channel uncertainty. It is
%          either a K x K diagonal matrix with the error covariance at each
%          of the users or left out under perfect channel knowledge
%
%OUTPUT:
%sumRate = Vector with the maximum sum rates under block-diagonalization
%          with different total transmit powers
%W       = N x MK precoding matrix 
%Wnorm   = N x MK precoding matrix without power allocation


%Initialize the normalized precoding matrix
Wnorm = zeros(N,M*K);

%Initialize the precoding matrix for different total transmit powers
W = zeros(N,M*K,length(P));

%Set the error covariance matrix to an empty matrix (perfect CSI) if it was
%not given as input.
if nargin < 6
    E = [];
end

%Placeholder for the sum rate with different total transmit powers
sumRate = zeros(1,length(P));

%Placeholder for the eigenvalues of the precoded channels
lambda = zeros(M*K,1);


%Go through all the users
for j = 1:K
    
    Hj = H((j-1)*M+1:j*M,:); %Notation for the channel to User j
    
    %Compute right null space of the combined channel to the other users
    if K > 1
        Htilde = H([1:(j-1)*M j*M+1:K*M],:);
        Vtilde0 = null(Htilde); %Compute null space
    else
        Vtilde0 = eye(N); %Null space is the whole space if K=1
    end
    
    %Compute the singular value decomposition of the part of the channel to
    %user j that lies in the null space of the other users' channels
    [~,S,V] = svd(Hj*Vtilde0,'econ');
    
    %Store the eigenvalues (i.e., the squared singular values), which are
    %used for computing power allocation below
    lambda((j-1)*M+1:j*M,1) = diag(S).^2;
    
    %Store the normalized directions of the precoding matrix
    Wnorm(:,(j-1)*M+1:j*M) = Vtilde0*V;
end


%Go through all total transmit powers
for powerIndex = 1:length(P)

    lambdaInv = sort(1./(lambda+1e-8)); %Invert and sort the eigenvalues
    alpha_candidates = (P(powerIndex)+cumsum(lambdaInv))./(1:M*K)'; %Compute different values on the Lagrange multiplier (i.e., waterlevel) given that 1,2,...,MK of the singular directions get non-zero power
    optimalIndex = alpha_candidates-lambdaInv(1:end,1)>0 & alpha_candidates-[lambdaInv(2:end,1); Inf]<0; %Find the true Lagrange multiplier alpha by checking which one of the candidates that only turns on the singular directions that are supposed to be on
    waterLevel = alpha_candidates(optimalIndex); %Extract the optimal Lagrange multiplier (i.e., waterlevel in the waterfilling analogy)
    
    powerAllocation = waterLevel-1./lambda; %Compute power allocation
    powerAllocation(powerAllocation<0) = 0; %Make sure that inactive singular directions receive zero power
    
    %Compute precoding with optimal power allocation
    W(:,:,powerIndex) = Wnorm*diag(sqrt(powerAllocation));
    
    %Compute the sum rate
    if isempty(E)
        
        %Perfect channel knowledge
        sumRate(powerIndex) = real(log2(det(eye(M*K)+(H*(W(:,:,powerIndex)*W(:,:,powerIndex)')*H'))));
        
    else
        
        %Imperfect channel knowledge, characterized by the parameter E
        sumRate(powerIndex) = real(log2(det(eye(M*K)+(P(powerIndex)*(K-1)/K)*E+(H*(W(:,:,powerIndex)*W(:,:,powerIndex)')*H'))))-real(log2(det(eye(M*K)+(P(powerIndex)*(K-1)/K)*E)));
    
    end
    
end

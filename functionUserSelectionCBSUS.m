function [Omega,HH] = functionUserSelectionCBSUS(H,N,M,K,P,E)
%This is an implementation of the Capacity-based suboptimal user selection
%(CBSUS) algorithm, which was proposed in the paper:
%
%Z. Shen, R. Chen, J. Andrews, R. Heath, and B. Evans, “Low complexity
%user selection algorithms for multiuser MIMO systems with block
%diagonalization,” IEEE Trans. Signal Process., vol. 54, no. 9, pp. 3658–
%3663, 2006.
%
%This is a greedy user selection algorithm that selects users sequentially
%to maximize the sum rate with block-diagonalization. The algorithm will in
%general select fewer users than is maximally supported. The original
%algorithm is modified to 1) look for scheduling sets of any size (i.e.,
%not stopping when the sum rate is not increasing anymore) and 2) support
%imperfect CSI at the transmitter side, where the channel uncertainty is
%treated as additive noise. This is used in the article:
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
%INPUT:
%H       = MK x N channel matrix
%N       = Number of transmit antennas
%M       = Number of receive antenans per user
%K       = Number of users
%P       = Total transmit power
%E       = Optional parameter that determine the channel uncertainty. It is
%          either a K x K diagonal matrix with the error covariance at each
%          of the users or left out under perfect channel knowledge
%
%OUTPUT:
%Omega   = Scheduling set with indices of the selected users
%HH      = Channel matrix for the selected users


%Initate the scheduling set
Omega = [];

%Initiate set of unselected users
Upsilon = 1:K;

%Initiate channel matrix of selected users
HH = [];

%Maximal sum rate that has been achieved so far
Ctemp = 0;

%Set the error covariance matrix to an empty matrix (perfect CSI) if it was
%not given as input.
if nargin < 6
    E = [];
end
EE = []; %Initiate the block diagonal matrix of the error covariance matrix for selected users

%Placeholder for the sum rates with different numbers of selected users
maxSumRates = zeros(floor(N/M),1);


for nbrOfActiveUsers = 1:floor(N/M)
    
    %Placeholder for the sum rates when each of unselected users are added
    sumRatesCandidates = zeros(length(Upsilon),1);
    
    %Go through all of the unselected users
    for index = 1:length(Upsilon)
        
        k = Upsilon(index); %Candidate user
        
        Htemp = [HH; H((k-1)*M+1:k*M,:)]; %Compute channel of selected users if the candidate user is added
        
        %Compute the sum rate if User k is selected
        if isempty(E) == true
            
            %Perfect CSI as considered in "Low complexity user selection
            %algorithms for multiuser MIMO systems with block
            %diagonalization" by Z. Shen et al.
            sumRatesCandidates(index) = functionBlockDiagonalization(Htemp,N,M,nbrOfActiveUsers,P); %Compute sum rate when candidate user is added
            
        else
            
            %Imperfect CSI as proposed in the article, see footnotes 12 and 13.
            sumRatesCandidates(index) = functionBlockDiagonalization(Htemp,N,M,nbrOfActiveUsers,P,blkdiag(EE,E((k-1)*M+1:k*M,(k-1)*M+1:k*M)));
            
        end
        
    end
    
    %Find the candidate user that maximizes the sum rate
    [maxSumRates(nbrOfActiveUsers),userIndex] = max(sumRatesCandidates);
    
    
    
    Omega = [Omega; Upsilon(userIndex)]; %Add user to scheduling set
    HH = [HH; H((Upsilon(userIndex)-1)*M+1:Upsilon(userIndex)*M,:)]; %Add user's channel to the channel matrix of selected users
    
    %Update covariance matrix of channel uncertainty
    if isempty(E)==0
        EE = blkdiag(EE,E((Upsilon(userIndex)-1)*M+1:Upsilon(userIndex)*M,(Upsilon(userIndex)-1)*M+1:Upsilon(userIndex)*M));
    end
    
    Upsilon = Upsilon([1:userIndex-1 userIndex+1:length(Upsilon)]); %Remove selected user from set of unselected users
    
end

%Determine the number of selected users that maximizes the sum rate
[~,numberOfSelected] = max(maxSumRates);

%Prepare scheduling set for output
Omega = Omega(1:numberOfSelected);

%Prepare channel matrix of selected users for output
HH = HH(1:numberOfSelected*M,:);

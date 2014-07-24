function sumRate = functionSumrateComputation(H,W,M,K)
%Compute sum rates achieved by different precoding matrices. This is used
%in the article:
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
%INPUT:
%H       = MK x N channel matrix
%W       = N x MK x L matrix with L precoding matrices
%M       = Number of receive antennas per user
%K       = Number of users
%
%OUTPUT:
%sumRate = L x 1 vector with the sum rates achieved by each precoding matrix


%Extract number of precoding matrices
numberOfPrecodingMatrices = size(W,3);

%Placeholder for the user rates
rates = zeros(numberOfPrecodingMatrices,K);

%Go through all precoding matrices
for precodingIndex = 1:numberOfPrecodingMatrices
    
    %Go through all users
    for k = 1:K
        
        Hk = H((k-1)*M+1:k*M,:); %Channel to User k

        G = Hk*W(:,:,precodingIndex); %Effective channels for all signals
        G2 = Hk*W(:,(k-1)*M+1:k*M,precodingIndex); %Effective channels for useful signals
        
        %Compute achievable information rate based on Eq. (4)
        rates(precodingIndex,k) = real(log2(det(eye(M)+G*G'))-log2(det(eye(M)+G*G'-G2*G2')));
    
    end
    
end

%Compute achievable sum rates
sumRate = sum(rates,2);

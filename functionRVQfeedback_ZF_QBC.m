function [Hquantized,Heffective] = functionRVQfeedback_ZF_QBC(H,codebookSize,N,M,K)
%This is an implementation of the channel quantization based on random
%vector quantization (RVQ) and the quantization-based combining (QBC)
%procedure, which was proposed in the paper:
%
%N. Jindal, “Antenna combining for the MIMO downlink channel,” IEEE Trans.
%Wireless Commun., vol. 7, no. 10, pp. 3834–3844, 2008.
%
%QBC is a procedure to select the codeword and receive combining vector
%jointly to minimize the quantization error. This is used in the article:
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
%H            = MK x N channel matrix
%codebookSize = Number of entries in the codebook
%N            = Number of transmit antennas
%M            = Number of receive antenans per user
%K            = Number of users
%
%OUTPUT:
%Hquantized   = K x N matrix where the k:th row is quantized and normalized
%               channel to User k
%Heffective   = K x N matrix where the k:th row is the true channel after
%               the receive combining has been applied 

%Placeholder for quantized channel
Hquantized = zeros(K,N);

%Placeholder for true effective channel after receive combining
Heffective = zeros(K,N);

%Go through all users
for k = 1:K
    
    %Generate a random codebook
    RVQcodebook = randn(N,codebookSize)+1i*randn(N,codebookSize); %Generate random directions
    RVQcodebook = RVQcodebook./repmat(sqrt(sum(abs(RVQcodebook).^2,1)),[N 1]); %Normalize the random directions
    
    Hk = H((k-1)*M+1:k*M,:); %Channel of User k
    
    Q = orth(Hk'); %Compute an orthogonal basis of the column span of the channel (Step 1 in Section IV.B of [8])
    
    [~,index] = max(sum(abs(Q'*RVQcodebook).^2,1)); %Find the codebook vector closest to the channel subspace (Step 2 in Section IV.B of [8])
    
    codeword = RVQcodebook(:,index); %Extract the closest codeword
    
    sproj = (Q*Q')*codeword; %Determine the direction of the effective channel (Step 3 in Section IV.B of [8])
    
    %Compute the receive combining vector (Step 4 in Section IV.B of [8])
    gammak = Hk'\sproj;
    gammak = gammak/norm(gammak);

    Hquantized(k,:) = codeword'; %Store the quantized channel direction
    
    Heffective(k,:) = gammak'*Hk; %Store the true effective channel after receive combining
    
end

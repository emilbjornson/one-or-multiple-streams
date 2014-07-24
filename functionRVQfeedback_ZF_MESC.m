function [Hquantized,Heffective] = functionRVQfeedback_ZF_MESC(H,codebookSize,N,M,K,P)
%This is an implementation of the channel quantization based on random
%vector quantization (RVQ) and the maximum expected SINR combiner (MESC),
%which was proposed in the paper:
%
%M. Trivellato, F. Boccardi, and H. Huang, “On transceiver design and
%channel quantization for downlink multiuser MIMO systems with limited
%feedback,” IEEE J. Sel. Areas Commun., vol. 26, no. 8, pp. 1494–1504, 2008.
%
%MESC is a procedure to select the codeword and receive combining vector
%jointly to maximize an SINR-like metric. This is used in the article:
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
%P            = Total transmit power
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

%Average transmit power per stream
Peffective = P/N;

%Go through all users
for k = 1:K
    
    %Generate a random codebook
    RVQcodebook = randn(N,codebookSize)+1i*randn(N,codebookSize); %Generate random directions
    RVQcodebook = RVQcodebook./repmat(sqrt(sum(abs(RVQcodebook).^2,1)),[N 1]); %Normalize the random directions
    
    Hk = H((k-1)*M+1:k*M,:); %Channel of User k
    
    %Compute the best codeword based on the MESC metric in Eq. (11) in [9].
    gamma_ki = zeros(codebookSize,1);
    
    for i = 1:codebookSize %Go through all codwords
        
        Ak = Peffective*(Hk*(RVQcodebook(:,i)*RVQcodebook(:,i)')*Hk'); %Defined in Eq. (13) in [9]
        Bk = eye(M)+Peffective*Hk*(eye(N)-RVQcodebook(:,i)*RVQcodebook(:,i)')*Hk'; %Defined in Eq. (14) in [9]
        gamma_ki(i) = real(max(eig(Ak,Bk))); %Compute the value of the metric in Eq. (12) in [9].
        
    end
    
    [~,index] = max(gamma_ki); %Find the codebook vector that maximizes the metric
    
    %Compute the B-matrix from Eq. (14) in [9] for the selected codeword
    Bk = eye(M)+Peffective*Hk*(eye(N)-RVQcodebook(:,index)*RVQcodebook(:,index)')*Hk';
    
    %Compute the corresponding receive combining vector based on Eq. (15)
    %in [9]
    u = Bk\(Hk*RVQcodebook(:,index));
    u = u/norm(u);
    
    Hquantized(k,:) = RVQcodebook(:,index)'; %Store the quantized channel direction
    
    Heffective(k,:) = u'*Hk; %Store the true effective channel after receive combining
    
end

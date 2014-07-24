function Hquantized = functionRVQfeedback_BD(H,codebooksize,N,M,K)
%This is an implementation of the channel quantization of N x M channels,
%based on the principles of random vector quantization (RVQ). The same
%metrics are used as in the paper:
%
%N. Ravindran and N. Jindal, “Limited feedback-based block diagonalization
%for the MIMO broadcast channel,” IEEE J. Sel. Areas Commun., vol. 26, no.
%8, pp. 1473–1482, 2008.
%
%This is used in the article:
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
%codebookSize = Number of entries in the codebook per channel dimension
%N            = Number of transmit antennas
%M            = Number of receive antenans per user
%K            = Number of users
%
%OUTPUT:
%Hquantized   = MK x N matrix where each block of M rows are the quantized
%               normalized subspace spanned by the channel of User k


%Placeholder for quantized channel
Hquantized = zeros(N,N);

%Go through all users
for k = 1:K
    
    %Generate a random codebook
    RVQcodebook = zeros(N,M*codebooksize); %Placeholder for codebook

    randomRealizations = randn(N,M*codebooksize)+1i*randn(N,M*codebooksize); %Random realizations

    for n = 1:codebooksize
        RVQcodebook(:,(n-1)*M+1:n*M,:) = orth(randomRealizations(:,(n-1)*M+1:n*M,:)); %Compute an orthogonal basis for the rank-M subspace spanned by the random realization
    end
    
    Hk = H((k-1)*M+1:k*M,:); %Channel of User k
    
    %Select the codeword that minimizes the chordal distance to the true
    %channel, as proposed in Eq. (2) in [7].
    Q = orth(Hk');
    [~,ind] = max(sum(reshape(sum(abs(Q'*RVQcodebook).^2,1),[M codebooksize]),1));
    
    %Store the quantized channel subspace
    Hquantized((k-1)*M+1:k*M,:) = RVQcodebook(:,(ind-1)*M+1:ind*M)';
    
end

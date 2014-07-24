function activeStreams = functionGreedyStreamAllocation(H,N,M,K,P,E,maxStreamsPerUser)
%Compute the achievable sum rate by a greedy data stream allocation, where
%the streams are allocated sequentially to the user that would make the
%largest improvement in the sum rate. The precoding is based on
%block-diagonalization, with different dimensions of the block. The
%implementation resembles the multiuser eigenmode transmission (MET)
%proposed in the paper:
%
%F. Boccardi and H. Huang, “A near-optimum technique using linear
%precoding for the MIMO broadcast channel,” in Proc. IEEE ICASSP,
%vol. 3, 2007, pp. 17–20.
%
%The implemention goes beyond the MET scheme in an attempt to find an even
%better solution. The main difference is that we continue to run the
%algorithm until all N streams have been allocated, to avoid local optima
%where one extra stream cannot improve the sum rate but multiple extra
%streams might give an improvement. Moreover, we have extended the
%algorithm to take imperfect CSI into account and to, consequently, have a
%maximal number of streams per user (i.e., number of eigenmodes that have
%been estimated). This implemenation is used in the paper:
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
%N       = Number of transmit antennas and maximal number of streams
%M       = Number of receive antenans per user
%K       = Number of users
%P       = Vector with different total transmit powers
%E       = Optional parameter that determine the channel uncertainty. It is
%          either a K x K diagonal matrix with the error covariance at each
%          of the users or left out under perfect channel knowledge
%maxStreamsPerUser = Optional K x 1 vector with the maximal number of data
%streams that each user can be allocated.
%
%OUTPUT:
%activeStreams = KM x 1 vector with indices of the active data streams
%                selected by the algorithm


%Set the error covariance matrix to an empty matrix (perfect CSI) if it was
%not given as input. Otherwise, the diagonal elements are extracted for
%later use.
if nargin<6
    E = [];
else
    E = diag(E);
end

%The maximal number of streams allocated per user is either given as an
%input parameter or set to its maximal value.
if nargin<7
    maxStreamsPerUser = M*ones(K,1);
end


%Matrix over all users and their maximal number of data streams. Each
%column represents the stream allocation for the corresponding of active
%data streams
activeStreams = zeros(K*M,N);

%Preallocate the vector with the achievable sum rate for each number of
%active streams.
achievableSumRate = zeros(N,1);

%Allocate one stream at a time, up to the maximal number of N
for streamNumber = 1:N
    
    %Prepare to store the sum rate that is achieved when User j is
    %allocated the next stream
    performancePerUser = zeros(K,1);
    
    %Prepare to store the index of the stream that the User j wants to
    %have allocated.
    streamIndex = zeros(K,1);
    
    %Go through all users
    for j = 1:K
        
        %Extract the indices of the streams that User j has already been allocated
        activeStreamUserj = activeStreams(((j-1)*M+1:j*M),streamNumber);
        
        %Check if the user is allowed to have one more stream
        if sum(activeStreamUserj) < maxStreamsPerUser(j)
            
            %Compute the indices of the streams that the user might be
            %allocated in this iteration
            candidate_streams = (1:maxStreamsPerUser(j));
            candidate_streams = candidate_streams(activeStreamUserj(1:maxStreamsPerUser(j))==0);
            
            %Prepare to store the sum rate that each user achieves when
            %they is allocated each of potential streams
            performanceForDifferentStreams = zeros(maxStreamsPerUser(j),1);
            
            %Go through the candidate streams and compute the achievable
            %sum rate when this one is allocated
            for m = candidate_streams
                
                %Update the vector with stream allocations
                activeStreamsCandidate = activeStreams(:,streamNumber);
                activeStreamsCandidate((j-1)*M+m) = 1;
                
                %Compute the achievable performance
                performanceForDifferentStreams(m) = functionGreedyZeroForcing(H,N,M,K,P,activeStreamsCandidate,E);
                
            end
            
            %Store the highest achievable performance among the different
            %streams that User j can be allocated, and the corresponding index
            [performancePerUser(j),streamIndex(j)] = max(performanceForDifferentStreams);
        end

    end
    
    %Determine which user that gets the next stream, based on that it gives
    %the largest sum rate in this stage
    [achievableSumRate(streamNumber),ind] = max(performancePerUser);
    
    %Update the list of active streams
    activeStreams((ind-1)*M+streamIndex(ind),streamNumber) = 1;

    %Prepare for the next iteration by initializing which streams that are
    %active from before
    if streamNumber<N
        activeStreams(:,streamNumber+1) = activeStreams(:,streamNumber);
    end
   
end

%Find the number of active streams that gives the highest sum rate
[~,optimalNumberOfStreams] = max(achievableSumRate);

%Prepare the output
activeStreams = activeStreams(:,optimalNumberOfStreams);

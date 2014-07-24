%This Matlab script can be used to generate Figure 4, 5, 6, 7, 8, and 9 in 
%the following article:
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
%Please note that the channels are generated randomly, thus the results
%will not be exactly the same as in the paper.

%Initialization
close all;
clear all;

%%Simulation parameters

rng('shuffle'); %Initiate the random number generators with a random seed
%%If rng('shuffle'); is not supported by your Matlab version, you can use
%%the following commands instead:
%randn('state',sum(100*clock));

N = 8; %Number of transmit antennas at base station
M = 4; %Number of receive antennas per user

%Range of number of users
Kvalues = [N:M:50 52];
Kmax = max(Kvalues); %Maximal number of users

%Range of spatial correlation factors
correlationFactors = [0 0.4 0.8];

%Number of realizations in the Monte Carlo simulations
nbrOfMonteCarloRealizations = 500;


%%Select which figure to simulate
simulateFigure = 9;

if simulateFigure == 4
    
    %Is perfect CSI available at the BS or only imperfect CSI from channel
    %estimation in the uplink?
    %true: Imperfect CSI
    %false: Perfect CSI
    imperfectChannelEstimation = false;
    
    %Will all users be located at the same distance (called the cell edge)
    %or the users distributed randomly in a cell?
    %0: All users at cell edge
    %1: Distribued uniformly in the cell
    %2: Distribued uniformly in the cell and priotize "strong" users
    variableLargeScaleFading = 0;
    
    %Range of SNR values (each specified at the cell edge)
    PdB = [10 20]; %dB scale
    P = 10.^(PdB/10); %Linear scale
    
elseif simulateFigure == 5 || simulateFigure == 6
    
    %Is perfect CSI available at the BS or only imperfect CSI from channel
    %estimation in the uplink?
    %true: Imperfect CSI
    %false: Perfect CSI
    imperfectChannelEstimation = false;
    
    %Will all users be located at the same distance (called the cell edge)
    %or the users distributed randomly in a cell?
    %0: All users at cell edge
    %1: Distribued uniformly in the cell
    %2: Distribued uniformly in the cell and priotize "strong" users
    variableLargeScaleFading = 1;
    
    %Store and plot the number of streams that MET allocate to different
    %users in the cell (as is done in Figure 6).
    storeStreamAllocation = true;
    
    %Range of SNR values (each specified at the cell edge)
    PdB = 20; %dB scale
    P = 10.^(PdB/10); %Linear scale
    
elseif simulateFigure == 7
    
    %Is perfect CSI available at the BS or only imperfect CSI from channel
    %estimation in the uplink?
    %true: Imperfect CSI
    %false: Perfect CSI
    imperfectChannelEstimation = true;
    
    %Will all users be located at the same distance (called the cell edge)
    %or the users distributed randomly in a cell?
    %0: All users at cell edge
    %1: Distribued uniformly in the cell
    %2: Distribued uniformly in the cell and priotize "strong" users
    variableLargeScaleFading = 0;
    
    %Range of SNR values (each specified at the cell edge)
    PdB = [10 20]; %dB scale
    P = 10.^(PdB/10); %Linear scale
    
elseif simulateFigure == 8
    
    %Is perfect CSI available at the BS or only imperfect CSI from channel
    %estimation in the uplink?
    %true: Imperfect CSI
    %false: Perfect CSI
    imperfectChannelEstimation = true;
    
    %Will all users be located at the same distance (called the cell edge)
    %or the users distributed randomly in a cell?
    %0: All users at cell edge
    %1: Distribued uniformly in the cell
    %2: Distribued uniformly in the cell and priotize "strong" users
    variableLargeScaleFading = 1;
    
    %Store and plot the number of streams that MET allocate to different
    %users in the cell (as is done in Figure 6).
    storeStreamAllocation = false;
    
    %Range of SNR values (each specified at the cell edge)
    PdB = 20; %dB scale
    P = 10.^(PdB/10); %Linear scale
    
elseif simulateFigure == 9
    
    %Is perfect CSI available at the BS or only imperfect CSI from channel
    %estimation in the uplink?
    %true: Imperfect CSI
    %false: Perfect CSI
    imperfectChannelEstimation = true;
    
    %Will all users be located at the same distance (called the cell edge)
    %or the users distributed randomly in a cell?
    %0: All users at cell edge
    %1: Distribued uniformly in the cell
    %2: Distribued uniformly in the cell and priotize "strong" users
    variableLargeScaleFading = 2;
    
    %Store and plot the number of streams that MET allocate to different
    %users in the cell (as is done in Figure 6).
    storeStreamAllocation = false;
    
    %Range of SNR values (each specified at the cell edge)
    PdB = 20; %dB scale
    P = 10.^(PdB/10); %Linear scale
    
end







%Define a circular cell with heterogeneous channel conditions
if variableLargeScaleFading > 0
    
    cellRadius = 250; %Cell radius (in meters)
    minimalUserDistance = 35; %Minimal distance between a user and the BS (in meters)
    pathlossCoefficient = 3.5; %Exponent in the distance dependent pathloss model
    
    shadowFadingStddB = 8; %Standard deviation of the shadow fading (in dB scale)
    shadowFadingStd = 10.^(shadowFadingStddB/10); %Standard deviation of the shadow fading (in linear scale)
    
    
    %Prepare to store the number of streams that are allocated in different
    %parts of the cell
    if storeStreamAllocation == true
        
        cellCenterMaxDistance = 100; %Users at a distance <= 100 meters from the BS is considered as cell center users
        cellEdgeMinDistancce = 200; %Users at a distance >= 200 meters from the BS is considered as cell edge users
        
        %Placeholders for storing the number of streams allocated per user
        %with the greedy MET algorithm
        streamAllocation_Average = zeros(M,nbrOfMonteCarloRealizations,length(P),length(Kvalues),length(correlationFactors));
        streamAllocation_CellCenter = zeros(M,nbrOfMonteCarloRealizations,length(P),length(Kvalues),length(correlationFactors));
        streamAllocation_CellEdge = zeros(M,nbrOfMonteCarloRealizations,length(P),length(Kvalues),length(correlationFactors));
        
    end
    
end


%Compute the gain ||h||^2 from Lemma 1 by Monte-Carlo simulations
if imperfectChannelEstimation == true
    
    meanh2 = zeros(M,length(correlationFactors));
    
    nbrOfMonteCarloRealizations_channelgain = 10000; %Number of realizations
    
    Huncorrelated_channelgain = (randn(M*nbrOfMonteCarloRealizations_channelgain,N)+1i*randn(M*nbrOfMonteCarloRealizations_channelgain,N))/sqrt(2);
    
    for m = 1:length(correlationFactors)
        
        Rrx = toeplitz( correlationFactors(m).^(0:M-1) );
        
        for k = 1:nbrOfMonteCarloRealizations_channelgain
            Hcorrelated = sqrtm(Rrx)*Huncorrelated_channelgain((k-1)*M+1:k*M,:);
            meanh2(:,m) = meanh2(:,m) + svd(Hcorrelated).^2/nbrOfMonteCarloRealizations_channelgain;
        end
        
    end
    
    %Remove the Monte-Carlo realizations to save memory
    clear Huncorrelated_channelgain;
    clear Hcorrelated;
    
end


%Placeholders for storing the simulation results
sumRateBD = zeros(nbrOfMonteCarloRealizations,length(P),length(Kvalues),length(correlationFactors));
sumRateZFC = zeros(nbrOfMonteCarloRealizations,length(P),length(Kvalues),length(correlationFactors));

if imperfectChannelEstimation == false || variableLargeScaleFading == 2
    sumRateMET = zeros(nbrOfMonteCarloRealizations,length(P),length(Kvalues),length(correlationFactors));
end

%Go through all Monte Carlo realizations
for n = 1:nbrOfMonteCarloRealizations
    
    %Output of the simulation progress
    disp(['Progress: ' num2str(n) ' out of ' num2str(nbrOfMonteCarloRealizations) ' realizations.']);
    
    %Generate channel matrices for current Monte Carlo realizations
    Huncorrelated = (randn(M*Kmax,N)+1i*randn(M*Kmax,N))/sqrt(2); %Channels for the maximal number of users (before spatial correlation has been applied)
    Nmax = (randn(M*Kmax,N)+1i*randn(M*Kmax,N))/sqrt(2); %Noise matrix for the maximal number of users
    thetaRx = 2*pi*rand(Kmax,1);  %Angle to each user as seen from the base station
    
    
    %Initialize homogenous or heterogenous channel variances
    if variableLargeScaleFading==1 || variableLargeScaleFading==2
        
        shadowFadingRealizations=10.^(randn(Kmax,1)*shadowFadingStd/10); %Generate shadow fading realizations
        relativeDistances = sqrt(rand(Kmax,1)*(1^2-(minimalUserDistance/cellRadius)^2)+(minimalUserDistance/cellRadius)^2); %Generate random user locations with a minimal distance
        channelVariances = 1./(relativeDistances.^pathlossCoefficient)./shadowFadingRealizations; %Computing the channel variances with distance-dependent pathloss and shadow fading
        channelVariances = channelVariances/exp((log(10)*shadowFadingStd/10)^2/2); %Normalize the shadow fading such that the channel variance is one at the cell edge
        
    elseif variableLargeScaleFading == 0
        
        channelVariances = ones(Kmax,1); %Set all channel variances to one in the case with heterogenous channel conditions
        
    end
    
    
    %Go through different number of users
    for l = 1:length(Kvalues)
        
        %Extract out the current total number of users (from which a subset
        %will be selected for transmission).
        Ktotal = Kvalues(l);
        
        if variableLargeScaleFading == 2
            effectiveDistances = sort(channelVariances(1:Ktotal),'descend');
        elseif variableLargeScaleFading == 0 || variableLargeScaleFading == 1
            effectiveDistances = channelVariances(1:Ktotal);
        end
        
        
        %Go through all spatial correlation factors
        for m = 1:length(correlationFactors)
            
            %Extract out the spatial correlation factor at the receiver side
            rRx = correlationFactors(m);
            
            %Compute correlation matrix (based on exponential correlation)
            %and also extract the eigenvalues
            Rr = eye(M*Ktotal);
            eigenvalues = zeros(Ktotal*M,1);
            
            for k = 1:Ktotal
                Rr((k-1)*M+1:k*M,(k-1)*M+1:k*M) = effectiveDistances(k)*toeplitz( ( abs(rRx) * exp(1i*thetaRx(k)) ).^(0:M-1) ); %Exponential correlation model from Eq. (16).
                eigenvalues((k-1)*M+1:k*M,1) = eig(Rr((k-1)*M+1:k*M,(k-1)*M+1:k*M)); %Store eigenvalues of the correlation matrix
            end
            
            
            %Compute channel matrix with spatial correlation
            H = sqrtm(Rr)*Huncorrelated(1:M*Ktotal,:);
            
            %Extract the noise matrix for the current number of users
            Noise = Nmax(1:M*Ktotal,:);
            
            
            
            %Scheme: Block diagonalization (BD) with user selection
            %(4 streams per selected user).
            if imperfectChannelEstimation == false %Assuming perfect CSI
                
                %Go through all transmit powers
                for powerIndex = 1:length(P)
                    [scheduledBD,H_BD] = functionUserSelectionCBSUS(H,N,M,Ktotal,P(powerIndex)); %Select users based on the CBSUS algorithm
                    sumRateBD(n,powerIndex,l,m) = functionBlockDiagonalization(H_BD,N,M,length(scheduledBD),P(powerIndex));  %Compute sum rate for the selected users
                end
                
            elseif imperfectChannelEstimation == true %Assuming imperfect CSI, caused by channel estimation
                
                %Go through all transmit powers
                for powerIndex = 1:length(P)
                    
                    %Assuming that Ktotal is the number of channel
                    %directions that are estimated, Ktotal/M is the number
                    %of active users with BD.
                    KtotalBD = Ktotal/M;
                    
                    %Compute linear MMSE channel estimate of the channel
                    %and the corresponding error covariance matrix
                    Hestimate = (inv(Rr(1:KtotalBD*M,1:KtotalBD*M)) + P(powerIndex)*eye(KtotalBD*M))\(P(powerIndex)*H(1:KtotalBD*M,:)+sqrt(P(powerIndex))*Noise(1:KtotalBD*M,:));
                    Eest = inv(inv(Rr(1:KtotalBD*M,1:KtotalBD*M)) + P(powerIndex)*eye(KtotalBD*M));
                    
                    %Select users and compute approximate BD precoding
                    %using the estimated channel matrix
                    [scheduledBD,H_BD] = functionUserSelectionCBSUS(Hestimate,N,M,KtotalBD,P(powerIndex),Eest); %Select users based on the modified CBSUS algorithm (as described in Footnote 12)
                    [~,precodingBD] = functionBlockDiagonalization(H_BD,N,M,length(scheduledBD),P(powerIndex)); %Compute approximate BD precoding for the selected users
                    
                    %Extract out the true channels of the selected users
                    Hscheduled = zeros(length(scheduledBD)*M,N);
                    for j = 1:length(scheduledBD)
                        Hscheduled((j-1)*M+1:j*M,:) = H((scheduledBD(j)-1)*M+1:scheduledBD(j)*M,:);
                    end
                    
                    %Compute the achievable sum rate using the approximate
                    %BD precoding matrix
                    sumRateBD(n,powerIndex,l,m) = functionSumrateComputation(Hscheduled,precodingBD,M,length(scheduledBD));
                    
                end
                
            end
            
            
            
            %Scheme: Multi-user Eigenmode Transmission (MET), which gives a
            %greedy stream allocation
            
            if imperfectChannelEstimation == false %Assuming perfect CSI
                
                %Prepare to compute the effective channels with singular
                %value decomposition based precoding
                Heffective = zeros(M*Ktotal,N);
                
                %Go through all users
                for k = 1:Ktotal
                    Hk = H((k-1)*M+1:k*M,:); %Channel of User k
                    [U,S,V] = svd(Hk); %Compute singular value decomposition of the channel of User k
                    Heffective((k-1)*M+1:k*M,:) = S*V'; %Compute the effective orthogonalized channel when a singular value decomposition has been applied at the receiver side
                end
                
                %Go through all transmit powers
                for powerIndex = 1:length(P)
                    
                    %Determine which data streams that should be active
                    %using the MET scheme.
                    activeStreams = functionGreedyStreamAllocation(Heffective,N,M,Ktotal,P(powerIndex));
                    
                    %Compute the MET precoding for the active data streams.
                    [~,precoding,~] = functionGreedyZeroForcing(Heffective,N,M,Ktotal,P(powerIndex),activeStreams);
                    
                    %Go through all potential users
                    for k = 1:Ktotal
                        
                        nbrOfStreams = sum(activeStreams((k-1)*M+1:k*M)); %Compute the number of active data streams for the current users
                        
                        %Compute the achievable rate if the user has been
                        %allocated at least one data stream.
                        if nbrOfStreams>0
                            
                            Hk = H((k-1)*M+1:k*M,:); %Current user channel
                            
                            signal = Hk*(precoding(:,(k-1)*M+1:k*M)*precoding(:,(k-1)*M+1:k*M)')*Hk'; %Correlation matrix for the received signal
                            interference = eye(M)+Hk*(precoding*precoding')*Hk'-signal; %Correlation matrix for the received interference
                            
                            %Store the rate achieved by the user
                            sumRateMET(n,powerIndex,l,m) = sumRateMET(n,powerIndex,l,m) + real(log2(det(eye(M)+interference\signal)));
                            
                            
                            %Store the number of streams that each of the
                            %users are allocated with the MET scheme
                            if variableLargeScaleFading == 1 && storeStreamAllocation == true
                                
                                %Store the result for all users
                                streamAllocation_Average(nbrOfStreams,n,powerIndex,l,m) = streamAllocation_Average(nbrOfStreams,n,powerIndex,l,m)+1;
                                
                                %Separate the results for cell center and
                                %cell edge users
                                if relativeDistances(k) <= cellCenterMaxDistance/cellRadius
                                    streamAllocation_CellCenter(nbrOfStreams,n,powerIndex,l,m) = streamAllocation_CellCenter(nbrOfStreams,n,powerIndex,l,m)+1;
                                elseif relativeDistances(k) > cellEdgeMinDistancce/cellRadius
                                    streamAllocation_CellEdge(nbrOfStreams,n,powerIndex,l,m) = streamAllocation_CellEdge(nbrOfStreams,n,powerIndex,l,m)+1;
                                end
                                
                            end
                            
                        end
                        
                    end
                    
                end
                

            elseif imperfectChannelEstimation == true && variableLargeScaleFading == 2 %Assuming imperfect CSI, caused by channel estimation
                
                %Prepare to compute the effective channels with singular
                %value decomposition based precoding
                Heffective = zeros(Ktotal*M,N);
                
                %Go through all users
                for k = 1:Ktotal
                    Hk = H((k-1)*M+1:k*M,:); %Channel of User k
                    [U,S,V] = svd(Hk); %Compute singular value decomposition of the channel of User k
                    Heffective((k-1)*M+1:k*M,:) = S*V'; %Compute the effective orthogonalized channel when a singular value decomposition has been applied at the receiver side
                end
                
                
                %Set a threshold such that only the Ktotal statistically 
                %strongest eigendirection are estimated and used for
                %potential transmission
                sortedEigenvalues = sort(eigenvalues,1,'descend');
                threshold = sortedEigenvalues(Ktotal);
                
                %Use the threshold to compute the maximal number of streams
                %per user in the MET-like stream allocation scheme
                maxStreamsPerUser = zeros(Ktotal,1);
                for k = 1:Ktotal
                    maxStreamsPerUser(k) = sum(eigenvalues((k-1)*M+1:k*M,1)>=threshold);
                end
                
                %Go through all transmit powers
                for powerIndex = 1:length(P)
                    
                    %Compute linear MMSE channel estimate of the channel
                    %and the corresponding error covariance matrix
                    IoverLambdaMeanh2 = diag(repmat(1./meanh2(1:M,m)',[1 Ktotal]) ./ kron(effectiveDistances(1:Ktotal)',ones(1,M))); %Precompute 1/(eigenvalues*||h||^2)
                    Hestimate = (IoverLambdaMeanh2 + P(powerIndex)*eye(M*Ktotal))\(P(powerIndex)*Heffective+sqrt(P(powerIndex))*Noise(1:M*Ktotal,:));
                    Eest = inv(IoverLambdaMeanh2 + P(powerIndex)*eye(M*Ktotal));
                    
                    
                    %Compute the MET-like stream allocation and precoding
                    %by greedy stream allocation.
                    activeStreams = functionGreedyStreamAllocation(Hestimate,N,M,Ktotal,P(powerIndex),Eest,maxStreamsPerUser);
                    [~,precoding,~] = functionGreedyZeroForcing(Hestimate,N,M,Ktotal,P(powerIndex),activeStreams);
                    
                    for k = 1:Ktotal
                        
                        nbrOfStreams = sum(activeStreams((k-1)*M+1:k*M));  %Compute the number of active data streams for the current users
                        
                        %Compute the achievable rate if the user has been
                        %allocated at least one data stream.
                        if nbrOfStreams>0
                            
                            Hk = H((k-1)*M+1:k*M,:); %Current user channel
                            
                            signal = Hk*(precoding(:,(k-1)*M+1:k*M)*precoding(:,(k-1)*M+1:k*M)')*Hk'; %Correlation matrix for the received signal
                            interference = eye(M)+Hk*(precoding*precoding')*Hk'-signal; %Correlation matrix for the received interference
                            
                            %Store the rate achieved by the user
                            sumRateMET(n,powerIndex,l,m) = sumRateMET(n,powerIndex,l,m) + real(log2(det(eye(M)+interference\signal)));
                            
                        end

                    end
                    
                end
                
            end
            
            
            
            %Scheme: Zero-forcing with combining (ZFC) and user selection
            
            %Compute effective channel when the users apply receive
            %combining along the strong singular direction.
            Heffective = zeros(Ktotal,N);
            
            %Go through all users
            for k = 1:Ktotal
                Hk = H((k-1)*M+1:k*M,:); %Channel of User k
                [U,S,V] = svd(Hk,'econ'); %Compute singular value decomposition of the channel of User k
                Heffective(k,:) = U(:,1)'*Hk; %Compute and store the effective channel after receive combining
            end
            
            if imperfectChannelEstimation == false %Assuming perfect CSI
                
                %Go through all transmit powers
                for powerIndex = 1:length(P)
                    
                    %Select users based on the CBSUS algorithm
                    [scheduledZF,H_ZF] = functionUserSelectionCBSUS(Heffective,N,1,Ktotal,P(powerIndex));
                    
                    %Compute ZFC precoding for the selected users
                    [~,W_ZFC] = functionBlockDiagonalization(H_ZF,N,1,length(scheduledZF),P(powerIndex));
                    
                    %Placeholder for the effective channel after MMSE receive combining
                    HeffectiveMMSE = zeros(length(scheduledZF),N);
                    
                    %Go through the scheduled users
                    for k = 1:length(scheduledZF)
                        Hk = H((scheduledZF(k)-1)*M+1:scheduledZF(k)*M,:); %Channel of scheduled User k
                        
                        %Compute MMSE receive combiner as described in Remark 1
                        cMMSE = (eye(M)+Hk*(W_ZFC*W_ZFC')*Hk')\(Hk*W_ZFC(:,k));
                        
                        %Compute the effective channel with the MMSE receive combiner
                        if norm(cMMSE) > 0
                            cMMSE = cMMSE/norm(cMMSE);
                            HeffectiveMMSE(k,:) = cMMSE'*Hk;
                        end
                        
                    end
                    
                    %Compute sum rate for ZF with MMSE receive combining
                    channelGains = abs(HeffectiveMMSE*W_ZFC).^2;
                    signalPowers = diag(channelGains);
                    interferencePowers = sum(channelGains,2)-signalPowers;
                    sumRateZFC(n,powerIndex,l,m) = sum(log2(1+signalPowers./(1+interferencePowers)));
                    
                end
                
                
            elseif imperfectChannelEstimation == true %Assuming imperfect CSI, caused by channel estimation
                
                %Go through all transmit powers
                for powerIndex=1:length(P)
                    
                    %Compute linear MMSE channel estimate of the channel
                    %and the corresponding error covariance matrix
                    IoverLambdaMeanh2 = diag(1./effectiveDistances)/meanh2(1,m); %Precompute 1/(eigenvalues*||h||^2)
                    Hestimate = (IoverLambdaMeanh2 + P(powerIndex)*eye(Ktotal))\(P(powerIndex)*Heffective(1:Ktotal,:)+sqrt(P(powerIndex))*Noise(1:Ktotal,:));
                    Eest = inv(IoverLambdaMeanh2 + P(powerIndex)*eye(Ktotal));
                    
                    %Select users based on the modified CBSUS algorithm (as
                    %described in Footnote 12)
                    [scheduledZF,H_ZF] = functionUserSelectionCBSUS(Hestimate,N,1,Ktotal,P(powerIndex),Eest);
                    
                    %Compute approximate ZF precoding for the selected users
                    [~,W_ZFC] = functionBlockDiagonalization(H_ZF,N,1,length(scheduledZF),P(powerIndex));
                    
                    %Placeholder for the effective channel after MMSE receive combining
                    HeffectiveMMSE = zeros(length(scheduledZF),N);
                    
                    %Go through the scheduled users
                    for k = 1:length(scheduledZF)
                        Hk = H((scheduledZF(k)-1)*M+1:scheduledZF(k)*M,:); %Channel of scheduled User k
                        
                        %Compute MMSE receive combiner as described in Remark 1
                        cMMSE = (eye(M)+Hk*(W_ZFC*W_ZFC')*Hk')\(Hk*W_ZFC(:,k));
                        
                        %Compute the effective channel with the MMSE receive combiner
                        if norm(cMMSE) > 0
                            cMMSE = cMMSE/norm(cMMSE);
                            HeffectiveMMSE(k,:) = cMMSE'*Hk;
                        end
                        
                    end
                    
                    %Compute sum rate for ZF with MMSE receive combining
                    channelGains = abs(HeffectiveMMSE*W_ZFC).^2;
                    signalPowers = diag(channelGains);
                    interferencePowers = sum(channelGains,2)-signalPowers;
                    sumRateZFC(n,powerIndex,l,m) = sum(log2(1+signalPowers./(1+interferencePowers)));
                    
                end
                
            end
            
        end
        
    end
    
end


%Prepare to plot the results
%Compute the average sum rates for all transmission schemes
averageSumRateBD = reshape(mean(sumRateBD,1),[length(P) length(Kvalues) length(correlationFactors)]);
averageSumRateZFC = reshape(mean(sumRateZFC,1),[length(P) length(Kvalues) length(correlationFactors)]);

if imperfectChannelEstimation == false || variableLargeScaleFading == 2
    averageSumRateMET = reshape(mean(sumRateMET,1),[length(P) length(Kvalues) length(correlationFactors)]);
end


%Prepare the markers that are used to differentiate the spatial correlation
%factors in the plots.
plotMarkers = ['d','s','o']; %One marker for each spatial correlation factor
distanceMarkers = 3; %Fraction of simulation points that are marked


%Plot all the main results in one figure.
figure(simulateFigure); hold on; box on;

for m = 1:length(correlationFactors)
    
    for powerIndex = 1:length(PdB)
        
        if imperfectChannelEstimation == false || variableLargeScaleFading == 2
            plot(Kvalues,averageSumRateMET(powerIndex,:,m),'k--','LineWidth',1);
        end
        
        plot(Kvalues,averageSumRateZFC(powerIndex,:,m),'r-','LineWidth',1);
        plot(Kvalues,averageSumRateBD(powerIndex,:,m),'b-.','LineWidth',1);
        
        if imperfectChannelEstimation == false || variableLargeScaleFading == 2
            plot(Kvalues(2:distanceMarkers:end-1),averageSumRateMET(powerIndex,2:distanceMarkers:end-1,m),['k' plotMarkers(m)],'LineWidth',1);
        end
        plot(Kvalues(2:distanceMarkers:end-1),averageSumRateZFC(powerIndex,2:distanceMarkers:end-1,m),['r' plotMarkers(m)],'LineWidth',1);
        plot(Kvalues(2:distanceMarkers:end-1),averageSumRateBD(powerIndex,2:distanceMarkers:end-1,m),['b' plotMarkers(m)],'LineWidth',1);
    end
    
end


if simulateFigure == 4
    axis([8 50 0 65]);
elseif simulateFigure == 5 || simulateFigure == 6
    axis([8 50 50 105]);
elseif simulateFigure == 7
    axis([8 50 0 65]);
elseif simulateFigure == 8
    axis([8 50 30 95]);
elseif simulateFigure == 9
    axis([8 50 50 105]);
end

if imperfectChannelEstimation == false || variableLargeScaleFading == 2
    legend('MET (Greedy stream allocation)','1 stream/selected user (ZFC)','4 streams/selected user (BD)','Location','SouthEast');
elseif imperfectChannelEstimation == true && variableLargeScaleFading == 2
    legend('Greedy allocation (MET-like)','1 stream/selected user (ZFC)','4 streams/selected user (BD)','Location','SouthEast');
else
    legend('1 stream/selected user (ZFC)','4 streams/selected user (BD)','Location','SouthEast');
end

xlabel('Total Number of Users');
ylabel('Average Sum Rate [bit/channel use]');




%Plot Figure 6 from the paper, which explains some of the results from
%Figure 5. For simplicity, it is divided into three different figures (one
%for each spatial correlation factor).
if simulateFigure == 5 || simulateFigure == 6
    
    kIndex = 4; %Selects the results for K=20
    
    labels = cell(1,4);
    labels{1} = '1 stream';
    labels{2} = '2 streams';
    labels{3} = '3 streams';
    labels{4} = '4 streams';
    
    names = cell(1,3);
    names{1} = ['\rho = ' num2str(correlationFactors(1))];
    names{2} = ['\rho = ' num2str(correlationFactors(2))];
    names{3} = ['\rho = ' num2str(correlationFactors(3))];
    
    %Extract the number of streams that an active user get, in the
    %whole cell, at the cell center and at the cell edge.
    streamDistribution_Average = zeros(4,length(correlationFactors));
    streamDistribution_CellCenter = zeros(4,length(correlationFactors));
    streamDistribution_CellEdge = zeros(4,length(correlationFactors));
    
    for m = 1:length(correlationFactors)
        
        streamDistribution_Average(:,m) = mean(streamAllocation_Average(:,:,1,kIndex,m),2)/sum(mean(streamAllocation_Average(:,:,1,kIndex,m),2));
        streamDistribution_CellCenter(:,m) = mean(streamAllocation_CellCenter(:,:,1,kIndex,m),2)/sum(mean(streamAllocation_CellCenter(:,:,1,kIndex,m),2));
        streamDistribution_CellEdge(:,m) = mean(streamAllocation_CellEdge(:,:,1,kIndex,m),2)/sum(mean(streamAllocation_CellEdge(:,:,1,kIndex,m),2));
        
    end
    
    %Plot the results as bars
    figure(6);
    figureHandle = bar(streamDistribution_Average','stacked');
    colors = fliplr(colormap(hot));
    colormap(colors(5:64,:));
    axis([0.5 3.5 0 1.2]);
    set(gca,'xticklabel', names);
    legend(figureHandle,labels,'Orientation','Horizontal');
    title('Whole Cell');
    
    figure(7);
    figureHandle = bar(streamDistribution_CellCenter','stacked');
    colors = fliplr(colormap(hot));
    colormap(colors(5:64,:));
    axis([0.5 3.5 0 1.2]);
    set(gca,'xticklabel', names);
    legend(figureHandle,labels,'Orientation','Horizontal');
    title('Cell Center');
    
    figure(8);
    figureHandle = bar(streamDistribution_Average','stacked');
    colors = fliplr(colormap(hot));
    colormap(colors(5:64,:));
    axis([0.5 3.5 0 1.2]);
    set(gca,'xticklabel', names);
    legend(figureHandle,labels,'Orientation','Horizontal');
    title('Cell Edge');
    
end

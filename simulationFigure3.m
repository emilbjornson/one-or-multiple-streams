%This Matlab script can be used to generate Figure 3 in the article:
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
M = 2; %Number of receive antennas per user

Ktotal = N; %Number of users in this simulation

%SNR values (we consider an infinitely large SNR so these values are only
%used as dummy values in the precoding optimization).
PdB = 50; %dB scale
P = 10.^(PdB/10); %Linear scale

%Number of realizations in the Monte Carlo simulations
nbrOfMonteCarloRealizations = 5000;

%Channel realizations for Monte Carlo simulations (before channel
%correlation has been applied).
Huncorrelated = (randn(M*Ktotal,N,nbrOfMonteCarloRealizations)+1i*randn(M*Ktotal,N,nbrOfMonteCarloRealizations))/sqrt(2);
thetasRx = 2*pi*rand(nbrOfMonteCarloRealizations,Ktotal); %Angle to the base station as seen from each user
thetasTx = 2*pi*rand(nbrOfMonteCarloRealizations,Ktotal); %Angle to each user as seen from the base station

%Range of spatial correlation factors
correlationFactors = [0:0.01:1-1e-3 1-1e-3];


%Placeholders for storing the simulation results
rateOffsetBD_Rx = zeros(nbrOfMonteCarloRealizations,length(correlationFactors)); %Block-diagonalization, correlation only at receiver side
rateOffsetZF_Rx = zeros(nbrOfMonteCarloRealizations,length(correlationFactors)); %Zero-forcing, correlation only at receiver side

rateOffsetBD_Tx = zeros(nbrOfMonteCarloRealizations,length(correlationFactors)); %Block-diagonalization, correlation only at transmitter side
rateOffsetZF_Tx = zeros(nbrOfMonteCarloRealizations,length(correlationFactors)); %Zero-forcing, correlation only at transmitter side

rateOffsetBD_RxTx = zeros(nbrOfMonteCarloRealizations,length(correlationFactors)); %Block-diagonalization, correlation at both transmitter and receiver sides
rateOffsetZF_RxTx = zeros(nbrOfMonteCarloRealizations,length(correlationFactors)); %Zero-forcing, correlation at both transmitter and receiver sides


%Go through all Monte Carlo realizations
for n = 1:nbrOfMonteCarloRealizations
    
    %Write out the progress at every 10 iterations
    if mod(n,10)==0
        disp(['Progress: ' num2str(n) ' out of ' num2str(nbrOfMonteCarloRealizations) ' realizations.']);
    end
    
    %Go through all correlation factors
    for c = 1:length(correlationFactors)
        
        
        %Placeholders for the correlated channels
        HcorrRx = zeros(M*Ktotal,N); %Case 1: Only correlation at receiver
        HcorrTx = zeros(M*Ktotal,N); %Case 2: Only correlation at transmitter
        HcorrRxTx = zeros(M*Ktotal,N); %Case 3: Correlation at both transmitter and receivers
        
        %Compute the correlated channel realizations for each user
        for k = 1:Ktotal
            %Compute correlation matrices in accordance to the
            %exponential correlation model in Eq. (16) for the current
            %spatial correlation factor. Each user has its own theta.
            Rr = toeplitz( ( abs(correlationFactors(c)) * exp(1i*thetasRx(n,k)) ).^(0:M-1) );
            Rt = toeplitz( ( abs(correlationFactors(c)) * exp(1i*thetasTx(n,k)) ).^(0:N-1) );
            
            %Case 1: Correlation only at receiver side
            HcorrRx((k-1)*M+1:k*M,:) = sqrtm(Rr) * Huncorrelated((k-1)*M+1:k*M,:,n);
            
            %Case 2: Correlation only at transmitter side
            HcorrTx((k-1)*M+1:k*M,:) = Huncorrelated((k-1)*M+1:k*M,:,n) * sqrtm(Rt);
            
            %Case 3: Correlation at both transmitter and receiver sides
            HcorrRxTx((k-1)*M+1:k*M,:) = sqrtm(Rr) * Huncorrelated((k-1)*M+1:k*M,:,n) * sqrtm(Rt);
        end
        
        
        
        %Case 1: Compute high-SNR rate offsets
        [~,~,Wnorm_BD] = functionBlockDiagonalization(HcorrRx(1:N,:),N,M,N/M,P); %Compute optimal precoding matrices for block-diagonalization
        rateOffsetBD_Rx(n,c) = real(log2(det(HcorrRx(1:N,:)*(Wnorm_BD*Wnorm_BD')*HcorrRx(1:N,:)'))); %Compute and store the rate offset term in Eq. (10)
        
        %Compute effective channel when the users apply receive
        %combining along the strong singular direction.
        Heffective = zeros(N,N); %Placeholder for effective channel
        for k = 1:Ktotal
            [U,S,V] = svd(HcorrRx((k-1)*M+1:k*M,:),'econ'); %Compute singular value decomposition for User k
            Heffective(k,:) = U(:,1)'*HcorrRx((k-1)*M+1:k*M,:); %Compute and store the effective channel after receive combining
        end
        
        [~,~,Wnorm_ZF] = functionBlockDiagonalization(Heffective,N,1,N,P); %Compute optimal precoding matrices for zero-forcing with receive combining
        rateOffsetZF_Rx(n,c) = real(log2(det(Heffective*(Wnorm_ZF*Wnorm_ZF')*Heffective'))); %Compute and store the rate offset term in Eq. (10)
        
        
        
        %Case 2: Compute high-SNR rate offsets
        [~,~,Wnorm_BD] = functionBlockDiagonalization(HcorrTx(1:N,:),N,M,N/M,P); %Compute optimal precoding matrices for block-diagonalization
        rateOffsetBD_Tx(n,c) = real(log2(det(HcorrTx(1:N,:)*(Wnorm_BD*Wnorm_BD')*HcorrTx(1:N,:)'))); %Compute and store the rate offset term in Eq. (10)
        
        %Compute effective channel when the users apply receive
        %combining along the strong singular direction.
        Heffective = zeros(N,N); %Placeholder for effective channel
        for k = 1:Ktotal
            [U,S,V] = svd(HcorrTx((k-1)*M+1:k*M,:),'econ'); %Compute singular value decomposition for User k
            Heffective(k,:) = U(:,1)'*HcorrTx((k-1)*M+1:k*M,:); %Compute and store the effective channel after receive combining
        end
        
        [~,~,Wnorm_ZF] = functionBlockDiagonalization(Heffective,N,1,N,P); %Compute optimal precoding matrices for zero-forcing with receive combining
        rateOffsetZF_Tx(n,c) = real(log2(det(Heffective*(Wnorm_ZF*Wnorm_ZF')*Heffective'))); %Compute and store the rate offset term in Eq. (10)
        
        
        
        %Case 3: Compute high-SNR rate offsets
        [~,~,Wnorm_BD] = functionBlockDiagonalization(HcorrRxTx(1:N,:),N,M,N/M,P); %Compute optimal precoding matrices for block-diagonalization
        rateOffsetBD_RxTx(n,c) = real(log2(det(HcorrRxTx(1:N,:)*(Wnorm_BD*Wnorm_BD')*HcorrRxTx(1:N,:)'))); %Compute and store the rate offset term in Eq. (10)
        
        %Compute effective channel when the users apply receive
        %combining along the strong singular direction.
        Heffective = zeros(N,N); %Placeholder for effective channel
        for k = 1:Ktotal
            [U,S,V] = svd(HcorrRxTx((k-1)*M+1:k*M,:),'econ'); %Compute singular value decomposition for User k
            Heffective(k,:) = U(:,1)'*HcorrRxTx((k-1)*M+1:k*M,:); %Compute and store the effective channel after receive combining
        end
        
        [~,~,Wnorm_ZF] = functionBlockDiagonalization(Heffective,N,1,N,P); %Compute optimal precoding matrices for zero-forcing with receive combining
        rateOffsetZF_RxTx(n,c) = real(log2(det(Heffective*(Wnorm_ZF*Wnorm_ZF')*Heffective'))); %Compute and store the rate offset term in Eq. (10)
        
    end
end


%Placeholders for upper and lower bounds on the high-SNR difference between
%block-diagonalization and zero-forcing with receive combining
upperBounds = zeros(length(correlationFactors),1);
lowerBounds = zeros(length(correlationFactors),1);

%Go through all correlation factors
for c = 1:length(correlationFactors)
    
    %Compute the eigenvalues of the receive correlation matrix for the
    %current spatial correlation factor. Note that the eigenvalues are
    %independent of theta, which thus has been set to zero.
    Rrx = toeplitz( ( abs(correlationFactors(c)) ).^(0:M-1) );
    eigenvaluesRx = eig(Rrx);
    
    %Compute upper bound based on Eq. (11) in Theorem 1
    upperBounds(c) = N*(log2(exp(1))/M*sum((M-1:-1:1)./(1:M-1)) + log2(geomean(eigenvaluesRx)) - log2(max(eigenvaluesRx)));
    
    %Compute expected value of the squared spectral norm of the channel,
    %which is the effective strength of the channel under zero-forcing.
    effectiveChannelNorm = zeros(nbrOfMonteCarloRealizations,Ktotal);
    for n = 1:nbrOfMonteCarloRealizations
        for k = 1:Ktotal
            effectiveChannelNorm(n,k) = norm(sqrtm(Rrx)*Huncorrelated(1+M*(k-1):M*k,:,n),2).^2; %Compute squared spectral norm
        end
    end
    expectedValue = mean(effectiveChannelNorm(:));
    
    %Compute upper bound based on Eq. (11) in Theorem 1
    lowerBounds(c) = N*(log2(exp(1))/M*sum((M-1:-1:1)./(1:M-1)) + log2(geomean(eigenvaluesRx)) + psi(N)/log(2) - log2(expectedValue) );
    
end

%Compute the expected asymptotic differences
asymptoticDifference_Rx = mean(rateOffsetBD_Rx-rateOffsetZF_Rx);
asymptoticDifference_Tx = mean(rateOffsetBD_Tx-rateOffsetZF_Tx);
asymptoticDifference_RxTx = mean(rateOffsetBD_RxTx-rateOffsetZF_RxTx);


%Plot the expected asymptotic differences
figure; hold on; box on;

plot(correlationFactors,upperBounds,'r--','LineWidth',1);
plot(correlationFactors,asymptoticDifference_Rx,'r','LineWidth',1);
plot(correlationFactors,asymptoticDifference_Tx,'b-.','LineWidth',1);
plot(correlationFactors,asymptoticDifference_RxTx,'k:','LineWidth',1);
plot(correlationFactors,lowerBounds,'r--','LineWidth',1);

axis([0 1 -30 10]);

legend('Rx-Corr (Bounds)','Rx-Corr','Tx-Corr','Both','Location','SouthWest');

xlabel('Spatial Correlation Factor');
ylabel('Expected Asymptotic Difference (BD-ZF)')

%This Matlab script can be used to generate Figure 11 in the article:
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

N = 6; %Number of transmit antennas at base station
M = 2; %Number of receive antennas per user

Ktotal = N; %Number of users in this simulation

%Number of quantization bit per channel dimension (5 bit/user with ZFC and
%10 bit/user with BD).
B = 5; 

%Range of SNR values
PdB = 5:5:30; %dB scale
P = 10.^(PdB/10); %Linear scale

%Number of realizations in the Monte Carlo simulations
nbrOfMonteCarloRealizations = 1000;

%Channel realizations for Monte Carlo simulations
Huncorrelated = (randn(M*Ktotal,N,nbrOfMonteCarloRealizations)+1i*randn(M*Ktotal,N,nbrOfMonteCarloRealizations))/sqrt(2);
 

%Placeholders for storing the simulation results
sumRateSU = zeros(nbrOfMonteCarloRealizations,length(P));
sumRateBD = zeros(nbrOfMonteCarloRealizations,length(P));
sumRateZF_QBC_MMSE = zeros(nbrOfMonteCarloRealizations,length(P));
sumRateZF_QBC = zeros(nbrOfMonteCarloRealizations,length(P));


%Go through all Monte Carlo realizations
for n = 1:nbrOfMonteCarloRealizations
    
    %Write out the progress at every 10 iterations
    if mod(n,10)==0
        disp(['Progress: ' num2str(n) ' out of ' num2str(nbrOfMonteCarloRealizations) ' realizations.']);
    end
    
    
    %Current channel realization
    H = Huncorrelated(:,:,n); %True channel
    
    %Quantized channel under block diagonalization or single-user
    %transmission, based on random vector quantization (RVQ)
    HquantizedBD = functionRVQfeedback_BD(H(1:Ktotal,:),2^(B*M),N,M,Ktotal/M); 
    
    %Quantized channel under zero-forcing with combining, based on 
    %random vector quantization (RVQ)
    [HquantizedZFC,Heffective]=functionRVQfeedback_ZF_QBC(H,2^B,N,M,Ktotal);
    
    
    %Scheme: Single-user transmission (1 user, M antennas on user)
    [~,W_SU] = functionBlockDiagonalization(HquantizedBD(1:M,:),N,M,1,P); %Compute precoding based on quantized channel
    sumRateSU(n,:) = functionSumrateComputation(H(1:M,:),W_SU,M,1); %Compute sum rate achieved 
    
    
    %Scheme: Block-diagonalization (N/M users, M antennas/user)
    [~,W_BD] = functionBlockDiagonalization(HquantizedBD,N,M,N/M,P); %Compute precoding based on quantized channel
    sumRateBD(n,:) = functionSumrateComputation(H(1:N,:),W_BD,M,N/M)'; %Compute sum rate that is achieved 
    
    
    %Scheme: Zero-forcing with combining (N users, 1 antennas/user)
    [~,W_ZFC] = functionBlockDiagonalization(HquantizedZFC,N,1,N,P); %Compute precoding based on quantized channel
    sumRateZF_QBC(n,:) = functionSumrateComputation(Heffective,W_ZFC,1,N)'; %Compute sum rate that is achieved 
    
    %Go through all transmit powers
    for powerIndex = 1:length(P)
        
        %Placeholder for the effective channel after MMSE receive combining
        HeffectiveMMSE = zeros(N,N);
        
        %Go through all users
        for k = 1:N
            
            Hk = H((k-1)*M+1:k*M,:); %Channel to User k
            
            %Compute MMSE receive combiner as described in Remark 1
            cMMSE = (eye(M)+Hk*(W_ZFC(:,:,powerIndex)*W_ZFC(:,:,powerIndex)')*Hk')\(Hk*W_ZFC(:,k,powerIndex)); 
            
            %Compute the effective channel with MMSE receive combiner
            if norm(cMMSE) > 0
                cMMSE = cMMSE/norm(cMMSE);
                HeffectiveMMSE(k,:) = cMMSE'*Hk;
            end
            
        end
        
        %Compute sum rate for ZF with MMSE receive combining
        channelGains = abs(HeffectiveMMSE*W_ZFC(:,:,powerIndex)).^2;
        signalPowers = diag(channelGains);
        interferencePowers = sum(channelGains,2)-signalPowers;
        sumRateZF_QBC_MMSE(n,powerIndex) = sum(log2(1+signalPowers./(1+interferencePowers)));
        
    end
    
end


%Compute the average sum rates for all transmission schemes
averageSumRateSU = mean(sumRateSU,1);
averageSumRateBD = mean(sumRateBD,1);
average_sumrateZF_QBC_MMSE = mean(sumRateZF_QBC_MMSE,1);
average_sumrateZF_QBC = mean(sumRateZF_QBC,1);


%Plot the results
figure; hold on; box on;

plot(PdB,averageSumRateSU,'k:','LineWidth',1);
plot(PdB,averageSumRateBD,'b-.','LineWidth',1);
plot(PdB,average_sumrateZF_QBC_MMSE,'r','LineWidth',1);
plot(PdB,average_sumrateZF_QBC,'k--','LineWidth',1);

legend('Single-User Transmission','BD','ZFC (QBC+MMSE)','ZFC (QBC)','Location','Best');
xlabel('Total Transmit Power (dB)');
ylabel('Average Sum Rate [bit/channel use]');

axis([5 30 0 22]);

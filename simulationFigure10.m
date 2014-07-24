%This Matlab script can be used to generate Figure 10 in the article:
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

N = 4; %Number of transmit antennas at base station
M = 2; %Number of receive antennas per user

Ktotal = 8; %Number of users in this simulation

%Range of numbers of quantization bit per channel dimension
quantizationBits = 2:11;


%Compute the range of SNR values from the number of quantization bits,
%based on the formula in [7, Eq. (17)] that maintain a 3 dB gap between BD
%with perfect and quantized CSI. Note that we use another notation thatn in
%[7], that is, M and N have the opposite meanings.

T = M*(N-M); %Parameter defined below Eq. (8) in [7]
C_NM = prod(factorial(N-(1:M))./factorial(M-(1:M)))/factorial(T); %Another parameter defined below Eq. (8) in [7]
C_NMprime = M^(M*(N-M))*C_NM;%Parameter defined below Eq. (17) in [7]
PdB = ( M*quantizationBits-log2(C_NMprime) ) / (M*(N-M)/(10*log10(2))); %Compute power in dB scale according to Eq. (17) in [7]. Note that M*quantizationBits is the total number of bits per user in our notation
P = 10.^(PdB/10); %Linear scale


%Number of realizations in the Monte Carlo simulations
nbrOfMonteCarloRealizations = 100;

%Channel realizations for Monte Carlo simulations
Huncorrelated = (randn(M*Ktotal,N,nbrOfMonteCarloRealizations)+1i*randn(M*Ktotal,N,nbrOfMonteCarloRealizations))/sqrt(2);


%Placeholders for storing the simulation results
sumRateSU = zeros(nbrOfMonteCarloRealizations,length(P));
sumRateBD = zeros(nbrOfMonteCarloRealizations,length(P));
sumRateBDperfectCSI = zeros(nbrOfMonteCarloRealizations,length(P));
sumRateZFC = zeros(nbrOfMonteCarloRealizations,length(P));


%Go through all Monte Carlo realizations
for n = 1:nbrOfMonteCarloRealizations
    
    %Write out the progress at every iterations
    disp(['Progress: ' num2str(n) ' out of ' num2str(nbrOfMonteCarloRealizations) ' realizations.']);
    
    
    %Current channel realization
    H = Huncorrelated(:,:,n); %True channel
    
    
    %Scheme: Block-diagonalization with perfect CSI
    for powerIndex = 1:length(P)
        
        [scheduledBD,H_BD] = functionUserSelectionCBSUS(H,N,M,Ktotal,P(powerIndex)); %Compute the sum rate maximizing user selection
        sumRateBDperfectCSI(n,powerIndex) = functionBlockDiagonalization(H_BD,N,M,length(scheduledBD),P(powerIndex)); %Compute the corresponding sum rate
        
    end
    
    
    for powerIndex = 1:length(P)
        
        %The number of quantization bits increase with the SNR according to
        %[7, Eq. (17)], to maintain an approximate 3 dB from BD with
        %perfect CSI.
        B = quantizationBits(powerIndex);
        
        %Quantized channel under block diagonalization or single-user
        %transmission, based on random vector quantization (RVQ)
        HquantizedBD = functionRVQfeedback_BD(H(1:Ktotal,:),2^(B*M),N,M,Ktotal/M);
        
        
        %Scheme: Single-user transmission with imperfect CSI to a randomly selected/first user (1 user, M antennas on user)
        [~,W_SU] = functionBlockDiagonalization(HquantizedBD(1:M,:),N,M,1,P(powerIndex));
        sumRateSU(n,powerIndex) = functionSumrateComputation(H(1:M,:),W_SU,M,1);
        
        
        %Scheme: Block diagonalization with imperfect CSI and user selection
        
        D_BD = (gamma(1/T)/T)*(2^(-B*M)*C_NM)^(-1/T); %Compute average quantization distortion as in Eq. (21).
        E_quantBD = (N/T)*D_BD * eye(Ktotal); %Compute approximative quantization error covariance matrix as described in Footnote 13
        
        [scheduledBD,Hquant_BD] = functionUserSelectionCBSUS(HquantizedBD,N,M,Ktotal/M,P(powerIndex),E_quantBD); %Select users based on the CBSUS algorithm using quantized channels
        [~,W_BD] = functionBlockDiagonalization(Hquant_BD,N,M,length(scheduledBD),P(powerIndex)); %Compute precoding for the selected users based on quantized channels
        
        %Pick out the true channels for the selected users
        Hselected = zeros(length(scheduledBD)*M,N);
        for m = 1:length(scheduledBD)
            Hselected((m-1)*M+1:m*M,:) = H((scheduledBD(m)-1)*M+1:scheduledBD(m)*M,:);
        end
        
        sumRateBD(n,powerIndex) = functionSumrateComputation(Hselected,W_BD,M,length(scheduledBD)); %Compute sum rate that is achieved
        
        
        
        %Quantized channel under zero-forcing with combining, based on
        %random vector quantization (RVQ)
        HquantizedZFC = functionRVQfeedback_ZF_MESC(H,2^B,N,M,Ktotal,P(powerIndex));
        
        
        %Scheme: Zero-forcing with imperfect CSI and user selection
        
        D_QBC = 2^(-B/(N-M)) * nchoosek(N-1,M-1)^(-1/(N-M)); %Compute average quantization distortion as in Eq. (24).
        
        E_quantZFC = (N-M+1)/(N-1) * D_QBC * eye(Ktotal); %Compute approximative quantization error covariance matrix as described in Footnote 13. Note that G = (N-M+1) for uncorrelated channels; see Lemma 3 in [8].
        
        [scheduledZFC,HquantizedZFC] = functionUserSelectionCBSUS(HquantizedZFC,N,1,Ktotal,P(powerIndex),E_quantZFC); %Select users based on the CBSUS algorithm using quantized channels
        
        [~,W_ZFC] = functionBlockDiagonalization(HquantizedZFC,N,1,length(scheduledZFC),P(powerIndex)); %Compute precoding for the selected users based on quantized channels
        
        %Placeholder for the effective channel after MMSE receive combining
        HeffectiveMMSE = zeros(length(scheduledZFC),N);
        
        for k = 1:length(scheduledZFC)
            
            Hk = H((scheduledZFC(k)-1)*M+1:scheduledZFC(k)*M,:); %Channel to selected user k
            
            %Compute MMSE receive combiner as described in Remark 1
            cMMSE = (eye(M)+Hk*(W_ZFC*W_ZFC')*Hk')\(Hk*W_ZFC(:,k));
            
            %Compute the effective channel with MMSE receive combiner
            if norm(cMMSE)>0
                cMMSE=cMMSE/norm(cMMSE);
                HeffectiveMMSE(k,:)=cMMSE'*Hk;
            end
            
        end
        
        
        %Compute sum rate for ZF with MMSE receive combining
        channelGains = abs(HeffectiveMMSE*W_ZFC).^2;
        signalPowers = diag(channelGains);
        interferencePowers = sum(channelGains,2)-signalPowers;
        sumRateZFC(n,powerIndex) = sum(log2(1+signalPowers./(1+interferencePowers)));
        
    end
    
end


%Compute the average sum rates for all transmission schemes
averageSumRateSU = mean(sumRateSU,1);
averageSumRateBD = mean(sumRateBD,1);
averageSumRateBDperfectCSI = mean(sumRateBDperfectCSI,1);
averageSumRateZFC = mean(sumRateZFC,1);


%Plot the results
figure; hold on; box on;

plot(PdB,averageSumRateBDperfectCSI,'b--','LineWidth',1);
plot(PdB,averageSumRateZFC,'r','LineWidth',1);
plot(PdB,averageSumRateBD,'b-.','LineWidth',1);
plot(PdB,averageSumRateSU,'k:','LineWidth',1);

legend('BD, Perfect CSI','ZFC','BD','Single-User Transmission','Location','Best');
xlabel('Total Transmit Power (dB)');
ylabel('Average Sum Rate [bit/channel use]');

axis([0 15 0 18]);

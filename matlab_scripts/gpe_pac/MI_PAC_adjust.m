%% function for calculating the PAC modulation index
function MI = MI_PAC_adjust(raw_phase,raw_amp)

num_phabins = size(raw_phase,1);
num_ampbins = size(raw_amp{1},1);
frepair = cell(num_phabins,num_ampbins);
for i = 1:num_phabins
    for j = 1:num_ampbins
        frepair{i,j} = [raw_phase(i,:);raw_amp{i}(j,:)];
        
    end
end

avg_frepair = cell(num_phabins,num_ampbins);
nBins = 18;
binEdges = linspace(-pi,pi,nBins+1);
%binCenters = binEdges(1:end-1)-diff(binEdges)/2;

for i = 1:num_phabins
    for j = 1:num_ampbins
        betaphase = frepair{i,j}(1,:);
        [~,binIdx]=histc(betaphase,binEdges);
        amp = frepair{i,j}(2,:);
        avg_amp = zeros(1,nBins);
        avg_pha = zeros(1,nBins);
        for bin = 1:nBins
            if any(binIdx(:)==bin)
                avg_amp(bin)=mean(amp(binIdx(:)==bin));
                avg_pha(bin)=mean(betaphase(binIdx(:)==bin));
            end
        end
        avg_frepair{i,j}(1,:)= avg_pha;
        avg_frepair{i,j}(2,:)=avg_amp;
    end
end

MI = zeros(num_phabins,num_ampbins);
for i = 1:num_phabins
    for j = 1:num_ampbins
        phaseBin = avg_frepair{i,j}(1,:);
        ampBin = avg_frepair{i,j}(2,:);
        ampSum = sum(ampBin);
        ampNorm = ampBin./ampSum;
        ampQ = ones(1,nBins)./nBins;
        if any(ampNorm==0)
            ampNorm(ampNorm==0)=eps;
        end
        distKL =sum(ampNorm.*log(ampNorm./ampQ));
        MI(i,j) = distKL./log(nBins);
    end
end

end
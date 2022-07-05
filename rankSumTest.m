function PseudoCAP=rankSumTest(PseudoCAP)
%Wilcoxon rank sum test for Individual functional categories
%'PseudoCAP'is a data struct contains gene information(espression level and noise) of individual functional categories
for i=1:length(PseudoCAP)
    if PseudoCAP(i).EXP_medianGFP>PseudoCAP(i).EXP_otherGFP_median
        PseudoCAP(i).P_EXP_GFP=-log10(ranksum(PseudoCAP(i).EXP_meanGFP,PseudoCAP(i).EXP_otherGFP,'tail','right'));
    else
        PseudoCAP(i).P_EXP_GFP=log10(ranksum(PseudoCAP(i).EXP_meanGFP,PseudoCAP(i).EXP_otherGFP,'tail','left'));
    end
    if PseudoCAP(i).STA_medianGFP>PseudoCAP(i).STA_otherGFP_median
        PseudoCAP(i).P_STA_GFP=-log10(ranksum(PseudoCAP(i).STA_meanGFP,PseudoCAP(i).STA_otherGFP,'tail','right'));
    else
        PseudoCAP(i).P_STA_GFP=log10(ranksum(PseudoCAP(i).STA_meanGFP,PseudoCAP(i).STA_otherGFP,'tail','left'));
    end
    %%
    if PseudoCAP(i).EXP_medianNoise>PseudoCAP(i).EXP_otherNoise_median
        PseudoCAP(i).P_EXP_Noise=-log10(ranksum(PseudoCAP(i).EXP_NoiseTotalGFP,PseudoCAP(i).EXP_otherNoise,'tail','right'));
    else
        PseudoCAP(i).P_EXP_Noise=log10(ranksum(PseudoCAP(i).EXP_NoiseTotalGFP,PseudoCAP(i).EXP_otherNoise,'tail','left'));
    end
    if PseudoCAP(i).STA_medianNoise>PseudoCAP(i).STA_otherNoise_median
        PseudoCAP(i).P_STA_Noise=-log10(ranksum(PseudoCAP(i).STA_NoiseTotalGFP,PseudoCAP(i).STA_otherNoise,'tail','right'));
    else
        PseudoCAP(i).P_STA_Noise=log10(ranksum(PseudoCAP(i).STA_NoiseTotalGFP,PseudoCAP(i).STA_otherNoise,'tail','left'));
    end
      %%
    if PseudoCAP(i).EXP_medianIntNoise>PseudoCAP(i).EXP_otherIntNoise_median
        PseudoCAP(i).P_EXP_IntNoise=-log10(ranksum(PseudoCAP(i).EXP_NoiseIntGFP,PseudoCAP(i).EXP_otherIntNoise,'tail','right'));
    else
        PseudoCAP(i).P_EXP_IntNoise=log10(ranksum(PseudoCAP(i).EXP_NoiseIntGFP,PseudoCAP(i).EXP_otherIntNoise,'tail','left'));
    end
    if PseudoCAP(i).STA_medianIntNoise>PseudoCAP(i).STA_otherIntNoise_median
        PseudoCAP(i).P_STA_IntNoise=-log10(ranksum(PseudoCAP(i).STA_NoiseIntGFP,PseudoCAP(i).STA_otherIntNoise,'tail','right'));
    else
        PseudoCAP(i).P_STA_IntNoise=log10(ranksum(PseudoCAP(i).STA_NoiseIntGFP,PseudoCAP(i).STA_otherIntNoise,'tail','left'));
    end
    %%
%     if PseudoCAP(i).EXP_medianResiduals>PseudoCAP(i).EXP_otherResiduals_median
%         PseudoCAP(i).P_EXP_Residuals=-log10(ranksum(PseudoCAP(i).EXP_TotalResiduals,PseudoCAP(i).EXP_otherResiduals,'tail','right'));
%     else
%         PseudoCAP(i).P_EXP_Residuals=log10(ranksum(PseudoCAP(i).EXP_TotalResiduals,PseudoCAP(i).EXP_otherResiduals,'tail','left'));
%     end
%     if PseudoCAP(i).STA_medianResiduals>PseudoCAP(i).STA_otherResiduals_median
%         PseudoCAP(i).P_STA_Residuals=-log10(ranksum(PseudoCAP(i).STA_TotalResiduals,PseudoCAP(i).STA_otherResiduals,'tail','right'));
%     else
%         PseudoCAP(i).P_STA_Residuals=log10(ranksum(PseudoCAP(i).STA_TotalResiduals,PseudoCAP(i).STA_otherResiduals,'tail','left'));
%     end
    %%
    if PseudoCAP(i).EXP_medianIntResiduals>PseudoCAP(i).EXP_otherIntResiduals_median
        PseudoCAP(i).P_EXP_IntResiduals=-log10(ranksum(PseudoCAP(i).EXP_IntResiduals,PseudoCAP(i).EXP_otherIntResiduals,'tail','right'));
    else
        PseudoCAP(i).P_EXP_IntResiduals=log10(ranksum(PseudoCAP(i).EXP_IntResiduals,PseudoCAP(i).EXP_otherIntResiduals,'tail','left'));
    end
    if PseudoCAP(i).STA_medianIntResiduals>PseudoCAP(i).STA_otherIntResiduals_median
        PseudoCAP(i).P_STA_IntResiduals=-log10(ranksum(PseudoCAP(i).STA_IntResiduals,PseudoCAP(i).STA_otherIntResiduals,'tail','right'));
    else
        PseudoCAP(i).P_STA_IntResiduals=log10(ranksum(PseudoCAP(i).STA_IntResiduals,PseudoCAP(i).STA_otherIntResiduals,'tail','left'));
    end
end
end


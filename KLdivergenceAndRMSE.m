% Calculating KL divergence and RMSE values of  each strains
parfor iData=1:length(AllData)
    xi=linspace(min(AllData(iData).EXP_sfGFP),max(AllData(iData).EXP_sfGFP),102);
    xi=xi(2:101);
    [f,xi]=ksdensity(AllData(iData).EXP_sfGFP,xi,'Support','positive','Function','cdf');
    pd1=fitdist(AllData(iData).EXP_sfGFP,'Gamma');
    cdf_Gamma=cdf(pd1,xi);
    pd2=fitdist(AllData(iData).EXP_sfGFP,'Lognormal');
    cdf_Lognormal=cdf(pd2,xi);
    AllData(iData).EXP_RMSE_Gamma=calculateRMSE(f,cdf_Gamma);
    AllData(iData).EXP_RMSE_Lognormal=calculateRMSE(f,cdf_Lognormal);
    AllData(iData).EXP_KLdivergence_Gamma=calculateKLdivergence(f,cdf_Gamma);
    AllData(iData).EXP_KLdivergence_Lognormal=calculateKLdivergence(f,cdf_Lognormal);    
    %%
    xi=linspace(min(AllData(iData).STA_sfGFP),max(AllData(iData).STA_sfGFP),102);
    xi=xi(2:101);
    [f,xi]=ksdensity(AllData(iData).STA_sfGFP,xi,'Support','positive','Function','cdf');
    pd1=fitdist(AllData(iData).STA_sfGFP,'Gamma');
    cdf_Gamma=cdf(pd1,xi);
    pd2=fitdist(AllData(iData).STA_sfGFP,'Lognormal');
    cdf_Lognormal=cdf(pd2,xi);
    AllData(iData).STA_RMSE_Gamma=calculateRMSE(f,cdf_Gamma);
    AllData(iData).STA_RMSE_Lognormal=calculateRMSE(f,cdf_Lognormal);
    AllData(iData).STA_KLdivergence_Gamma=calculateKLdivergence(f,cdf_Gamma);
    AllData(iData).STA_KLdivergence_Lognormal=calculateKLdivergence(f,cdf_Lognormal);
    disp(iData)
end

function RMSEvalue=calculateRMSE(x1,x2)
RMSEvalue=sqrt(mean((x1-x2).^2));
end

function KLvalue=calculateKLdivergence(x1,x2)
x1=x1/sum(x1);
x2=x2/sum(x2);
KLvalue=sum(x1.*(log2(x1)-log2(x2)));
end


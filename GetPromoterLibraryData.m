function GetPromoterLibraryData(strainNum,montageNum)
% This function reads the information from the mask image and obtains data such as expressions and noise.
%Input parameter 'strainNum' should be two numbers, the first is the number of strains taken in a round, the second is the number of rounds taken, 
%input parameter 'montageNum' should be the number of fields of view containing the stationary
% and logarithmic phases [montageNumSTA,montageNumEXP],GeneInfoFile should be an excel file containing the id information of the gene.
%Example; 24 strains were taken at once, then 2 rounds were taken with a stable phase field of view of 16 and a logarithmic phase field of view of 25, GetPromoterLibraryData([24,2],[16,25]);
clc
disp('Please select the folder where the data is located')
pause(1)
dirSaveFile=uigetdir();
disp(['Results storage folder:',dirSaveFile])
dirFileList=dir(dirSaveFile);
for iFile=1:numel(dirFileList)-2
    if ~isempty(strfind(dirFileList(iFile+2).name,'EXP'))||~isempty(strfind(dirFileList(iFile+2).name,'exp'))...
            ||~isempty(strfind(dirFileList(iFile+2).name,'DS'))||~isempty(strfind(dirFileList(iFile+2).name,'ds'))&&dirFileList(iFile+2).isdir==1
        dirFileEXP=[dirSaveFile,'\',dirFileList(iFile+2).name];%Folders for log phase data
        disp(['Folder where log period data is located：',32,dirFileEXP])
    elseif ~isempty(strfind(dirFileList(iFile+2).name,'STA'))||~isempty(strfind(dirFileList(iFile+2).name,'sta'))...
            ||~isempty(strfind(dirFileList(iFile+2).name,'WD'))||~isempty(strfind(dirFileList(iFile+2).name,'wd'))&&dirFileList(iFile+2).isdir==1
        dirFileSTA=[dirSaveFile,'\',dirFileList(iFile+2).name];%Folder for stationary phase data
        disp(['Folder where stable period data is located：',32,dirFileSTA])
    end
end
if exist([dirSaveFile,'\GeneInfo.xlsx'],'file')==0
    disp('No files with gene name information found')
    return
end
GeneInfoFile=[dirSaveFile,'\GeneInfo.xlsx'];%excel spreadsheet of genetic information
strainDataFile=[dirSaveFile,'\montageData'];
mkdir([dirSaveFile,'\montageData'])
[~,GeneInfo,~]=xlsread(GeneInfoFile);
if size(GeneInfo,1)~=strainNum(1)*strainNum(2)
    disp('there is something wrong in the GeneInfoFile')
    return
end
[~]=scanBasicExperimentResultcwh(dirFileEXP,'exp',strainNum,montageNum(2),strainDataFile);
[~]=scanBasicExperimentResultcwh(dirFileSTA,'sta',strainNum,montageNum(1),strainDataFile);
plotPromoterLibraryData(strainDataFile,strainNum(1)*strainNum(2),GeneInfo)
end

function dirAllFile=scanBasicExperimentResultcwh(dirAllFile,phase,strainNum,montageNum,strainDataFile)
%Obtaining information on the fluorescence intensity of single bacteria
allNameList=dir(dirAllFile);
load([dirAllFile,'\mip.mat']);
 fluoImageInitializedcwh(dirAllFile,mip,phase);
DataHere=[];
for iField=1:numel(allNameList)-2
    if ~(strcmp(allNameList(iField+2).name(1:5),'feild') || strcmp(allNameList(iField+2).name(1:5),'felid')) && allNameList(iField+2).isdir==0%strcmp(allNameList(iField+2).name(end-2:end),'mat')
        continue
    end
    fieldNum=str2num(allNameList(iField+2).name(end-3:end));
    if ~isempty(fieldNum)
        disp(allNameList(iField+2).name)
        DataHere=[DataHere,fieldNum];
        dirFile=[dirAllFile,'\',allNameList(iField+2).name];
        dirImage=[dirFile,'\Tracking'];
        nameList=dir(dirImage);
        n=0;
        for i=1:numel(nameList)-3
            if strcmp(nameList(i+3).name(1:5),'image')
                n=n+1;
                temp=load([dirImage,'\',nameList(i+3).name]);
                maskImage=temp.imageTracking;
                maskImage=bwmorph(maskImage,'open');
                maskImage=bwmorph(maskImage,'spur');
                maskImage=bwmorph(maskImage,'majority');
                maskImage=bwmorph(maskImage,'hbreak');
                cc=regionprops(maskImage,'MajorAxisLength','MinorAxisLength','PixelIdxList','Area');
                if strcmp(phase,'exp')==1
                    for iCC=1:size(cc,1)
                        if cc(iCC).MajorAxisLength/cc(iCC).MinorAxisLength<=2 ||cc(iCC).MajorAxisLength<=25 ||cc(iCC).MajorAxisLength>=100 || cc(iCC).MinorAxisLength<= 9 ...
                                || cc(iCC).MinorAxisLength>=19 ||cc(iCC).Area>=1100
                            image1=maskImage;
                            image1(cc(iCC).PixelIdxList)=0;
                            maskImage=image1;
                        end
                    end
                else
                    for iCC=1:size(cc,1)
                        if cc(iCC).MajorAxisLength<15 ||cc(iCC).MajorAxisLength>50 || cc(iCC).MinorAxisLength< 12.5 ...
                                || cc(iCC).MinorAxisLength>20 ||cc(iCC).Area>=650
                            image1=maskImage;
                            image1(cc(iCC).PixelIdxList)=0;
                            maskImage=image1;
                        end
                    end
                end
                maskImage=bwareaopen(maskImage,200);
                maskImages(:,:,2*n-1)=maskImage;
                maskImages(:,:,2*n)=maskImage;
            end
        end
        bacNum=getBasicFigure(maskImages,dirFile,mip,fieldNum);
        clear maskImages
    end
end
clc
%%
if strainNum(2)==1
    if montageNum*strainNum(1)~=numel(DataHere)%The total number of fields should be multiplied by the number of strains multiplied by the number of fields of a single strain
        disp('there is something wrong about the field')
    else
        DataHere=sort(DataHere);
        getmontageData1(dirAllFile,montageNum,DataHere,strainDataFile,phase);%Collect data for the first image only
    end
else
    if montageNum*strainNum(1)~=numel(DataHere)
        disp('there is something wrong about the field')
    else
        DataHere=sort(DataHere);
        getmontageData2(dirAllFile,montageNum,DataHere,strainDataFile,phase);
    end
end
if exist([dirAllFile,'\imageStack'],'dir')==0
    image2stack(strainNum,montageNum,dirAllFile);%Save images as stack
end
end
function image2stack(strainNum,montageNum,dirfile)
%Save images as stacks
% dirfile=uigetdir();
mkdir([dirfile,'\imageStack'])
if strainNum(2)==1
    for iStrain=1:strainNum(1)
        mkdir([dirfile,'\imageStack\Strain',num2str(iStrain,'%03.f')])
        for i=(iStrain-1)*montageNum+1:iStrain*montageNum
            imageCyOFP=import_tiff_stack([dirfile,'\feild',num2str(i,'%04.f'),'\CyOFP\imageCyOFP00001.tif']);
            imagesfGFP=import_tiff_stack([dirfile,'\feild',num2str(i,'%04.f'),'\sfGFP\imagesfGFP00001.tif']);
            imwrite(imageCyOFP,[dirfile,'\imageStack\Strain',num2str(iStrain,'%03.f'),'\CyOFP.tif'],'WriteMode','append')
            imwrite(imagesfGFP,[dirfile,'\imageStack\Strain',num2str(iStrain,'%03.f'),'\sfGFP.tif'],'WriteMode','append')
            clear imageCyOFP
            clear imagesfGFP
        end
    end
else
    for iStrain=1:strainNum(1)
        mkdir([dirfile,'\imageStack\Strain',num2str(iStrain,'%03.f')])
        for i=(iStrain-1)*montageNum+1:iStrain*montageNum
            imageCyOFP=import_tiff_stack([dirfile,'\feild',num2str(i,'%04.f'),'\CyOFP\imageCyOFP00001.tif']);
            imagesfGFP=import_tiff_stack([dirfile,'\feild',num2str(i,'%04.f'),'\sfGFP\imagesfGFP00001.tif']);
            imwrite(imageCyOFP,[dirfile,'\imageStack\Strain',num2str(iStrain,'%03.f'),'\CyOFP.tif'],'WriteMode','append')
            imwrite(imagesfGFP,[dirfile,'\imageStack\Strain',num2str(iStrain,'%03.f'),'\sfGFP.tif'],'WriteMode','append')
            clear imageCyOFP
            clear imagesfGFP
        end
    end
    for iStrain=strainNum(1)+1:2*strainNum(1)
        mkdir([dirfile,'\imageStack\Strain',num2str(iStrain,'%03.f')])
        for i=(iStrain-25)*montageNum+1:(iStrain-24)*montageNum
            imageCyOFP=import_tiff_stack([dirfile,'\feild',num2str(i,'%04.f'),'\CyOFP\imageCyOFP00002.tif']);
            imagesfGFP=import_tiff_stack([dirfile,'\feild',num2str(i,'%04.f'),'\sfGFP\imagesfGFP00002.tif']);
            imwrite(imageCyOFP,[dirfile,'\imageStack\Strain',num2str(iStrain,'%03.f'),'\CyOFP.tif'],'WriteMode','append')
            imwrite(imagesfGFP,[dirfile,'\imageStack\Strain',num2str(iStrain,'%03.f'),'\sfGFP.tif'],'WriteMode','append')
            clear imageCyOFP
            clear imagesfGFP
        end
    end
end
end

function bacNum=getBasicFigure(maskImages,dirFile,mip,fieldNum)
%Acquisition of fluorescence intensity of SfGFP images and CyOFP images
fluoChannel{1,1}='CyOFP';
fluoChannel{2,1}='sfGFP';
fluoChannel{3,1}='mScarletI';
fluoChannel{4,1}='TDsmURFP';
fluoChannel{5,1}='CyPet';
fluoChannel{6,1}='Venus';
fluoChannel{7,1}='mAmetrine';
if isempty(maskImages)
    for iChannel=1:numel(fluoChannel)
        result=[];
        try
            load([dirFile,'\',fluoChannel{iChannel},'new\','frameInfo.mat'])
        catch err
            continue
        end
        switch iChannel
            case 1
                bacInfo{1}.meanCyOFP=1;
            case 2
                bacInfo{1}.meanGFP=1;
            case 3
                bacInfo{1}.meanmScalet=1;
            case 4
                bacInfo{1}.meanRFP=1;
            case 5
                bacInfo{1}.meanCyPet=1;
            case 6
                bacInfo{1}.meanVenus=1;
            case 7
                bacInfo{1}.meanmAmetrine=1;
        end
    end
    bacInfo{1}.majorLength=1;
    bacInfo{1}.minorLength=1;
    bacInfo{1}.time=1;
    if ~isempty(bacInfo)
        bacInfo{1}.tag=mip.feildTag.tag{fieldNum};
        bacInfo{1}.tagValue=mip.feildTag.tagValue(fieldNum);
    end
    save([dirFile,'\bacInfo&tree'],'bacInfo');
    bacNum=[];
    return
end
load([dirFile,'\CyOFP\frameInfo.mat']);
timeBegin=frameInfo(1,1:6);
for iTime=1:size(frameInfo,1)
    timeAll(iTime)=etime(frameInfo(iTime,1:6),timeBegin)/60;
end
maskImages=maskImages(:,:,2:2:end);
% plot bacNum vs time
for iImage=1:size(maskImages,3)
    image=maskImages(:,:,iImage);
    cc=bwconncomp(image);
    bacNum(iImage)=cc.NumObjects;
    cc=regionprops(image,'MajorAxisLength','MinorAxisLength','Area');
    majorLength=[];
    minorLength=[];
    Area=[];
    for iSmallBac=1:numel(cc)
        majorLength=[majorLength;cc(iSmallBac).MajorAxisLength];
        minorLength=[minorLength;cc(iSmallBac).MinorAxisLength];
        Area=[Area;cc(iSmallBac).Area];
    end
    bacInfo{iImage}.majorLength=majorLength;
    bacInfo{iImage}.minorLength=minorLength;
    bacInfo{iImage}.Area=Area;
    bacInfo{iImage}.time=timeAll(iImage);
end
for iChannel=1:numel(fluoChannel)
    bioTreeFrame=[];
    result=[];
    try
        load([dirFile,'\',fluoChannel{iChannel},'new\','frameInfo.mat'])
    catch err
        continue
    end
    for i=1:size(frameInfo,1)
        bioTreeTimer=etime(frameInfo(i,1:6),timeBegin)/60;
        diffTime=abs(timeAll-bioTreeTimer);
        [~,index]=find(diffTime==min(diffTime));
        bioTreeFrame(i)=index(1);
    end
    for iPic=1:numel(bioTreeFrame)
        if bioTreeFrame(iPic)==0
            continue
        end
        image=load([dirFile,'\',fluoChannel{iChannel},'new','\image',fluoChannel{iChannel},num2str((iPic),'%05.f'),'.mat']);
        image=image.imageI;
        maskI=maskImages(:,:,bioTreeFrame(iPic));
        cc=regionprops(maskI,image,'MeanIntensity');
        meanIntensity=[];
        for iSmallBac=1:numel(cc)
            meanIntensity=[meanIntensity;cc(iSmallBac).MeanIntensity];
        end
        switch iChannel
            case 1
                bacInfo{bioTreeFrame(iPic)}.meanCyOFP=meanIntensity;
            case 2
                bacInfo{bioTreeFrame(iPic)}.meanGFP=meanIntensity;
            case 3
                bacInfo{bioTreeFrame(iPic)}.meanmScalet=meanIntensity;
            case 4
                bacInfo{bioTreeFrame(iPic)}.meanRFP=meanIntensity;
            case 5
                bacInfo{bioTreeFrame(iPic)}.meanCyPet=meanIntensity;
            case 6
                bacInfo{bioTreeFrame(iPic)}.meanVenus=meanIntensity;
            case 7
                bacInfo{bioTreeFrame(iPic)}.meanmAmetrine=meanIntensity;
        end
        result=[result;timeAll(iPic)/60,mean(image(maskI))];
    end
    if ~isempty(bacInfo)
        bacInfo{1}.tag=mip.feildTag.tag{fieldNum};
        bacInfo{1}.tagValue=mip.feildTag.tagValue(fieldNum);
    end
end
fluoChannel{8,1}='RedMask';
fluoChannel{9,1}='BlueMask';
fluoChannel{10,1}='GreenMask';
for iChannel=8:10
    bioTreeFrame=[];
    result=[];
    try
        load([dirFile,'\',fluoChannel{iChannel},'\','frameInfo.mat'])
    catch err
        continue
    end
    for i=1:size(frameInfo,1)-1
        deltaTime(i)=abs(etime(frameInfo(i+1,1:6),frameInfo(i,1:6)));
    end
    deltaTime(end+1)=deltaTime(end);
    frameInfo(:,15)=frameInfo(:,15).*frameInfo(:,9)./deltaTime'/1000;
    for i=1:size(frameInfo,1)
        bioTreeTimer=etime(frameInfo(i,1:6),timeBegin)/60;
        diffTime=abs(timeAll-bioTreeTimer);
        [~,index]=find(diffTime==min(diffTime));
        bioTreeFrame(i)=index(1);
    end
    for iPic=1:numel(bioTreeFrame)
        if bioTreeFrame(iPic)==0
            continue
        end
        image=import_tiff_stack([dirFile,'\',fluoChannel{iChannel}(1:end-4),'Control\image',fluoChannel{iChannel}(1:end-4),'Control',num2str(1,'%05.f'),'.tif']);
        maskI=maskImages(:,:,bioTreeFrame(iPic));
        cc=regionprops(maskI,image,'MeanIntensity');
        meanIntensity=[];
        for iSmallBac=1:numel(cc)
            meanIntensity=[meanIntensity;cc(iSmallBac).MeanIntensity];
        end
        switch iChannel
            case 8
                bacInfo{bioTreeFrame(iPic)}.redControl=meanIntensity*16/frameInfo(iPic,14)*frameInfo(iPic,15);
            case 9
                bacInfo{bioTreeFrame(iPic)}.blueControl=meanIntensity*16/frameInfo(iPic,14)*frameInfo(iPic,15);
            case 10
                bacInfo{bioTreeFrame(iPic)}.greenControl=meanIntensity*16/frameInfo(iPic,14)*frameInfo(iPic,15);
        end
    end
end
save([dirFile,'\bacInfo&tree'],'bacInfo');
end


function trackingFluoChannelReGet(dirFieldFile,fluoHere)
%Converting fluorescence intensity information into protein concentration
for i=1:numel(fluoHere)
    fileNum(i)=numel(dir([dirFieldFile,'\',fluoHere{i}]));
end
if max(fileNum)==min(fileNum)
    return
end
trackingNum=find(fileNum==max(fileNum));
fluoNum=find(fileNum==min(fileNum));
fluoFrameInfo=load([dirFieldFile,'\',fluoHere{fluoNum(1)},'\frameInfo.mat']);
fluoFrameInfo=fluoFrameInfo.frameInfo;
trackingFrameInfo=load([dirFieldFile,'\',fluoHere{trackingNum(1)},'\frameInfo.mat']);
trackingFrameInfo=trackingFrameInfo.frameInfo;
movefile([dirFieldFile,'\',fluoHere{trackingNum(1)}],[dirFieldFile,'\',fluoHere{trackingNum(1)},'Pre']);
mkdir([dirFieldFile,'\',fluoHere{trackingNum(1)}]);
for i=1:size(fluoFrameInfo,1)
    timeI=fluoFrameInfo(i,1:6);
    trackingTime=trackingFrameInfo(:,1:6);
    for j=1:size(trackingTime,1)
        deltaTime(j)=etime(trackingTime(j,:),timeI);
    end
    focusTime=find(abs(deltaTime)==min(abs(deltaTime)));
    copyfile([dirFieldFile,'\',fluoHere{trackingNum(1)},'Pre\image',fluoHere{trackingNum(1)},num2str(focusTime,'%05.f'),'.tif'],[dirFieldFile,'\',fluoHere{trackingNum(1)},'\image',fluoHere{trackingNum(1)},num2str(i,'%05.f'),'.tif']);
    frameInfo(i,:)=trackingFrameInfo(focusTime,:);
end
save([dirFieldFile,'\',fluoHere{trackingNum(1)},'\frameInfo.mat'],'frameInfo');
end

function fluoImageInitializedcwh(dirFile,mip,phase)
fluoChannel=mip.Calib.fluo_protein_info;
fieldNameList=dir(dirFile);
sita=144.9;
load([dirFile,'\LumencorField.mat']);
fluoHere=[fluoChannel(4),fluoChannel(5)];
backGround=LumencorField.CyanField;
for iField=1:numel(fieldNameList)-2
    %     fieldNum=0;
    dirFieldFile=[dirFile,'\',fieldNameList(iField+2).name];
    if numel(fieldNameList(iField+2).name)>=9 && strcmp(fieldNameList(iField+2).name(end-8:end-4),'feild')
        %         fieldNum=fieldNum+1;
        mkdir([dirFieldFile,'\',fluoHere{1},'new']);
        copyfile([dirFieldFile,'\',fluoHere{1},'\frameInfo.mat'],[dirFieldFile,'\',fluoHere{1},'new\frameInfo.mat'])
        mkdir([dirFieldFile,'\',fluoHere{2},'new']);
        copyfile([dirFieldFile,'\',fluoHere{2},'\frameInfo.mat'],[dirFieldFile,'\',fluoHere{2},'new\frameInfo.mat'])
        trackingFluoChannelReGet(dirFieldFile,fluoHere);
        for i=1:numel(dir([dirFieldFile,'\',fluoHere{1}]))-3
            for iFluo=1:numel(fluoHere)
                load([dirFieldFile,'\',fluoHere{iFluo},'\frameInfo.mat']);
                image=import_tiff_stack([dirFieldFile,'\',fluoHere{iFluo},'\image',fluoHere{iFluo},num2str(i,'%05.f'),'.tif']);
                image=double(image)-110;
                image=image./backGround;
                image=image/frameInfo(i,14)*frameInfo(i,13)/frameInfo(i,9);
                image=image(:);
                newMatrix(iFluo,:)=double(image);
                clear image
            end
            switch phase
                case 'sta'
                    florescenceCorrectionMatrix=[1,0.0254;0.0743,1];%String light correction matrix from the slope of the y=a*x fit after normalisation of the fluorescence of a single fluorescent strain
                case 'exp'
                    florescenceCorrectionMatrix=[1,0.0765;0.0904,1];%String light correction matrix from the slope a of the y=a*x fit after fluorescence normalization of single fluorescent strains, 20190329
            end
            newMatrix=florescenceCorrectionMatrix\newMatrix;
            caliMatrix=[1.89234564955976e-05,0;0,3.70463569321944e-06];%Changed the string light factor to 0
            %             newMatrix=inv(caliMatrix)*newMatrix;
            newMatrix=caliMatrix\newMatrix;
            for iFluo=1:numel(fluoHere)
                imageI=reshape(newMatrix(iFluo,:),2048,2048)*sita/1000;
                save([dirFieldFile,'\',fluoHere{iFluo},'new\image',fluoHere{iFluo},num2str(i,'%05.f'),'.mat'],'imageI')
            end
        end
    end
end
end

function imageStack= import_tiff_stack( fname )
%Reading images
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
warning off all
infoImage=imfinfo(fname);
frameNum=size(infoImage,1);
imageWith=infoImage(1).Width;
imageHeight=infoImage(1).Height;
imageBit=infoImage(1).BitsPerSample;
imageBit=strcat('uint',num2str(imageBit));
imageStack=zeros(imageHeight,imageWith,imageBit);
imageCurrent=Tiff(fname,'r');
for iframe=1:frameNum
    imageCurrent.setDirectory(iframe);
    imageStack(:,:,iframe)= imageCurrent.read();
end
imageCurrent.close();
end

function getmontageData1(dirFile,montageNum,DataHere,saveFile,phase)
%Calculate the expression and noise of each strain
for iStrain=1:(numel(DataHere)/montageNum)
    if strcmp(phase,'sta')==1
        subSavefile=[saveFile,'\Strain',num2str(iStrain,'%03.f'),'\stationary phase'];
    else
        subSavefile=[saveFile,'\Strain',num2str(iStrain,'%03.f'),'\exponential phase'];
    end
    mkdir(subSavefile)
    montageData.bacNum=0;
    montageData.majorLength=[];
    montageData.minorLength=[];
    montageData.Area=[];
    montageData.CyOFP=[];
    montageData.sfGFP=[];
    montageData.FLRatio=[];
    for iField=2:montageNum
        if iField==2
            copyfile([dirFile,'\feild',num2str(DataHere((iStrain-1)*montageNum+iField),'%04.f'),'\sfGFP\imagesfGFP00001.tif'],[subSavefile,'\imagesfGFP.tif'])
        end
        load([dirFile,'\feild',num2str(DataHere((iStrain-1)*montageNum+iField),'%04.f'),'\bacInfo&tree.mat']);
        disp(DataHere((iStrain-1)*montageNum+iField))
        montageData.bacNum=montageData.bacNum+length(bacInfo{1}.meanCyOFP);
        montageData.majorLength=[montageData.majorLength;bacInfo{1}.majorLength];
        montageData.minorLength=[montageData.minorLength;bacInfo{1}.minorLength];
        montageData.Area=[montageData.Area;bacInfo{1}.Area];
        montageData.CyOFP=[montageData.CyOFP;bacInfo{1}.meanCyOFP];
        montageData.sfGFP=[montageData.sfGFP;bacInfo{1}.meanGFP];
        montageData.FLRatio=[montageData.FLRatio;bacInfo{1}.meanGFP./bacInfo{1}.meanCyOFP];
        clear bacInfo
    end
    clc
    montageData.meanCyOFP=mean(abs(montageData.CyOFP));
    montageData.meansfGFP=mean(abs(montageData.sfGFP));
    montageData.meanFLRatio=mean(abs(montageData.FLRatio));
    %Noise calculation
    %External noise
    NoiseExt=mean(montageData.CyOFP.*montageData.sfGFP)/(montageData.meanCyOFP*montageData.meansfGFP)-1;%NoiseExt=(<GR>/(<G><R>))-1
    montageData.NoiseExt=NoiseExt;
    %CyOFP
    stdCyOFP=std(montageData.CyOFP);
    NoiseTotalCyOFP=(stdCyOFP/montageData.meanCyOFP)^2;
    NoiseIntCyOFP=NoiseTotalCyOFP-NoiseExt;%NoiseTotal^2=NoiseExt^2+NoiseInt^2
    montageData.NoiseTotalCyOFP=NoiseTotalCyOFP;
    montageData.NoiseIntCyOFP=NoiseIntCyOFP;
    %sfGFP
    stdsfGFP=std(montageData.sfGFP);
    NoiseTotalsfGFP=(stdsfGFP/montageData.meansfGFP)^2;
    NoiseIntsfGFP=NoiseTotalsfGFP-NoiseExt;
    montageData.NoiseTotalsfGFP=NoiseTotalsfGFP;
    montageData.NoiseIntsfGFP=NoiseIntsfGFP;
    save([subSavefile,'\result'],'montageData')
    
    montageData=GetPositiveData(montageData);
    if montageData.isPositive
        montageData.meanCyOFP=mean(montageData.CyOFP);
        montageData.meansfGFP=mean(montageData.sfGFP);
        montageData.meanFLRatio=mean(montageData.FLRatio);
        %External noise
        NoiseExt=mean(montageData.CyOFP.*montageData.sfGFP)/(montageData.meanCyOFP*montageData.meansfGFP)-1;%NoiseExt=(<GR>/(<G><R>))-1,
        montageData.NoiseExt=NoiseExt;
        %CyOFP
        stdCyOFP=std(montageData.CyOFP);
        NoiseTotalCyOFP=(stdCyOFP/montageData.meanCyOFP)^2;
        NoiseIntCyOFP=NoiseTotalCyOFP-NoiseExt;%NoiseTotal^2=NoiseExt^2+NoiseInt^2
        montageData.NoiseTotalCyOFP=NoiseTotalCyOFP;
        montageData.NoiseIntCyOFP=NoiseIntCyOFP;
        %sfGFP
        stdsfGFP=std(montageData.sfGFP);
        NoiseTotalsfGFP=(stdsfGFP/montageData.meansfGFP)^2;
        NoiseIntsfGFP=NoiseTotalsfGFP-NoiseExt;
        montageData.NoiseTotalsfGFP=NoiseTotalsfGFP;
        montageData.NoiseIntsfGFP=NoiseIntsfGFP;
    end
    save([subSavefile,'\resultNew'],'montageData')
    clear montageData
end
end

function getmontageData2(dirFile,montageNum,DataHere,saveFile,phase)
%Calculate the expression and noise of each strain
subStrainNum=numel(DataHere)/montageNum;
for iStrain=1:(numel(DataHere)/montageNum)
    if strcmp(phase,'sta')==1
        subSavefile1=[saveFile,'\Strain',num2str(iStrain,'%03.f'),'\stationary phase'];
        subSavefile2=[saveFile,'\Strain',num2str(iStrain+subStrainNum,'%03.f'),'\stationary phase'];
    else
        subSavefile1=[saveFile,'\Strain',num2str(iStrain,'%03.f'),'\exponential phase'];
        subSavefile2=[saveFile,'\Strain',num2str(iStrain+subStrainNum,'%03.f'),'\exponential phase'];
    end
    mkdir(subSavefile1)
    mkdir(subSavefile2)
    %1
    montageData1.bacNum=0;
    montageData1.majorLength=[];
    montageData1.minorLength=[];
    montageData1.Area=[];
    montageData1.CyOFP=[];
    montageData1.sfGFP=[];
    montageData1.FLRatio=[];
    %2
    montageData2.bacNum=0;
    montageData2.majorLength=[];
    montageData2.minorLength=[];
    montageData2.Area=[];
    montageData2.CyOFP=[];
    montageData2.sfGFP=[];
    montageData2.FLRatio=[];
    for iField=2:montageNum
        if iField==2
            copyfile([dirFile,'\feild',num2str(DataHere((iStrain-1)*montageNum+iField),'%04.f'),'\sfGFP\imagesfGFP00001.tif'],[subSavefile1,'\imagesfGFP.tif'])
            copyfile([dirFile,'\feild',num2str(DataHere((iStrain-1)*montageNum+iField),'%04.f'),'\sfGFP\imagesfGFP00002.tif'],[subSavefile2,'\imagesfGFP.tif'])
        end
        load([dirFile,'\feild',num2str(DataHere((iStrain-1)*montageNum+iField),'%04.f'),'\bacInfo&tree.mat']);
        disp(DataHere((iStrain-1)*montageNum+iField))
        %1
        montageData1.bacNum=montageData1.bacNum+length(bacInfo{1}.meanCyOFP);
        montageData1.majorLength=[montageData1.majorLength;bacInfo{1}.majorLength];
        montageData1.minorLength=[montageData1.minorLength;bacInfo{1}.minorLength];
        montageData1.Area=[montageData1.Area;bacInfo{1}.Area];
        montageData1.CyOFP=[montageData1.CyOFP;bacInfo{1}.meanCyOFP];
        montageData1.sfGFP=[montageData1.sfGFP;bacInfo{1}.meanGFP];
        montageData1.FLRatio=[montageData1.FLRatio;bacInfo{1}.meanGFP./bacInfo{1}.meanCyOFP];
        %2
        montageData2.bacNum=montageData2.bacNum+length(bacInfo{2}.meanCyOFP);
        montageData2.majorLength=[montageData2.majorLength;bacInfo{2}.majorLength];
        montageData2.minorLength=[montageData2.minorLength;bacInfo{2}.minorLength];
        montageData2.Area=[montageData2.Area;bacInfo{2}.Area];
        montageData2.CyOFP=[montageData2.CyOFP;bacInfo{2}.meanCyOFP];
        montageData2.sfGFP=[montageData2.sfGFP;bacInfo{2}.meanGFP];
        montageData2.FLRatio=[montageData2.FLRatio;bacInfo{2}.meanGFP./bacInfo{2}.meanCyOFP];
        clear bacInfo
    end
    clc
    %1
    montageData1.meanCyOFP=mean(montageData1.CyOFP);
    montageData1.meansfGFP=mean(montageData1.sfGFP);
    montageData1.meanFLRatio=mean(montageData1.FLRatio);
    %2
    montageData2.meanCyOFP=mean(montageData2.CyOFP);
    montageData2.meansfGFP=mean(montageData2.sfGFP);
    montageData2.meanFLRatio=mean(montageData2.FLRatio);
    %Noise calculation
    %External noise
    %1
    NoiseExt1=mean(montageData1.CyOFP.*montageData1.sfGFP)/(montageData1.meanCyOFP*montageData1.meansfGFP)-1;
    montageData1.NoiseExt=NoiseExt1;
    %2
    NoiseExt2=mean(montageData2.CyOFP.*montageData2.sfGFP)/(montageData2.meanCyOFP*montageData2.meansfGFP)-1;
    montageData2.NoiseExt=NoiseExt2;
    %
    %CyOFP
    %1
    stdCyOFP1=std(montageData1.CyOFP);
    NoiseTotalCyOFP1=(stdCyOFP1/montageData1.meanCyOFP)^2;%CyOFP total noise N=CV=std/mean
    NoiseIntCyOFP1=NoiseTotalCyOFP1-NoiseExt1;%NoiseTotal^2=NoiseExt^2+NoiseInt^2
    montageData1.NoiseTotalCyOFP=NoiseTotalCyOFP1;
    montageData1.NoiseIntCyOFP=NoiseIntCyOFP1;
    %2
    stdCyOFP2=std(montageData2.CyOFP);
    NoiseTotalCyOFP2=(stdCyOFP2/montageData2.meanCyOFP)^2;%CyOFP total noise N=CV=std/mean
    NoiseIntCyOFP2=NoiseTotalCyOFP2-NoiseExt2;%NoiseTotal=NoiseExt+NoiseInt
    montageData2.NoiseTotalCyOFP=NoiseTotalCyOFP2;
    montageData2.NoiseIntCyOFP=NoiseIntCyOFP2;
    %sfGFP
    %1
    stdsfGFP1=std(montageData1.sfGFP);
    NoiseTotalsfGFP1= (stdsfGFP1/montageData1.meansfGFP)^2;
    NoiseIntsfGFP1=NoiseTotalsfGFP1-NoiseExt1;
    montageData1.NoiseTotalsfGFP=NoiseTotalsfGFP1;
    montageData1.NoiseIntsfGFP=NoiseIntsfGFP1;
    montageData=montageData1;
    save([subSavefile1,'\result'],'montageData')
    clear montageData
    %2
    stdsfGFP2=std(montageData2.sfGFP);
    NoiseTotalsfGFP2=(stdsfGFP2/montageData2.meansfGFP)^2;
    NoiseIntsfGFP2=NoiseTotalsfGFP2-NoiseExt2;
    montageData2.NoiseTotalsfGFP=NoiseTotalsfGFP2;
    montageData2.NoiseIntsfGFP=NoiseIntsfGFP2;
    montageData=montageData2;
    save([subSavefile2,'\result'],'montageData')
    clear montageData
    %1
    montageData1=GetPositiveData(montageData1);
    if montageData1.isPositive
        montageData1.meanCyOFP=mean(montageData1.CyOFP);
        montageData1.meansfGFP=mean(montageData1.sfGFP);
        montageData1.meanFLRatio=mean(montageData1.FLRatio);
        NoiseExt1=mean(montageData1.CyOFP.*montageData1.sfGFP)/(montageData1.meanCyOFP*montageData1.meansfGFP)-1;%NoiseExt=(<GR>/(<G><R>))-1,CyOFP and sfGFP external noise equal
        montageData1.NoiseExt=NoiseExt1;
        %CyOFP
        stdCyOFP1=std(montageData1.CyOFP);
        NoiseTotalCyOFP1=(stdCyOFP1/montageData1.meanCyOFP)^2;%CyOFP total noise N=CV=std/mean
        NoiseIntCyOFP1=NoiseTotalCyOFP1-NoiseExt1;%NoiseTotal^2=NoiseExt^2+NoiseInt^2
        montageData1.NoiseTotalCyOFP=NoiseTotalCyOFP1;
        montageData1.NoiseIntCyOFP=NoiseIntCyOFP1;
        %sfGFP
        stdsfGFP1=std(montageData1.sfGFP);
        NoiseTotalsfGFP1=(stdsfGFP1/montageData1.meansfGFP)^2;
        NoiseIntsfGFP1=NoiseTotalsfGFP1-NoiseExt1;
        montageData1.NoiseTotalsfGFP=NoiseTotalsfGFP1;
        montageData1.NoiseIntsfGFP=NoiseIntsfGFP1;
    end
    montageData=montageData1;
    save([subSavefile1,'\resultNew'],'montageData')
    clear montageData
    %2
    montageData2=GetPositiveData(montageData2);
    if montageData2.isPositive
        montageData2.meanCyOFP=mean(montageData2.CyOFP);
        montageData2.meansfGFP=mean(montageData2.sfGFP);
        montageData2.meanFLRatio=mean(montageData2.FLRatio);
        NoiseExt2=mean(montageData2.CyOFP.*montageData2.sfGFP)/(montageData2.meanCyOFP*montageData2.meansfGFP)-1;%NoiseExt=(<GR>/(<G><R>))-1,CyOFP and sfGFP external noise equal
        montageData2.NoiseExt=NoiseExt2;
        %CyOFP
        stdCyOFP2=std(montageData2.CyOFP);
        NoiseTotalCyOFP2=(stdCyOFP2/montageData2.meanCyOFP)^2;%CyOFP total noise N=CV=std/mean
        NoiseIntCyOFP2=NoiseTotalCyOFP2-NoiseExt2;%NoiseTotal=NoiseExt+NoiseInt
        montageData2.NoiseTotalCyOFP=NoiseTotalCyOFP2;
        montageData2.NoiseIntCyOFP=NoiseIntCyOFP2;
        %sfGFP
        stdsfGFP2=std(montageData2.sfGFP);
        NoiseTotalsfGFP2=(stdsfGFP2/montageData2.meansfGFP)^2;
        NoiseIntsfGFP2=NoiseTotalsfGFP2-NoiseExt2;
        montageData2.NoiseTotalsfGFP=NoiseTotalsfGFP2;
        montageData2.NoiseIntsfGFP=NoiseIntsfGFP2;
    end
    montageData=montageData2;
    save([subSavefile2,'\resultNew'],'montageData')
    clear montageData
end
end

function montageData=GetPositiveData(montageData)
montageDataNew.bacNum=0;
montageDataNew.majorLength=[];
montageDataNew.minorLength=[];
montageDataNew.Area=[];
montageDataNew.CyOFP=[];
montageDataNew.sfGFP=[];
montageDataNew.FLRatio=[];
montageDataNew.isPositive=0;
for iBac=1:montageData.bacNum
    if montageData.sfGFP(iBac)>=0
        montageDataNew.bacNum=montageDataNew.bacNum+1;
        montageDataNew.majorLength=[montageDataNew.majorLength;montageData.majorLength(iBac)];
        montageDataNew.minorLength=[montageDataNew.minorLength;montageData.minorLength(iBac)];
        montageDataNew.Area=[montageDataNew.Area;montageData.Area(iBac)];
        montageDataNew.CyOFP=[montageDataNew.CyOFP;montageData.CyOFP(iBac)*montageData.Area(iBac)*montageData.minorLength(iBac)*0.065^3*1e-15*6.022e23*1e-6*2/3];
        montageDataNew.sfGFP=[montageDataNew.sfGFP;montageData.sfGFP(iBac)*montageData.Area(iBac)*montageData.minorLength(iBac)*0.065^3*1e-15*6.022e23*1e-6*2/3];
        montageDataNew.FLRatio=[montageDataNew.FLRatio;montageData.FLRatio(iBac)];
    end
end
if montageDataNew.bacNum/montageData.bacNum>=0.8
    montageDataNew.isPositive=1;
end
montageData=montageDataNew;
end



function out = scatplot_ys(x,y,method,radius)
% x,y are taken to be logarithmic, focusing more on the magnitude of change in logarithmic coordinates,
% and later the graph is transformed back to the original xy, x=10.^x;y=10.^y;

% x=log10(x);
% y=log10(y);
% Usage out = scatplot_modifiedys(x,y,'re');
% The default parameter method does not need to be entered and is changed to 'circles' and radius does not need to be entered and is defined by default
N=100;%Size of the grid
n=5; %Matrix for two-dimensional convolution filtering
po=3; %plot type
ms=10;% markersize
%New method 'rectangle' for density, with radius based on 5% of the current data point as the detection range

%The main program is mainly based on the function scatplot to get
%%
% Scatter plot with color indicating data density
%
% USAGE:
%   out = scatplot(x,y,method,radius,N,n,po,ms)
%   out = scatplot(x,y,dd)
%
% DESCRIPTION:
%   Draws a scatter plot with a colorscale
%   representing the data density computed
%   using three methods
%
% INPUT VARIABLES:
%   x,y - are the data points
%   method - is the method used to calculate data densities:
%       'circles' - uses circles with a determined area
%               centered at each data point
%       'squares' - uses squares with a determined area
%               centered at each data point
%       'voronoi' - uses voronoi cells to determin data densities
%               default method is 'voronoi'
%   radius - is the radius used for the circles or squares
%       used to calculate the data densities if
%       (Note: only used in methods 'circles' and 'squares'
%           default radius is sqrt((range(x)/30)^2 + (range(y)/30)^2)
%   N - is the size of the square mesh (N x N) used to
%       filter and calculate contours
%       default is 100
%   n - is the number of coeficients used in the 2-D
%       running mean filter
%       default is 5
%       (Note: if n is length(2), n(2) is tjhe number of
%       of times the filter is applied)
%   po - plot options:
%       0 - No plot
%       1 - plots only colored data points (filtered)
%       2 - plots colored data points and contours (filtered)
%       3 - plots only colored data points (unfiltered)
%       4 - plots colored data points and contours (unfiltered)
%           default is 1
%   ms - uses this marker size for filled circles
%       default is 4
%
% OUTPUT VARIABLE:
%   out - structure array that contains the following fields:
%       dd - unfiltered data densities at (x,y)
%       ddf - filtered data densities at (x,y)
%       radius - area used in 'circles' and 'squares'
%               methods to calculate densities
%       xi - x coordenates for zi matrix
%       yi - y coordenates for zi matrix
%       zi - unfiltered data densities at (xi,yi)
%       zif - filtered data densities at (xi,yi)
%       [c,h] = contour matrix C as described in
%           CONTOURC and a handle H to a contourgroup object
%       hs = scatter points handles
%
%Copy-Left, Alejandro Sanchez-Barba, 2005
%%
if nargin==0
    out=scatplotdemo;
    return
end
% if nargin<3 | isempty(method)
%     method = 'vo';
% end
%Change the default radius method to circular 'circles'
if nargin<3 || isempty(method)
    method = 'ci';
end
if isnumeric(method)
    gsp(x,y,method,2)
    return
else
    method = method(1:2);
end
if nargin<4 || isempty(radius)
    radius = sqrt((range(x)/400)^2 + (range(y)/400)^2);
    %     radius = sqrt((range(x)/50)^2 + (range(y)/50)^2);
end
%%

%%
%Correct data if necessary
x = x(:);
y = y(:);
%Asuming x and y match
idat = isfinite(x);
x = x(idat);
y = y(idat);
% figure;
holdstate = ishold;
if holdstate==0
    cla
end
hold on
%--------- Caclulate data density ---------
dd = datadensity(x,y,method,radius);
%------------- Gridding -------------------
xi = repmat(linspace(min(x),max(x),N),N,1);
yi = repmat(linspace(min(y),max(y),N)',1,N);
% [xi,yi]=meshgrid(linspace(min(x),max(x),N),linspace(min(y),max(y),N));
zi = griddata(x,y,dd,xi,yi);
%----- Bidimensional running mean filter -----
zi(isnan(zi)) = 0;
coef = ones(n(1),1)/n(1);
zif = conv2(coef,coef,zi,'same');
if length(n)>1
    for k=1:n(2)
        zif = conv2(coef,coef,zif,'same');
    end
end
%-------- New Filtered data densities --------
ddf = griddata(xi,yi,zif,x,y);
%----------- Plotting --------------------
switch po
    case {1,2}
        if po==2
            [c,h] = contour(xi,yi,zif);
            out.c = c;
            out.h = h;
        end %if
        hs = gsp(x,y,ddf,ms);
        out.hs = hs;
        colorbar
    case {3,4}
        if po>3
            [c,h] = contour(xi,yi,zi);
            out.c = c;
        end %if
        hs = gsp(x,y,dd,ms);
        out.hs = hs;
        colorbar
end %switch
%------Relocate variables and place NaN's ----------
dd(idat) = dd;
dd(~idat) = NaN;
ddf(idat) = ddf;
ddf(~idat) = NaN;
%--------- Collect variables ----------------
out.dd = dd;
out.ddf = ddf;
out.radius = radius;
out.xi = xi;
out.yi = yi;
out.zi = zi;
out.zif = zif;
if ~holdstate
    hold off
end
return
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function out=scatplotdemo
po = 3;
method = 'circles';
N = [];
n = [];
ms = 5;
x = randn(2000,1);
y = randn(2000,1);
radius = sqrt((range(x)/50)^2 + (range(y)/50)^2);
out=scatplot_ys(x,y,method,radius);

return
end
%~~~~~~~~~~ Data Density ~~~~~~~~~~~~~~
function dd = datadensity(x,y,method,r)
%Computes the data density (points/area) of scattered points
%Striped Down version
%
% USAGE:
%   dd = datadensity(x,y,method,radius)
%
% INPUT:
%   (x,y) -  coordinates of points
%   method - either 'squares','circles', or 'voronoi'
%       default = 'voronoi'
%   radius - Equal to the circle radius or half the square width
Ld = length(x);
dd = zeros(Ld,1);
switch method %Calculate Data Density
    % ys add rectangle method
    case 're'  %---- Using rectangle ----
        r=0.1;
        for k=1:Ld
            dd(k) = sum( x>((1-r)*x(k)) & x<((1+r)*x(k)) & y>((1-r)*y(k)) & y<((1+r)*y(k)) );
            area=(r*2*x(k))*(r*2*y(k));
            dd(k)=dd(k)/area;
        end %for
    case 'sq'  %---- Using squares ----
        for k=1:Ld
            dd(k) = sum( x>(x(k)-r) & x<(x(k)+r) & y>(y(k)-r) & y<(y(k)+r) );
        end %for
        area = (2*r)^2;
        dd = dd/area;
    case 'ci'
        for k=1:Ld
            dd(k) = sum( sqrt((x-x(k)).^2 + (y-y(k)).^2) < r );
        end
        area = pi*r^2;
        dd = dd/area;
    case 'vo'  %----- Using voronoi cells ------
        [v,c] = voronoin([x,y]);
        for k=1:length(c)
            %If at least one of the indices is 1,
            %then it is an open region, its area
            %is infinity and the data density is 0
            if all(c{k}>1)
                a = polyarea(v(c{k},1),v(c{k},2));
                dd(k) = 1/a;
            end %if
        end %for
end %switch
return
end

function varargout = gsp(x,y,c,ms)
%Graphs scattered poits


x=x(~isnan(c));
y=y(~isnan(c));
c=c(~isnan(c));
%define the color
colorindex=2;
if colorindex==1

    color1=[0,0.45,0.74];
    color2=[0.95,0.98,1];
end
if colorindex==2

    color1=[0.85,0.33,0.1];
    color2=[1,0.92,0.89];
end
I=64;

mymap(:,1)=linspace(color2(1),color1(1),I);
mymap(:,2)=linspace(color2(2),color1(2),I);
mymap(:,3)=linspace(color2(3),color1(3),I);

map = colormap(mymap);
map = colormap(jet(256));
ind = fix((c-min(c))/(max(c)-min(c))*(size(map,1)-1))+1;
h = [];
h=scatter(x,y,ms,map(ind,:),'filled','MarkerFaceAlpha',0.5);
if nargout==1
    varargout{1} = h;
end
return
end
%%
function plotPromoterLibraryData(DataFile,DataNum,GeneInfo)
for iStrain=1:DataNum
    montageDataSTA=load([DataFile,'\Strain',num2str(iStrain,'%03.f'),'\stationary phase\result.mat']);
    montageDataEXP=load([DataFile,'\Strain',num2str(iStrain,'%03.f'),'\exponential phase\result.mat']);
    imageEXPGFP=import_tiff_stack([DataFile,'\Strain',num2str(iStrain,'%03.f'),'\exponential phase\imagesfGFP.tif']);
    imageSTAGFP=import_tiff_stack([DataFile,'\Strain',num2str(iStrain,'%03.f'),'\stationary phase\imagesfGFP.tif']);
    %%
    if numel(montageDataEXP.montageData.CyOFP)<800||numel(montageDataSTA.montageData.CyOFP)<800
        disp(['Strain',num2str(iStrain,'%03.f'),32,'is so little data'])
        continue
    end
    %图1
    if size(GeneInfo,2)>=3&&~isempty(GeneInfo{iStrain,3})
        StrainInfo=[GeneInfo{iStrain,1},32,32,GeneInfo{iStrain,2},32,32,GeneInfo{iStrain,3}];
    else
        StrainInfo=[GeneInfo{iStrain,1},32,32,GeneInfo{iStrain,2}];
    end
    disp(['Strain',num2str(iStrain,'%03.f')])
    figure1=figure('OuterPosition',[-7 33 1936 1056],'Visible','off');
    %sfGFP表达量
    meanSFGFP=[montageDataEXP.montageData.meansfGFP,montageDataSTA.montageData.meansfGFP];
    Ymax=abs(max(meanSFGFP));
    x=categorical({'EXP PHASE','STA PHASE'});
    subplot1=subplot(2,2,1);
    hold(subplot1,'on')
    bar(subplot1,x,meanSFGFP,0.2)
    ylim(subplot1,[0,Ymax*1.1])
    % axis tight
    ylabel('sfGFP(uM)')
    set(gca,'FontSize',17);
    %图片
    imageEXPGFP=imageEXPGFP(512:1535,512:1535);
    imageSTAGFP=imageSTAGFP(512:1535,512:1535);
    imageR=zeros(1024,1024);
    imageB=zeros(1024,1024);
    RGBEXP=cat(3,imageR,imageEXPGFP,imageB);
    RGBSTA=cat(3,imageR,imageSTAGFP,imageB);
    subplot2=subplot(2,2,2);
    hold(subplot2,'on')
    imshowpair(imadjust(RGBEXP,[0,0.0017,0;1,0.0100,1],[]),imadjust(RGBSTA,[0,0.0017,0;1,0.0100,1],[]),'montage')
    axis tight
    
    noiseAll=abs([montageDataEXP.montageData.NoiseIntsfGFP,montageDataEXP.montageData.NoiseIntCyOFP,...
        montageDataEXP.montageData.NoiseExt,montageDataEXP.montageData.NoiseExt,...
        montageDataEXP.montageData.NoiseTotalsfGFP,montageDataEXP.montageData.NoiseTotalCyOFP,...
        montageDataSTA.montageData.NoiseIntsfGFP,montageDataSTA.montageData.NoiseIntCyOFP,...
        montageDataSTA.montageData.NoiseExt,montageDataSTA.montageData.NoiseExt,...
        montageDataSTA.montageData.NoiseTotalsfGFP,montageDataSTA.montageData.NoiseTotalCyOFP]);
    maxNoise=max(noiseAll);
    %Log period noise
    x=categorical({'Intrinsic','extrinsic','total'});
    noise=[montageDataEXP.montageData.NoiseIntsfGFP,montageDataEXP.montageData.NoiseIntCyOFP;montageDataEXP.montageData.NoiseExt,montageDataEXP.montageData.NoiseExt;montageDataEXP.montageData.NoiseTotalsfGFP,montageDataEXP.montageData.NoiseTotalCyOFP];
    subplot3=subplot(2,2,3);
    hold(subplot3,'on')
    bar(subplot3,x,noise)
    ylim(subplot3,[0,maxNoise+0.1])
    % axis tight
    ylabel('nosie')
    title('exponential phase')
    set(gca,'FontSize',17);
    
    %Stable period noise
    x=categorical({'Intrinsic','extrinsic','total'});
    noise=[montageDataSTA.montageData.NoiseIntsfGFP,montageDataSTA.montageData.NoiseIntCyOFP;montageDataSTA.montageData.NoiseExt,montageDataSTA.montageData.NoiseExt;montageDataSTA.montageData.NoiseTotalsfGFP,montageDataSTA.montageData.NoiseTotalCyOFP];
    subplot4=subplot(2,2,4);
    hold(subplot4,'on')
    bar(subplot4,x,noise)
    ylim(subplot4,[0,maxNoise+0.1])
    % axis tight
    ylabel('nosie')
    title('stationary phase')
    set(gca,'FontSize',17);
    
    annotation(figure1,'textbox',...
        [0.628 0.513 0.0963 0.0471],...
        'String',{'EXP PHASE'},...
        'FontSize',17,...
        'EdgeColor','none');
    
    %  textbox
    annotation(figure1,'textbox',...
        [0.777 0.513 0.0963 0.0471],...
        'String',{'STA PHASE'},...
        'FontSize',17,...
        'EdgeColor','none');
    
    
    suptitle_cwh(figure1,StrainInfo,24)
    saveas(figure1,[DataFile,'\Strain',num2str(iStrain,'%03.f'),'\figure1.tif'])
    saveas(figure1,[DataFile,'\Strain',num2str(iStrain,'%03.f'),'\figure1.fig'])
    close all
    
    
    %%
    %figure 2
    %Logarithmic Period Scatter Plot
    figure2=figure('OuterPosition',[-7 33 1936 1056],'Visible','off');
    subplot1=subplot(2,4,1);
    hold(subplot1,'on');
    [~] = scatplot_ys(montageDataEXP.montageData.CyOFP,montageDataEXP.montageData.sfGFP);
    % axis tight
    xlim(subplot1,[10,80])
    y1=ylim;
    if y1(1)<=0
        ylim(subplot1,[0,inf])
    end
    xlabel('CyOFP(uM)')
    ylabel('sfGFP(uM)')
    title('exponential phase')
    set(gca,'FontSize',17);
    
    subplot2=subplot(2,4,2);
    hold(subplot2,'on');
    plotIndensityDistribution(montageDataEXP.montageData.FLRatio,montageDataEXP.montageData.meanFLRatio,'b');
    axis tight
    xlabel('GR Ratio')
    ylabel('probability density')
    title('exponential phase')
    set(gca,'FontSize',17);
    subplot3=subplot(2,4,3);
    hold(subplot3,'on');
    plotIndensityDistribution(montageDataEXP.montageData.sfGFP,montageDataEXP.montageData.meansfGFP,'c');
    Data1=montageDataEXP.montageData.sfGFP;
    if numel(Data1(Data1<=0))/numel(Data1)<=0.2
        plotgamfit(montageDataEXP.montageData.sfGFP,subplot3);
    end
    clear Data1
    axis tight
    xlabel('sfGFP(uM)')
    ylabel('probability density')
    title('exponential phase')
    set(gca,'FontSize',17);
    %CyOFP分布
    
    subplot4=subplot(2,4,4);
    hold(subplot4,'on');
    plotIndensityDistribution(montageDataEXP.montageData.CyOFP,montageDataEXP.montageData.meanCyOFP,'m');
    plotgamfit(montageDataEXP.montageData.CyOFP,subplot4);
    axis tight
    xlabel('CyOFP(uM)')
    ylabel('probability density')
    title('exponential phase')
    set(gca,'FontSize',17);
    subplot5=subplot(2,4,5);
    hold(subplot5,'on');
    [~] = scatplot_ys(montageDataSTA.montageData.CyOFP,montageDataSTA.montageData.sfGFP);
    xlim(subplot5,[100,900])
    % axis tight
    y1=ylim;
    if y1(1)<=0
        ylim(subplot5,[0,inf])
    end
    xlabel('CyOFP(uM)')
    ylabel('sfGFP(uM)')
    title('stationary phase')
    set(gca,'FontSize',17);
    subplot6=subplot(2,4,6);
    hold(subplot6,'on');
    plotIndensityDistribution(montageDataSTA.montageData.FLRatio,montageDataSTA.montageData.meanFLRatio,'b');
    axis tight
    xlabel('GR Ratio')
    ylabel('probability density')
    title('stationary phase')
    set(gca,'FontSize',17);
    subplot7=subplot(2,4,7);
    hold(subplot7,'on');
    plotIndensityDistribution(montageDataSTA.montageData.sfGFP,montageDataSTA.montageData.meansfGFP,'c');
    Data1=montageDataSTA.montageData.sfGFP;
    if numel(Data1(Data1<=0))/numel(Data1)<=0.2
        plotgamfit(montageDataSTA.montageData.sfGFP,subplot7);
    end
    clear Data1
    axis tight
    xlabel('sfGFP(uM)')
    ylabel('probability density')
    title('stationary phase')
    set(gca,'FontSize',17);
    subplot8=subplot(2,4,8);
    hold(subplot8,'on');
    plotIndensityDistribution(montageDataSTA.montageData.CyOFP,montageDataSTA.montageData.meanCyOFP,'m');
    plotgamfit(montageDataSTA.montageData.CyOFP,subplot8);
    axis tight
    xlabel('CyOFP(uM)')
    ylabel('probability density')
    title('stationary phase')
    set(gca,'FontSize',17);
    
    suptitle_cwh(figure2,StrainInfo,24)
    saveas(figure2,[DataFile,'\Strain',num2str(iStrain,'%03.f'),'\figure2.tif'])
    saveas(figure2,[DataFile,'\Strain',num2str(iStrain,'%03.f'),'\figure2.fig'])
    close all
end
end

function plotIndensityDistribution(Data,meanData,color)
stdData=std(Data);
xData=linspace(meanData-3*stdData,meanData+3*stdData,200);
yData=ksdensity(Data,xData);
plot(xData,yData,'LineStyle','none','Marker','o','MarkerIndices',1:5:length(xData),'MarkerEdgeColor',color)
end

function plotgamfit(Data,subplot)
DataNew=Data(Data>0);
DataNew=sort(DataNew);
phat=gamfit(DataNew);
y=gampdf(DataNew,phat(1),phat(2));
plot(subplot,DataNew,y,'Color','k','lineWidth',2)
text(subplot,0.6,0.8,['a=',num2str(phat(1)),char(13,10)','b=',num2str(phat(2))],'fontSize',15,'Color','k','Units','normalized');
end

function hout=suptitle_cwh(figure,str,fontsize)
%SUPTITLE puts a title above all subplots.
%
%	SUPTITLE('text') adds text to the top of the figure
%	above all subplots (a "super title"). Use this function
%	after all subplot commands.
%
%   SUPTITLE is a helper function for yeastdemo.

%   Copyright 2003-2014 The MathWorks, Inc.


% Warning: If the figure or axis units are non-default, this
% function will temporarily change the units.

% Parameters used to position the supertitle.

% Amount of the figure window devoted to subplots
plotregion = .92;

% Y position of title in normalized coordinates
titleypos  = .95;

% Fontsize for supertitle
fs = get(figure,'defaultaxesfontsize')+4;

% Fudge factor to adjust y spacing between subplots
fudge=1;

haold = gca;
figunits = get(figure,'units');

% Get the (approximate) difference between full height (plot + title
% + xlabel) and bounding rectangle.

if ~strcmp(figunits,'pixels')
    set(figure,'units','pixels');
    pos = get(figure,'position');
    set(figure,'units',figunits);
else
    pos = get(figure,'position');
end
ff = (fs-4)*1.27*5/pos(4)*fudge;

% The 5 here reflects about 3 characters of height below
% an axis and 2 above. 1.27 is pixels per point.

% Determine the bounding rectangle for all the plots

h = findobj(figure,'Type','axes');

oldUnits = get(h, {'Units'});
if ~all(strcmp(oldUnits, 'normalized'))
    % This code is based on normalized units, so we need to temporarily
    % change the axes to normalized units.
    set(h, 'Units', 'normalized');
    cleanup = onCleanup(@()resetUnits(h, oldUnits));
end

max_y=0;
min_y=1;
oldtitle = [];
numAxes = length(h);
thePositions = zeros(numAxes,4);
for i=1:numAxes
    pos=get(h(i),'pos');
    thePositions(i,:) = pos;
    if ~strcmp(get(h(i),'Tag'),'suptitle')
        if pos(2) < min_y
            min_y=pos(2)-ff/5*3;
        end
        if pos(4)+pos(2) > max_y
            max_y=pos(4)+pos(2)+ff/5*2;
        end
    else
        oldtitle = h(i);
    end
end

if max_y > plotregion
    scale = (plotregion-min_y)/(max_y-min_y);
    for i=1:numAxes
        pos = thePositions(i,:);
        pos(2) = (pos(2)-min_y)*scale+min_y;
        pos(4) = pos(4)*scale-(1-scale)*ff/5*3;
        set(h(i),'position',pos);
    end
end

np = get(figure,'nextplot');
set(figure,'nextplot','add');
if ~isempty(oldtitle)
    delete(oldtitle);
end
axes('pos',[0 1 1 1],'visible','off','Tag','suptitle');
ht=text(.5,titleypos-1,str);set(ht,'horizontalalignment','center','fontsize',fontsize,'FontWeight','bold');
set(figure,'nextplot',np);
axes(haold);
if nargout
    hout=ht;
end
end

function resetUnits(h, oldUnits)
% Reset units on axes object. Note that one of these objects could have
% been an old supertitle that has since been deleted.
valid = isgraphics(h);
set(h(valid), {'Units'}, oldUnits(valid));
end


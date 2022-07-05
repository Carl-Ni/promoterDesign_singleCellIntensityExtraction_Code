function Data=PVDfloMeanIntensity(imagefile)
%   Obtain information on the fluorescence intensity of the PVD
%   'imagefile' is  fluorescent images stack of PVD  
imageStack=import_tiff_stack(imagefile);
imageNum=size(imageStack,3);
channelNum=3;
montageNum=9;
kk=channelNum*montageNum;
sampleNum=imageNum/kk;
Data=[];
for i=1:sampleNum
    if sampleNum>1&&i==1
        maskImage=imageStack(:,:,2:channelNum:kk);
    else
        maskImage=imageStack(:,:,kk*(i-1)+3:channelNum:kk*i);
    end
    PVDIimage=imageStack(:,:,kk*(i-1)+1:channelNum:kk*i);
    maskImageNew=myImageProcessingFL(maskImage);
    meanIntensity=[];
    for m=1:montageNum
        subPVDIimage=PVDIimage(:,:,m);
        [L,Num]=bwlabel(maskImageNew(:,:,m));
        for k=1:Num
            meanIntensity=[meanIntensity,mean(subPVDIimage(L==k))];
        end
    end
    Data(i)=mean(meanIntensity);
end
end
        
function [afterProcessingImages,imageProcessingInfo]=myImageProcessingFL(beforeProcessingImages) 
%this function can segreate orignal images as you want and returen the mask images
imageType='uint16'; %here you can chenge your image type
cropInfo=true(size(beforeProcessingImages,1),size(beforeProcessingImages,2));
gaussianFilter=fspecial('gaussian',[5, 5],10); %here create Gaussian Blur Filter
edgeFilter=(ones(5,5)).*-1;edgeFilter(3,3)=24;  %here create edageFilter
afterProcessingImages=zeros(size(beforeProcessingImages),imageType);% here intiatlize the maskImages stacks
grayThresh=800;
areaThreshold=50;
maxIntensity=1000;

imageProcessingInfo.cropInfo=cropInfo;
imageProcessingInfo.grayThresh=grayThresh;
imageProcessingInfo.areaThreshold=areaThreshold;
imageProcessingInfo.maxIntensity=maxIntensity;
for iframe=1:size(beforeProcessingImages,3)
    afterProcessingImages(:,:,iframe)=imfilter(beforeProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),edgeFilter); %use edgeFilter process
    afterProcessingImages(:,:,iframe)=imfilter(afterProcessingImages(:,:,iframe),gaussianFilter); % use Guassian blur filter process
    afterProcessingImages(:,:,iframe)=im2bw(afterProcessingImages(:,:,iframe),grayThresh/(2^16-1));%1.2 for the image with few particles
    afterProcessingImages(:,:,iframe)=bwareaopen(afterProcessingImages(:,:,iframe),areaThreshold); %remove the objecvies which area is smaller than threshold
    cc=regionprops(logical(afterProcessingImages(:,:,iframe)),beforeProcessingImages(:,:,iframe),'MaxIntensity','PixelIdxList');
    for iCC=1:size(cc,1)
        if cc(iCC).MaxIntensity<=maxIntensity
            image=afterProcessingImages(:,:,iframe);
            image(cc(iCC).PixelIdxList)=0;
            afterProcessingImages(:,:,iframe)=image;
        end
    end
    afterProcessingImages(:,:,iframe)=imclearborder(afterProcessingImages(:,:,iframe));
    afterProcessingImages(:,:,iframe)=imfill(afterProcessingImages(:,:,iframe),'holes'); % fill holes process
end
afterProcessingImages=logical(afterProcessingImages);
end



function imageStack= import_tiff_stack( fname )
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

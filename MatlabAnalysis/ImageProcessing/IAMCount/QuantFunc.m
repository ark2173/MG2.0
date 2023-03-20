function [IBA1_Per,IAM_Per] = QuantFunc(NAME)
IBA1_Raw=uint8(imread(NAME,1));
DAPI_Raw=uint8(imread(NAME,2));
IAM_Raw=uint8(imread(NAME,3));

%
DAPI = im2bw(imadjust(DAPI_Raw));
IAM = IAM_Raw-imgaussfilt(IAM_Raw,10);
IAM = im2bw(IAM,0.2);
IBA1 = im2bw(imadjust(IBA1_Raw - imgaussfilt(IBA1_Raw,10)));
imshow([DAPI IBA1 IAM])

% Filtering out large artifacts
IBA1 = bwpropfilt(IBA1,'area',[5,50]);
DAPI = bwpropfilt(DAPI,'area',[3,50]);
IAM = bwpropfilt(IAM,'area',[2,50]);

% Detecting objects
IBA1_Bound = bwboundaries(IBA1);
DAPI_Bound = bwboundaries(DAPI);
IAM_Bound = bwboundaries(IAM);

SZ = size(DAPI)
FullImage=uint8(zeros(SZ(1),SZ(2),3));
FullImage(:,:,1)=IAM_Raw;
FullImage(:,:,2)=IBA1_Raw;
FullImage(:,:,3)=DAPI_Raw;

imshow(FullImage)
hold on
visboundaries(DAPI_Bound,'Color','b')
visboundaries(IBA1_Bound,'Color','g')
visboundaries(IAM_Bound,'Color','r')
hold off

IBA1_Per = length(IBA1_Bound)/length(DAPI_Bound)*100;
IAM_Per = length(IAM_Bound)/length(IBA1_Bound)*100
end
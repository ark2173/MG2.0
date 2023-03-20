function [IAM_Per] = QuantFunc_CAMBCSFB(NAME)
DAPI_Raw=uint8(imread(NAME,1));
IBA1_Raw=uint8(imread(NAME,2));
BCSFB_Raw=uint8(imread(NAME,3));
IAM_Raw=uint8(imread(NAME,4));

%
DAPI = im2bw(imadjust(DAPI_Raw-imgaussfilt(DAPI_Raw,10)));
IAM = IAM_Raw-imgaussfilt(IAM_Raw,10);
IAM = im2bw(imgaussfilt(imadjust(IAM),1),0.3);
IBA1 = im2bw(imadjust(IBA1_Raw - imgaussfilt(IBA1_Raw,10)));
BCSFB = imgaussfilt(BCSFB_Raw-imgaussfilt(BCSFB_Raw,100));
BCSFB = im2bw(imgaussfilt(imadjust(BCSFB),3),0.25);

% Filtering out large artifacts
IBA1 = bwpropfilt(IBA1,'area',[5,50]);
DAPI = bwpropfilt(DAPI,'area',[3,50]);
IAM = bwpropfilt(IAM,'area',[100,9e99]);
BCSFB = bwpropfilt(BCSFB,'area',[200,9e99]);

% applying BCSFB mask to others
IBA1 = IBA1.*BCSFB;
DAPI = DAPI.*BCSFB;
IAM = IAM.*BCSFB;
imshow([DAPI IBA1 IAM BCSFB])

IAM_Count = bwboundaries(IAM);

% IBA1 per BCSFB area
IAM_Per = 100*sum(sum(IAM))/sum(sum(BCSFB));
end
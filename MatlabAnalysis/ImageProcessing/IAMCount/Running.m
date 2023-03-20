clear all
cd('/media/alina/Backup1/RawData/Imaging/VisiumValidation/AD_IAMCount/')
addpath('/media/alina/Backup1/Programming/Matlab/ImageProcessing/IAMCount/')
DIR = dir;

for k=1:length(DIR)
    NAME=DIR(k).name;
    if(~startsWith(NAME,'.'))
        Final(k).name = NAME;
        [Final(k).IBA1_Per, Final(k).IAM_Per] = QuantFunc(NAME);
    end
end
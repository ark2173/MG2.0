clear all
cd('/media/alina/Backup1/RawData/Imaging/VisiumValidation/DAMBCSFB_AD/')
addpath('/media/alina/Backup1/Programming/Matlab/ImageProcessing/DAMBCSFB/')
DIR = dir;

for k=1:length(DIR)
    NAME=DIR(k).name;
    if(~startsWith(NAME,'.'))
        Final(k).name = NAME;
        [Final(k).IAM_Per] = QuantFunc_DAMBCSFB(NAME);
    end
end
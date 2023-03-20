clear all
cd('/media/alina//Papers/MG2.0/Data/ValidationImaging/VisiumValidation/CAMBCSFB_AD/')
addpath('/media/alina/Backup1/Papers/MG2.0/Scripts/MatlabAnalysis/ImageProcessing/CAMBCSFB/')
DIR = dir;

for k=1:length(DIR)
    NAME=DIR(k).name;
    if(~startsWith(NAME,'.'))
        Final(k).name = NAME;
        [Final(k).IAM_Per] = QuantFunc_CAMBCSFB(NAME);
    end
end
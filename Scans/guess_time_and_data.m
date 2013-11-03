function guess_time_and_data(nVec,nBin,nSec,nSims)

dub = 2;
approx_data = dub*nVec*nBin*nSims/1e6;

spm = 60;
approx_time = nSec*nSims/spm;

display(['This scan requires at least ' num2str(approx_data,'%0.2f') ' megabytes of memory and ' num2str(approx_time,'%0.2f') ' minutes to complete.']);
reply = input('Do you want to continue? y/n\n','s');
if strcmp(reply,'n'); error('Too much!'); end
function [ST] = BO_LUT()
%

% Build lookup table of savetags from each day of recording
ST = struct([]);

ST(end+1).recordingDate = '20181206';
ST(end).recordingSession='R01';
ST(end).bankid=1;
ST(end).Ntest=1;
ST(end).eyeID=3;     %% 1 for left eye, 2 for right eye; 3 for both eye
ST(end).Totaldepth=3;
ST(end).isDeep=0;
ST(end).Zero=[];
ST(end).L56=[];
ST(end).L4c=[];
ST(end).L4b=[];
ST(end).L23=[];
ST(end).subject = 'Quicky';
ST(end).probeVersion = 'option4';
ST(end).IncludedSession=0;

ST(end+1).recordingDate = '20181206';
ST(end).recordingSession='R02';
ST(end).bankid=2;
ST(end).Ntest=1;
ST(end).eyeID=3;     %% 1 for left eye, 2 for right eye; 3 for both eye
ST(end).Totaldepth=3;
ST(end).isDeep=0;
ST(end).subject = 'Quicky';
ST(end).probeVersion = 'option4';
ST(end).IncludedSession=0;

ST(end+1).recordingDate = '20181218';
ST(end).recordingSession='R01';
ST(end).bankid=1;
ST(end).Ntest=1;
ST(end).eyeID=3;     %% 1 for left eye, 2 for right eye; 3 for both eye
ST(end).prefeyeID=2;     %% 1 for left eye, 2 for right eye; 3 for both eye
ST(end).Totaldepth=3;
ST(end).isDeep=0;
ST(end).Zero=1070;
ST(end).L56=562;
ST(end).L4c=283;
ST(end).L4b=256;
ST(end).L23=691;
ST(end).ymax=2500;
ST(end).subject = 'Quicky';
ST(end).probeVersion = 'option3';
ST(end).IncludedSession=1;

ST(end+1).recordingDate = '20190124';
ST(end).recordingSession='R01';
ST(end).bankid=1;
ST(end).Ntest=1;
ST(end).eyeID=1;     %% 1 for left eye, 2 for right eye; 3 for both eye
ST(end).Totaldepth=3;
ST(end).isDeep=0;
ST(end).subject = 'Quicky';
ST(end).probeVersion = 'option3';
ST(end).IncludedSession=0;

ST(end+1).recordingDate = '20190305';
ST(end).recordingSession='R01';
ST(end).bankid=1;
ST(end).Ntest=1;
ST(end).eyeID=2;     %% 1 for left eye, 2 for right eye; 3 for both eye
ST(end).prefeyeID=2;     %% 1 for left eye, 2 for right eye; 3 for both eye
ST(end).Totaldepth=3;
ST(end).isDeep=0;

% ST(end).Zero=1070; %% orientation depth
% ST(end).L56=562;
% ST(end).L4c=283;
% ST(end).L4b=256;
% ST(end).L23=691;

ST(end).Zero=1170; %1150 %% bo depth
ST(end).L56=700; %660;
ST(end).L4c=310; %300;
ST(end).L4b=330; %350;
ST(end).L23=650;% 700;

ST(end).ymax=2500;
ST(end).csd='0305r1_CSD';
ST(end).subject = 'Quicky';
ST(end).probeVersion = 'option3';
ST(end).IncludedSession=1;

ST(end+1).recordingDate = '20190305';
ST(end).recordingSession='R02';
ST(end).bankid=1;
ST(end).Ntest=1;
ST(end).eyeID=1;
ST(end).prefeyeID=1;     %% 1 for left eye, 2 for right eye; 3 for both eye
ST(end).Totaldepth=3.5;
ST(end).isDeep=0;

% ST(end).Zero=990; %% orientation paper, orientation depth
% ST(end).L56=637;
% ST(end).L4c=312
% ST(end).L4b=308;
% ST(end).L23=652;

ST(end).Zero=1020; %1020; %% bo depth
ST(end).L56=660; %637
ST(end).L4c=310; %312
ST(end).L4b=350; %350;
ST(end).L23=660; %660;

ST(end).ymax=2500;
ST(end).csd='0305r2_CSD';
ST(end).subject = 'Quicky';
ST(end).probeVersion = 'option3';
ST(end).IncludedSession=1;

ST(end+1).recordingDate = '20190305';
ST(end).recordingSession='R03';
ST(end).bankid=1;
ST(end).Ntest=1;
ST(end).eyeID=2;     %% 1 for left eye, 2 for right eye; 3 for both eye
ST(end).prefeyeID=2;     %% 1 for left eye, 2 for right eye; 3 for both eye
ST(end).Totaldepth=3;
ST(end).isDeep=0;
ST(end).Zero=1380;
ST(end).L56=660; 
ST(end).L4c=300;
ST(end).L4b=450;
ST(end).L23=980;
ST(end).ymax=3200;
ST(end).subject = 'Quicky';
ST(end).probeVersion = 'option3';
ST(end).IncludedSession=1;

ST(end+1).recordingDate = '20190305';
ST(end).recordingSession='R04';
ST(end).bankid=2;
ST(end).Ntest=1;
ST(end).eyeID=3;
ST(end).prefeyeID=1;     %% 1 for left eye, 2 for right eye; 3 for both eye
ST(end).Totaldepth=3;
ST(end).isDeep=0;
ST(end).Zero=1190;
ST(end).L56=488;
ST(end).L4c=242;
ST(end).L4b=316;
ST(end).L23=607;
ST(end).ymax=2500;
ST(end).csd='0305r4_CSD';
ST(end).subject = 'Quicky';
ST(end).probeVersion = 'option3';
ST(end).IncludedSession=0;

ST(end+1).recordingDate = '20190320';
ST(end).recordingSession='R01';
ST(end).bankid=1;
ST(end).Ntest=1;
ST(end).eyeID=3;
ST(end).Totaldepth=3;
ST(end).isDeep=1;
ST(end).subject = 'Lannist';
ST(end).probeVersion = 'option3';
ST(end).IncludedSession=0;

ST(end+1).recordingDate = '20190320';
ST(end).recordingSession='R02';
ST(end).bankid=1;
ST(end).Ntest=1;
ST(end).eyeID=2;
ST(end).prefeyeID=2;     %% 1 for left eye, 2 for right eye; 3 for both eye
ST(end).Totaldepth=3;
ST(end).isDeep=1;
ST(end).Zero=800; %760; %% bo depth
ST(end).L56=400; %400;
ST(end).L4c=260; %250;
ST(end).L4b=360; %350;
ST(end).L23=500; %327

ST(end).ymax=1500;
ST(end).csd='0320r2_CSD';
ST(end).subject = 'Lannist';
ST(end).probeVersion = 'option3';
ST(end).IncludedSession=2;

ST(end+1).recordingDate = '20190320';
ST(end).recordingSession='R03';
ST(end).bankid=1;
ST(end).Ntest=1;
ST(end).eyeID=1;
ST(end).Totaldepth=3;
ST(end).isDeep=0;
ST(end).isDeep=0;
ST(end).Zero=2150;
ST(end).L56=750;
ST(end).L4c=280;
ST(end).L4b=349;
ST(end).L23=700;
ST(end).ymax=2500;
ST(end).subject = 'Lannist';
ST(end).probeVersion = 'option3';
ST(end).IncludedSession=0;

ST(end+1).recordingDate = '20190320';
ST(end).recordingSession='R03';
ST(end).bankid=1;
ST(end).Ntest=1;
ST(end).eyeID=1;
ST(end).Totaldepth=3;
ST(end).isDeep=1;
ST(end).isDeep=1;
ST(end).Zero=750;
ST(end).L56=374;
ST(end).L4c=328;
ST(end).L4b=349;
ST(end).L23=327;
ST(end).ymax=1500;
ST(end).subject = 'Lannist';
ST(end).probeVersion = 'option3';
ST(end).IncludedSession=0;

ST(end+1).recordingDate = '20190320';
ST(end).recordingSession='R04';
ST(end).bankid=1;
ST(end).Ntest=1;
ST(end).eyeID=1;
ST(end).prefeyeID=1;     %% 1 for left eye, 2 for right eye; 3 for both eye
ST(end).Totaldepth=3;
ST(end).isDeep=1;

% ST(end).Zero=550; %% orientation depth
% ST(end).L56=383;
% ST(end).L4c=250;
% ST(end).L4b=350;
% ST(end).L23=217;

ST(end).Zero=640;% 620;
ST(end).L56=420; %400;
ST(end).L4c=280;% 268;
ST(end).L4b=350; %350;
ST(end).L23=500; % 217

ST(end).ymax=1500;
ST(end).csd='0320r4_CSD';
ST(end).subject = 'Lannist';
ST(end).probeVersion = 'option3';
ST(end).IncludedSession=2;

ST(end+1).recordingDate = '20190320';
ST(end).recordingSession='R05';
ST(end).bankid=1;
ST(end).Ntest=1;
ST(end).eyeID=2;
ST(end).Totaldepth=3;
ST(end).isDeep=0;
ST(end).isDeep=0;
ST(end).Zero=2680;
ST(end).L56=674;
ST(end).L4c=280;
ST(end).L4b=340;
ST(end).L23=627;
ST(end).ymax=1500;
ST(end).subject = 'Lannist';
ST(end).probeVersion = 'option3';
ST(end).IncludedSession=0;

% ST(end+1).recordingDate = '20190320';
% ST(end).recordingSession='R05';
% ST(end).bankid=1;
% ST(end).Ntest=1;
% ST(end).eyeID=2;
% ST(end).Totaldepth=3;
% ST(end).isDeep=1;
% ST(end).isDeep=1;
% ST(end).Zero=750;
% ST(end).L56=374;
% ST(end).L4c=328;
% ST(end).L4b=349;
% ST(end).L23=327;
% ST(end).ymax=1500;
% ST(end).subject = 'Lannist';
% ST(end).probeVersion = 'option3';
% ST(end).IncludedSession=0;

ST(end+1).recordingDate = '20191009';
ST(end).recordingSession='R01';
ST(end).bankid=1;
ST(end).Ntest=1;
ST(end).eyeID=1;
ST(end).prefeyeID=1;     %% 1 for left eye, 2 for right eye; 3 for both eye
ST(end).Totaldepth=3.5;
ST(end).isDeep=0;
ST(end).Zero=1630;
ST(end).L56=760;
ST(end).L4c=400;
ST(end).L4b=600;
ST(end).L23=870;
ST(end).ymax=3500;
ST(end).subject = 'Davis';
ST(end).probeVersion = 'option3';
ST(end).IncludedSession=1;

ST(end+1).recordingDate = '20191009';
ST(end).recordingSession='R02';
ST(end).bankid=1;
ST(end).Ntest=1;
ST(end).eyeID=1;
ST(end).prefeyeID=1;     %% 1 for left eye, 2 for right eye; 3 for both eye
ST(end).Totaldepth=3.5;
ST(end).isDeep=0;
ST(end).Zero=1150;
ST(end).L56=1000;
ST(end).L4c=400;
ST(end).L4b=600;
ST(end).L23=1652;
ST(end).ymax=4000;
ST(end).subject = 'Davis';
ST(end).probeVersion = 'option3';
ST(end).IncludedSession=1;

% Create table
ST = struct2table(ST);
ST.index = [1 : size(ST,1)]';


rootpath='/Users/shushu/Dropbox/npix';
save(fullfile(rootpath,'result','bo_lut_based_on_ori_depth.mat'),'ST');
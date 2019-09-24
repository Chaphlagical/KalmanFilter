%
% Reads text file downloaded from
%
% http://navcen.uscg.gov/?pageName=gpsAlmanacs
%
% by clicking on "txt" in the line
%
% Current YUMA Almanac - .alm, .txt
%
% captured in Notepad and saved as a .txt file.
%
% The file is parsed and converted to a m-file containint the essential
% date for use in satellite simulations, which consists of the following
% data fields for each usable satellite 
%
% 1. Orbital Inclination [rad]
% 2. Right Ascension [rad]
% 3. Argument of perigee [rad]
% 4. Mean Anomaly [rad]
% 5. PRN number (ID)
%
% BECAUSE THE DATA FORMATTING FROM THE USCG SITE MAY CHANGE WITH TIME,
% THIS m-FILE PRODUCES MUCH DIAGNOSTIC OUTPUT TO MAKE IT EASIER TO
% MAKE THE APPROPRIATE CHANGES IN THE SCRIPT
%
% M. S.Grewal and A. P. Andrews
% Kaltman Filtering: Theory and Practice Using MATLAB, 4th Edition
% Wiley, 2014
%
%
% INPUT: 
%   The full name of the ASCII input text file downloaded with YUMA text
%        data, without the "txt" extension.
%   The full name of the ASCII output text file to be used as an
%        m-file
%
% OUTPUT: An m-file useable by the simulation programs, with the name
%   "YUMAdata.m"
%
close all;
clear all;
infile  = input('Input file name (without ''.txt'' extension): ','s');
infid   = fopen([infile,'.txt'],'r');
%
disp(['Opened ',[infile,'.txt'],' with id ',num2str(infid)]);
%
outfile  = input('Output file name (without ''.m'' extension): ','s');
outfid   = fopen([outfile,'.m'],'w');
%
disp(['Opened ',[outfile,'.txt'],' with id ',num2str(outfid)]);
%
fprintf(outfid,'%s\n',['% YUMA almanac data from ',[infile,'.txt']]);
fprintf(outfid,'%s\n','% file creation date and time '); 
fprintf(outfid,'%s\n',['% ',num2str(clock)]); 
fprintf(outfid,'%s\n','%'); 
fprintf(outfid,'%s\n',['% data extracted from ',infile,':']); 
fprintf(outfid,'%s\n','% inclination [rad], right ascension [rad], argument of perigee [rad],');
fprintf(outfid,'%s\n','%  mean anomaly [rad], PRNumber'); 
fprintf(outfid,'%s\n','%'); 
fprintf(outfid,'%s\n','% data save to global arrays'); 
fprintf(outfid,'%s\n','% right ascension (RA), argument of perigee + mean anomaly (PA)'); 
fprintf(outfid,'%s\n','%'); 
fprintf(outfid,'%s\n','global RA PA;'); 
fprintf(outfid,'%s\n','Yuma = ['); 
while ~feof(infid),
    l = fgetl(infid);
    disp(['                                   ',l]);
    pause(.1);
    done = 0;
    if length(l)~=0,
        if l(1:3)=='ID:',
            id = [',',num2str(str2num(l(29:30))),';'];
            disp(['ID =           ',id]);
            pause(.1)
        elseif l(1:6)=='Health',
            h  = l(29:31);
            disp(['Health =       ',h]);
            pause(.1);
        elseif l(1:25)=='Orbital Inclination(rad):',
            oi = l(28:40);
            disp(['Orb. Incl. =   ',oi]);
            pause(.1);
        elseif l(1:11)=='Right Ascen',
            ra = l(29:45);
            disp(['Rt. Ascen =    ',ra]);
            pause(.1);
        elseif l(1:19)=='Argument of Perigee',
            ap = l(28:39);
            disp(['Arg. Perigee = ',ap]);
            pause(.1);
        elseif l(1:9)=='Mean Anom',
            ma = l(28:45);
            disp(['Mean Anom. =   ',ma]);
            pause(.1);
        elseif l(1:4)=='week',
            wkno = l(30:32);
            done = 1;
        end;
        %
        if done,
            lout = [oi,',',ra,',',ap,',',ma,id];
            if h=='000',
                fprintf(outfid,'%s\n',lout);
            end;
        end;
    end;
end;
fclose(infid);
fprintf(outfid,'%s\n','%'); 
fprintf(outfid,'%s\n','% Need to insert ] before ; on last line.'); 
fprintf(outfid,'%s\n','%'); 
fprintf(outfid,'%s\n','[rows,cols] = size(Yuma);'); 
fprintf(outfid,'%s\n','for k=1:rows,'); 
fprintf(outfid,'%s\n','    RA(k) = Yuma(k,2);'); 
fprintf(outfid,'%s\n','    PA(k) = Yuma(k,3) + Yuma(k,4);'); 
fprintf(outfid,'%s\n','end;'); 
fprintf(outfid,'%s\n','%'); 
fprintf(outfid,'%s\n',['% Data for GPS week number ',wkno]);
fprintf(outfid,'%s\n','%'); 
%
fclose(outfid);
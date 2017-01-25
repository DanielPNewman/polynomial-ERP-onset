%Script to mark the latency of CPP onset on a single-trial-basis
%Writed by Rafael Abe ('RA'; abe_rafael@hotmail.com) and Dan Newman ('DN';
%dan.newman86@gmail.com). 30/01/2014 

close all
clear all

load NO1_C_resample_noart_erp %load an erp .mat file saved after running "runafew.m"


%CPP electrode:
% ch =[25]; %electrode Pz for brain products
ch =[31]; % electrode Pz for biosemi    %NOTE: may be worth defining
                                        % individualised ROI
                                        % electrodes for CPP
                                        %per participanbt later down the track to take into account
                                        % individual differences in
                                        % cap-fitting etc. Would give
                                        % better CPP signal

Fs=500; %sampling rate (Hz)

ts = -0.376*Fs:1.875*Fs;   % in sample points, the ERP epoch
t = ts*1000/Fs;

for trial=115:size(erp,3)
    if allrespLR(trial)~=0 %if correct response
        
        slope=0; r=0; positive_roots=0; peak=0; fitted_ERP_trace=0; clear counter
        
        ERP_function(trial,:)= polyfit(t,mean(erp(ch,:,trial),1),13);  %RA: Defining 'ERP_function' as the polynomial  that describes the single trial CPP curve
        slope=polyder(ERP_function(trial,:)); %RA: Defining the polynomial that describes the slope behaviour in each time point of the ERP_function curve. That will be the first derivative of ERP_function, as the first derivative gives the tangent of the curve in that time point
        r=roots(slope);  %RA: By solving the polynomial 'slope' it is possible determine the time points in which the 'ERP_function' reaches a Max or a Min. These roots are values of time (the roots can be real or imaginary numbers, we only want real numbers)
        
        for counter=1:length(r)
            if real(r(counter))>-376 && isreal(r(counter)) && real(r(counter))<allRT(trial)*1000/Fs
                peak=peak+1;
                positive_roots(peak)=r(counter);
            end
        end
        
        fitted_ERP_trace=polyval(ERP_function(trial,:),t);
        
        %Plot:
        figure
        plot(t,polyval(ERP_function(trial,:),t),'color','m'); hold on
        plot(t,erp(ch,:,trial),'color','k');
        line(xlim,[0,0],'Color','k');
        line([0,0],ylim,'Color','k');
        line([allRT(trial)*1000/Fs,allRT(trial)*1000/Fs],ylim);
        A=positive_roots>0 & positive_roots<allRT(trial)*1000/Fs-200; positive_roots=positive_roots(A); %DN: CPP onset must be a positive_root after stimulus onset and before RT (here I've stipulated that CPP onset can't be closer than 200ms from RT, but can remove this constraint and it has very little effect on results)
        for counter=1:length(positive_roots)
            if polyval(ERP_function(trial,:),positive_roots(counter))==min(polyval(ERP_function(trial,:),positive_roots)) %DN: CPP onset must be the positive_root at lowest amplitude of 'ERP_function' after stimulus onset and before RT. i.e. we wanted to define CPP onset as the first positive root after stimulus onset for which no proceeding positive root has a lower amplitude
                if   fitted_ERP_trace(t==2.*round(positive_roots(counter)/2))<fitted_ERP_trace(t==(2.*round(positive_roots(counter)/2))+2) &&  fitted_ERP_trace(t==2.*round(positive_roots(counter)/2))<fitted_ERP_trace(t==(2.*round(positive_roots(counter)/2))-2) %DN: Ensures that the CPP onset is marked as the positive_root when amplitude of 'ERP_function' reaches a Min rather than a Max.
                    line([positive_roots(counter),positive_roots(counter)],ylim,'Color','g');
                end
            end
        end
        
        
    end
end



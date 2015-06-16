function electrode_autoconfigure(~,~,use_spikes,custom_panel,handles)
%ELECTRODE_AUTOCONFIGURE   automatically determien electrode configuration
%   ELECTRODE_AUTOCONFIGURE() automatically determines the appropriate electrode configuration, using either spike
%   coincidence or voltage trace correlation.
%
%   Written by Marshall Crumiller
%   email: mcrumiller@gmail.com
%
%   Updates
%     2015-06-03: Created
%-----------------------------------------------------------------------------------------------------------------------
global data_v;
use_spikes=logical(get(use_spikes,'Value'));

% make sure user realizes that this will wipe data
choice = questdlg('Automatic detection of Electrode Groups may take a while, and will wipe the current sorting data. Are you sure you want to continue?',...
      'Change Electrode Configuration','Yes','No','No');

switch(choice)
    case 'Yes'
        % free up some memory for this
        return_button_Callback([],[],handles);
        clearvars -global -except data_v Ts;
        data_v = [];
        
        %xml_file = getappdata(handles.output,'xml_file');
        mat_file = getappdata(handles.output,'mat_file');
        
        % grab channels from all variables named data_*
        names = whos('-file',mat_file,'-regexp','data_(\d)*$');        
        names = {names(:).name};
        result = regexp(names,'data_(?<ch>\d*)','names');
        result = [result{:}];
        [channels,ind] = sort(str2double({result.ch}));
        names=names(ind);
        
        N = length(channels);
        
        % this will be our correlation matrix
        correlations = zeros(N);
        
        % determine spike coincidence for each channel, using a 5-std-dev-threshold
        thresh = 5;
        if(use_spikes)
            load(mat_file,'Ts');
            spikes = cell(N,1);
            h=waitbar(0,'Processing channels...','color',palette.background);
            set(findobj(h,'type','patch'), 'edgecolor',palette.text,'facecolor',palette.header);
            for i = 1:N
                waitbar(i/N,h);
                load(mat_file,names{i});
                eval(sprintf('data_v = %s''; clear %s;',names{i},names{i}));
                spikes{i} = M_threshold_spikes(handles,[thresh;-thresh],0,0);
            end
            delete(h);
            
            % determine the typical coincidence rate for Poisson spikes
            % with similar firing rates
            %duration=max(Ts(end),1e5);
            %firing_rate=mean(cellfun(@length,spikes))/(Ts(end)-Ts(1));
            %fake_spikes=num2cell(sort(rand(2,duration*firing_rate)*duration,2),2);
            
            %c = analyze_coincidence(fake_spikes); c=c(1,2);
            correlations = analyze_coincidence(spikes);
            
            % *** TEMPORARY FOR LGN RECORDING PURPOSES ONLY ***
            
            % generate proper channel list
            % note: temporary for LGN
            ind=cell2mat({[11 12 6 5],[13 14 4 3],[15 16 2 1],[25 26 24 23],[27 28 22 21],[29 39 20 19],[31 32 18 17],[10 9 7 8]});
            deletelocs=false(1,length(ind));
            for i = 1:length(ind)
                if(~any(channels==ind(i)))
                    deletelocs(i)=true;
                end
            end
            ind(deletelocs)=[];
            
            ind2=zeros(1,length(ind));
            for i = 1:length(ind2)
                ind2(i)=find(channels==ind(i));
            end
            
            fig_std;
            imagesc(correlations(ind2,ind2));
            set(gca,'xtick',1:length(ind2),'xticklabel',num2str(channels(ind2)'));
            set(gca,'ytick',1:length(ind2),'yticklabel',num2str(channels(ind2)'));
            
            
        % Calculate correlation pairs using voltage traces
        else
            h=waitbar(0,'Processing channels...','color',palette.background);
            set(findobj(h,'type','patch'), 'edgecolor',palette.text,'facecolor',palette.header);
            total_tests = N*(N-1)/2;
            ind=1;
            for i = 1:N-1
                load(mat_file,names{i});
                for j = i+1:N
                    waitbar(ind/total_tests,h);
                    load(mat_file,names{j});
                    eval(sprintf('correlations(i,j) = corr(%s'',%s'');',names{i},names{j}));
                    ind=ind+1;
                end
            end
            delete(h);
        
        end
        
        % *** TEMPORARY FOR LGN PURPOSES ONLY ***
        ind=cell2mat({[11 12 6 5],[13 14 4 3],[15 16 2 1],[25 26 24 23],[27 28 22 21],[29 39 20 19],[31 32 18 17],[10 9 7 8]});
        deletelocs=false(1,length(ind));
        for i = 1:length(ind)
            if(~any(channels==ind(i)))
                deletelocs(i)=true;
            end
        end
        ind(deletelocs)=[];
        
        ind2=zeros(1,length(ind));
        for i = 1:length(ind2)
            ind2(i)=find(channels==ind(i));
        end
        
        fig_std;
        imagesc(correlations(ind2,ind2));
        set(gca,'xtick',1:length(ind2),'xticklabel',num2str(channels(ind2)'));
        set(gca,'ytick',1:length(ind2),'yticklabel',num2str(channels(ind2)'));
        %
        set(handles.status_label,'String','Bloop!');
        delete(custom_panel);
    case 'No'
    otherwise
end
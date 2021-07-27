function [ dateeU_fig,top,bottom,plotMatrix ] = plotDateeU( dateeU,plot_type,tT_in,obs_data )
%plotDateeU Plots a series of date eU correlations & corresponding tT paths
%   input is a dateeU matrix with date in one column, eU in other column
%   and tT paths with time in one column, temperature in the other. Intent
%   is to be used on model output only



%how many dateeU subsets?
howMany=size(dateeU,2)/3;

%counter figures out how many plots per dateeU subset (i.e. grain sizes)
counter=zeros(howMany,1);

%index vector for finding where eU values "reset" in a given column. Used
%for comparing model output with different annealing models, max is set to
%three for now
keepTrack=zeros(2,1);
plotMatrix=zeros(length(unique(dateeU(:,2))),1);

%sizing parameters for the plots
max_dateeU = max(dateeU);
max_date = max_dateeU(1,1);
for i=1:3:size(max_dateeU)
    if(max_dateeU(1,i)>max_date)
        max_date = max_dateeU(1,i);
    end
end

max_time_temp = max(tT_in);
max_time = max_time_temp(1,1);
max_temp = max_time_temp(1,2);
for i=1:3:size(tT_in)
    if(max_time_temp(1,i)>max_time)
        max_time = max_time_temp(1,i);
    end
    
    if(max_time_temp(1,i+1)>max_temp)
        max_temp = max_time_temp(1,i+1);
    end
end

top_axis_size = [0 max_dateeU(1,2) 0 max_date];
bottom_axis_size = [0 max_time 0 max_temp];

color_options=magma(howMany);
%color_options=viridis(howMany);

if(plot_type=='grain_size')
    for i=1:howMany
        grain_size=dateeU(1,i*3);
    
        if(any(dateeU(:,i*3)>grain_size) && any(dateeU(:,i*3)<grain_size))
            greater=find(dateeU(:,i*3)>grain_size);
            less=find(dateeU(:,i*3)<grain_size);
            if(greater(1,1)<less(1,1))
                x1=dateeU(1:greater(1,1)-1,i*3-1);
                y1=dateeU(1:greater(1,1)-1,i*3-2);
                x2=dateeU(greater(1,1):greater(end,1),i*3-1);
                y2=dateeU(greater(1,1):greater(end,1),i*3-2);
                x3=dateeU(less(1,1):less(end,1),i*3-1);
                y3=dateeU(less(1,1):less(end,1),i*3-2);
                plotMatrix=[plotMatrix x1 y1 x2 y2 x3 y3];
            else
                x1=dateeU(1:less(1,1)-1,i*3-1);
                y1=dateeU(1:less(1,1)-1,i*3-2);
                x2=dateeU(less(1,1):less(end,1),i*3-1);
                y2=dateeU(less(1,1):less(end,1),i*3-2);
                x3=dateeU(greater(1,1):greater(end,1),i*3-1);
                y3=dateeU(greater(1,1):greater(end,1),i*3-2);
                plotMatrix=[plotMatrix x1 y1 x2 y2 x3 y3];
            end
            counter(i,1)=3;
        elseif(any(dateeU(:,i*3)>grain_size) || any(dateeU(:,i*3)<grain_size))
            different=find(dateeU(:,i*3)~=grain_size);            
            x1=dateeU(1:different(1,1)-1,i*3-1);
            y1=dateeU(1:different(1,1)-1,i*3-2);
            x2=dateeU(different(1,1):different(end,1),i*3-1);
            y2=dateeU(different(1,1):different(end,1),i*3-2);
            plotMatrix=[plotMatrix x1 y1 x2 y2];
            counter(i,1)=2;
        else
            [nonZero,~]=find(dateeU,1,'last');
            x1=dateeU(1:nonZero,i*3-1);
            y1=dateeU(1:nonZero,i*3-2);
            plotMatrix=[plotMatrix x1 y1];
            counter(i,1)=1;
        end
    end
    
    dateeU_fig=figure;
    hold on
    top=subplot(2,1,1);
    xlabel(top,'eU (ppm)')
    ylabel(top,'Date (Ma)')
    bottom=subplot(2,1,2);
    xlabel(bottom, 'Time (Ma)')
    tempLabel=['Temperature (' char(176) 'C)'];
    ylabel(bottom, tempLabel)
    position=2; %reset
    bottom.XDir='reverse';
    bottom.YDir='reverse';
    bottom.TickDir='out';
    top.TickDir='out';
    
    
    hold(top, 'on')
    hold(bottom, 'on')
    axis(top,top_axis_size)
    axis(bottom, bottom_axis_size)
    top.XAxis.LineWidth=1;
    top.YAxis.LineWidth=1;
    bottom.XAxis.LineWidth=1;
    bottom.YAxis.LineWidth=1;
    set(bottom,'XTick',0:bottom_axis_size(1,2)/10:bottom_axis_size(1,2))
    set(bottom, 'YTick',0:bottom_axis_size(1,4)/10:bottom_axis_size(1,4))
    set(top,'XTick',0:top_axis_size(1,2)/10:top_axis_size(1,2))
    set(top, 'YTick',0:top_axis_size(1,4)/10:top_axis_size(1,4))
    sz=25;
    if(obs_data~=0)
        err=2*obs_data(:,2);
        s=scatter(top,obs_data(:,3),obs_data(:,1),sz,'filled');
        s.MarkerFaceColor='r';
        s.MarkerEdgeColor='k';
        errorbar(top,obs_data(:,3), obs_data(:,1), err,'CapSize',0, 'LineStyle','none');
    end
    
    for j=1:howMany
        if(counter(j,1)==3)
            p=plot(top,plotMatrix(:,position),plotMatrix(:,position+1),...
                plotMatrix(:,position+2),plotMatrix(:,position+3),'--',...
                plotMatrix(:,position+4),plotMatrix(:,position+5),'--','LineWidth',1);
            set(p,'Color',color_options(j,:))
            position=position+6;
        elseif(counter(j,1)==2)
            p=plot(top,plotMatrix(:,position),plotMatrix(:,position+1),...
                plotMatrix(:,position+2),plotMatrix(:,position+3),'--','LineWidth',1);
            set(p,'Color',color_options(j,:))
            position=position+4;
        else
            p=plot(top,plotMatrix(:,position),plotMatrix(:,position+1),'LineWidth',1);
            set(p,'Color',color_options(j,:))
            position=position+2;
        end
        lastTime=find(tT_in(:,j*3-2),1,'last');
        p=plot(bottom,tT_in(1:lastTime,j*3-2),tT_in(1:lastTime,j*3-1),'LineWidth',1);
        set(p,'Color',color_options(j,:))
    end
    
    hold(top, 'off')
    hold(bottom, 'off')
    
elseif(plot_type=='model_comp')
    for i=1:howMany
        
        tracker=1;
        for k=2:size(dateeU,1)
            if(dateeU(k,i*3-1)<dateeU(k-1,i*3-1) && dateeU(k,i*3-1)~=0)
                keepTrack(tracker,1)=k;
                tracker=tracker+1;
            end
        end
        
        if(keepTrack(2,1)~=0)
            x1=dateeU(1:keepTrack(1,1)-1,i*3-1);
            y1=dateeU(1:keepTrack(1,1)-1,i*3-2);
            x2=dateeU(keepTrack(1,1):keepTrack(2,1)-1,i*3-1);
            y2=dateeU(keepTrack(1,1):keepTrack(2,1)-1,i*3-2);
            x3=dateeU(keepTrack(2,1):end,i*3-1);
            y3=dateeU(keepTrack(2,1):end,i*3-2);
            plotMatrix=[plotMatrix x1 y1 x2 y2 x3 y3];
            counter(i,1)=3;
        elseif(keepTrack(1,1)~=0)
            x1=dateeU(1:keepTrack(1,1)-1,i*3-1);
            y1=dateeU(1:keepTrack(1,1)-1,i*3-2);
            x2=dateeU(keepTrack(1,1):find(dateeU(:,i*3-2),1,'last'),i*3-1);
            y2=dateeU(keepTrack(1,1):find(dateeU(:,i*3-1),1,'last'),i*3-2);
            plotMatrix=[plotMatrix x1 y1 x2 y2];
            counter(i,1)=2;
        else
            x1=dateeU(1:find(dateeU(:,i*3-2),1,'last'),i*3-1);
            y1=dateeU(1:find(dateeU(:,i*3-1),1,'last'),i*3-2);
            plotMatrix=[plotMatrix x1 y1];
            counter(i,1)=1;
        end
        
    end
    
    dateeU_fig=figure;
    top=subplot(2,1,1);
    xlabel(top,'eU (ppm)')
    ylabel(top,'Date (Ma)')
    bottom=subplot(2,1,2);
    xlabel(bottom, 'Time (Ma)')
    tempLabel=['Temperature (' char(176) 'C)'];
    ylabel(bottom, tempLabel)
    position=2; %reset
    bottom.XDir='reverse';
    bottom.YDir='reverse';
    bottom.TickDir='out';
    top.TickDir='out';
    
   
    hold(top, 'on')
    hold(bottom, 'on')
    axis(top,top_axis_size)
    axis(bottom, bottom_axis_size)
    top.XAxis.LineWidth=1;
    top.YAxis.LineWidth=1;
    bottom.XAxis.LineWidth=1;
    bottom.YAxis.LineWidth=1;
    set(bottom,'XTick',0:bottom_axis_size(1,2)/10:bottom_axis_size(1,2))
    set(bottom, 'YTick',0:bottom_axis_size(1,4)/10:bottom_axis_size(1,4))
    set(top,'XTick',0:top_axis_size(1,2)/10:top_axis_size(1,2))
    set(top, 'YTick',0:top_axis_size(1,4)/10:top_axis_size(1,4))
    sz=40;
    if(obs_data~=0)
        err=obs_data(:,2);
        s=scatter(top,obs_data(:,3),obs_data(:,1),sz,'d','filled');
        s.MarkerFaceColor='b';
        errorbar(top,obs_data(:,3), obs_data(:,1), err,'CapSize',0, 'LineStyle','none');
    end
    
    for j=1:howMany
        if(counter(j,1)==3)
            p=plot(top,plotMatrix(:,position),plotMatrix(:,position+1),...
                plotMatrix(:,position+2),plotMatrix(:,position+3),'--',...
                plotMatrix(:,position+4),plotMatrix(:,position+5),':','LineWidth',1);
            set(p,'Color',color_options(j,:))
            position=position+6;
        elseif(counter(j,1)==2)
            p=plot(top,plotMatrix(:,position),plotMatrix(:,position+1),...
                plotMatrix(:,position+2),plotMatrix(:,position+3),'--','LineWidth',1);
            set(p,'Color',color_options(j,:))
            position=position+4;
        else
            p=plot(top,plotMatrix(:,position),plotMatrix(:,position+1),'LineWidth',1);
            set(p,'Color',color_options(j,:))
            position=position+2;
        end
        lastTime=find(tT_in(:,j*3-2),1,'last');
        p=plot(bottom,tT_in(1:lastTime,j*3-2),tT_in(1:lastTime,j*3-1),'LineWidth',1);
        set(p,'Color',color_options(j,:))
    end
    
    hold(top, 'off')
    hold(bottom, 'off')
    
end


end


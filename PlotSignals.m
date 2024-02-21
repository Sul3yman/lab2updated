function PlotSignals( plot_flag,Type,varargin)

if plot_flag
    
    
    if strcmp(Type,'Tx') %Plot Transmitter signals
        
        %Set figure placement for Tx plot
        h1=figure(1);
        scnsize = get(0,'ScreenSize');
        set(gcf,'PaperPositionMode','auto')
        set(h1, 'Position', [0 scnsize(4)-950 600 900])
        
        %Input args for Tx plot
        a=varargin{1};
        s=varargin{2};
        
        %Plot
        subplot(2,1,1)
            stem(a)
            title('Transmitted symbols')
            xlim([1 length(a)+1])
            ylabel('a[k]');
            xlabel('k')
            grid on
        subplot(2,1,2)
            stem(s)
            title('Transmitted samples')
            xlim([1 length(s)])
            ylabel('s[n]');
            xlabel('n')
            grid on
        
    elseif strcmp(Type,'Rx')    %Plot Receiver signals
        
        %Set figure placement for Rx plot
        h2=figure(2);
        scnsize = get(0,'ScreenSize');
        set(gcf,'PaperPositionMode','auto')
        set(h2, 'Position', [600 scnsize(4)-950 600 900])
        
        %Input args for Rx plot
        r=varargin{1};
        y=varargin{2};
        y_sampled=varargin{3};
        
        %Plot
        subplot(3,1,1)
            stem(r,'r')
            title('Received samples')
            xlim([1 length(r)])
            ylabel('r[n]');
            xlabel('n')
            grid on
        subplot(3,1,2)
            stem(y,'r')
            title('Output from MF')
            xlim([1 length(r)])
            ylabel('y[n]');
            xlabel('n')
            grid on
        subplot(3,1,3)
            stem(y_sampled,'r')
            title('Sampled MF')
            xlim([1 length(y_sampled)+1])
            ylabel('y[k]');
            xlabel('k')
            grid on
        
       
    end
    
    
end



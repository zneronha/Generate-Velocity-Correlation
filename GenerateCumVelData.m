function GenerateCumVelData
% Created 18 July 2017 for the purposes of determing velocity correlation
% as a function of distance
%Modified 19 July to check correlation values for close cells
%Modified 19 July to do real plots for sets of cell timeframes
%Further Updated to plot cumulativly
%updated August 1 to stop doing this
%Updated August 1 to plot over time

%VERSION 3.0: 28 DECEMBER 2017: UPDATED TO GENERATE ONE DATASET OVER THE
%WELL FOR IMPLEMENTATION IN CORRELATION FIT CODE

%CURRENT AS OF 26 JANUARY 2018

clearvars
close all
clc

%Ask user for parameters
prompt = 'Enter Wells to be analyzed?\n';
wells = input(prompt);
%define set of frames to be analyzed 
% set = [1:24;25:48;49:72;73:96;97:120;121:144;145:168;169:192;193:216;217:240;241:260];
% set = [1:
set = [1:10:249]'; %#ok<NBRAK>
numit = size(set,1);

cmap = cbrewer('seq','BuPu',size(set,1));


hlmat = nan(numit,numel(wells));
rmat = nan(numit,numel(wells));

tic
wellcounter = 0;

welldatastore = cell(max(wells),numit);
%loop through the wells
for well = wells
    
    
    figure;
    wellcounter = wellcounter + 1;
    Ycolorcounter = 1;
    %load the appropriate data
    filetoload = strcat('Z:\ENG_BBCancer_Shared\group\0Zach\Cluster data\EGF (E6) Data\egf(E6)newwell',...
        num2str(well),'.mat');
    load(filetoload);
    
    Tstoredata = nan(2,2);
    tcounter = 1;
    for itt = 1:numit
        frames = set(itt,:);
        %loop through the frames
        for frame = frames
            fprintf('Beginning frame %d!\n',frame);
            %store the x and y data
            frameX = storeX(:,frame); %#ok<NODEF>
            frameY = storeY(:,frame); %#ok<NODEF>
            frameVx = storevelX(:,frame); %#ok<NODEF>
            frameVy = storevelY(:,frame); %#ok<NODEF>

            %remove nan values
            m = isnan(frameX) == 0;
            frameX = frameX(m);
            frameY = frameY(m);
            frameVx = frameVx(m);
            frameVy = frameVy(m);

            %compile and normalze velocity
            frameVel = [frameVx frameVy];
            speed = sqrt((frameVel(:,1).^2)+(frameVel(:,2).^2));
            frameVel = frameVel./speed;

            %figure out matrix dimensions for preallocation
            counter = 1;
            for jj = 2:numel(frameX)-1
                counter = counter + jj;
            end

            %preallocate memory for speed
            %of the format DISTANCE CORRELATION
            StoreData = nan(counter,2);

            %cell counter
            cellcounter = 0;
            %loop through the cells
            for aa = 1:numel(frameX)-1
                cellposition = [frameX(aa) frameY(aa)];
                firstcellVEL = frameVel(aa,:);
                %compare to all the other cells  
                for bb = aa+1:numel(frameX)
                    %extract positions
                    tempX = frameX(bb);
                    tempY = frameY(bb);
                    secondcellVEL = frameVel(bb,:);
                    %calculate the distance and correlation
                    Dist = sqrt(((tempX-cellposition(1)).^2)+((tempY-cellposition(2)).^2));
                    Cor = dot(firstcellVEL,secondcellVEL);
                    %store the values
                    cellcounter = cellcounter + 1;
                    StoreData(cellcounter,:) = [Dist Cor]; 
                end
            end
            Tstoredata(tcounter:tcounter+size(StoreData,1)-1,:) = StoreData;
            tcounter = tcounter + size(StoreData,1);

        end

        %now bin the data for plotting
%         binsize = 30;
        %set up a matrix to modify
        tempdata = Tstoredata;
        binmat = [10:20:70 110:40:500];
        %set up matrix for plotting
        plotmat = zeros(numel(binmat),2);
        plotcounter = 1;
        
        for binceiling = binmat
            currentdata = tempdata((tempdata(:,1)<binceiling),:);
            dist = mean(currentdata(:,1));
%             dist = binceiling;
%             dist = binceiling-(binsize/2);
            cor = mean(currentdata(:,2));
            plotmat(plotcounter,:) = [dist cor];
            plotcounter = plotcounter + 1;
            tempdata = tempdata((tempdata(:,1)>binceiling),:);
        end
%         [~,~,HL,rsq] = ExponentialFit(plotmat);
%         hlmat(Ycolorcounter,wellcounter) = HL;
%         rmat(Ycolorcounter,wellcounter) = rsq;
        
%         xlim([0 500])
%         ylim([-.15 .6])
        plot(plotmat(:,1),plotmat(:,2),'LineWidth',2,'Color',cmap(Ycolorcounter,:));
        Ycolorcounter = Ycolorcounter + 1;
        hold on
        titlename = strcat('Well',num2str(well),'plot');
        title(titlename);
        
        welldatastore{well,itt} = plotmat;
    end
    
    

    
end
% %     legend('0-6 hrs','6-12 hrs','12-18 hrs','18-24 hrs','24-30 hrs','30-36 hrs','36-42 hrs','42-48 hrs','48-54 hrs','54-60 hrs','60-66 hrs','66-72 hrs');
%     title('EGF Starved Correlation Decay')
%     xlabel('Distance (Microns)');
%     ylabel('Mean Cummulative Correlation Coefficient');
% %     print(gcf,'Control Decay','-depsc','-r2000')
save('EGF(E6)Binned_Correlation_Data.mat','welldatastore');
toc
end

function [const1,const2,HL,rsq]=ExponentialFit(plotmat)
uu = isnan(plotmat(:,1))==0;
distance = plotmat(uu,1);
correlation = plotmat(uu,2);
linearfittype = fittype({'x','1'});

% const1 = -10;
% while const1<0 || const1>100
% [tt,g] = fit(distance,correlation,'exp1');
% const1 = tt.a;
% const2 = tt.b;
% HL = log(2)/(-1*const2);
% rsq = g.rsquare;
% end
% if rsq < 0.5
%     disp('Warning! Low RSQ Value!')
% end
plot(tt)
end
%% GND analysis

%%  Clear and close all

clear all; close all; clc

%% Define parameters

% Degree threshold of grain boundary
GrainAngle = 10;

% Number of threads for calculation
nthreads = 8;

% Poisson's ratio of material
poisson = 0.35;

% Threshold to plot GND type distribution map
threshold = 50;

%% Load EBSD data

[fname, pname] = uigetfile('*.crc', 'Choose ebsd .crc');
cd(pname)

ebsd = EBSD.load([pname '\' fname],'convertEuler2SpatialReferenceFrame','setting 2');
ebsd_raw = ebsd;
CS = ebsd.CS;

% x to east and y to down match the direction of SEM images
plotx2east; plotzIntoPlane

% plot raw map
figure
[~,mP] = plot(ebsd,ebsd.orientations);
print(gcf,[fname(1:end-4) '_RawMap'],'-dpng','-r400');
saveas(gcf,[fname(1:end-4) '_RawMap'],'fig');

mP.micronBar.visible = 'off';
print(gcf,[fname(1:end-4) '_RawMap_NoScaleBar'],'-dpng','-r400');

% plot ipf key
ipfKey = ipfColorKey(ebsd('Magnesium'));
figure;
plot(ipfKey)
print(gcf,[fname(1:end-4) '_ipfKey'],'-dpng','-r400');

%% Optional: choose specific region
% 
% ind_ebsd = find(ebsd.y < 130);
% ebsd = ebsd(ind_ebsd);
% 
% % plot raw map
% figure
% plot(ebsd,ebsd.orientations);
% print(gcf,[fname(1:end-4) '_Raw Chosen Region'],'-dpng','-r400');
% saveas(gcf,[fname(1:end-4) '_Raw Chosen Region'],'fig');

%% Denoise and reconstruct grains (not mandatory)
% Here fill, no grains, indexed pixels (all will be filled) are choosen for reconstruction

% reconstruct the grains
[grains,ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',GrainAngle*degree);

% remove very small grains
ebsd(grains(grains.grainSize<=3)) = [];
% redo grains reconstruction
[grains,ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',GrainAngle*degree);

% smooth and fill missing
%F = splineFilter;
F = halfQuadraticFilter;
%F = infimalConvolutionFilter;
ebsd = smooth(ebsd('indexed'),F,'fill');     % fill all
% ebsd = smooth(ebsd,F,'fill',grains);              % only use indexded part

% redo grains reconstruction
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',GrainAngle*degree);
grains = smooth(grains,5);

%%%% Plot orientation map
figure;
[~,mP] = plot(ebsd,ebsd.orientations);
hold on
plot(grains.boundary,'linewidth',1);
hold off

% Save maps
print(gcf,[fname(1:end-4) '_DenoiseMap'],'-dpng','-r400');
saveas(gcf,[fname(1:end-4) '_DenoiseMap'],'fig');

mP.micronBar.visible = 'off';
print(gcf,[fname(1:end-4) '_DenoiseMap_NoScaleBar'],'-dpng','-r400');

%% Calculate GND

% [disArray, systems]= GND_auto( ebsd, nthreads,poisson,cubicType ) ,

[disArray,systems]=GND_auto(ebsd,nthreads,poisson);
% Change the order of row
systems{1}([1 2 3]) = systems{1}([3 1 2]);
% Sum the number
GND=sum_dislocations(disArray,systems,ebsd);

% We use the geometric mean (geomean) because dislocation density is log-normally distributed.
fprintf(' --- Summary of GND Results, geometrical mean--- \n');
for i = 1:length(GND)
    fprintf('%10s - %15s: %f\n',ebsd(ebsd.phase==GND(i).phase).mineral,GND(i).name,geomean(GND(i).data(GND(i).data>0)));
end

fprintf(' --- Summary of GND Results, arithmetic mean--- \n');
for i = 1:length(GND)
    fprintf('%10s - %15s: %f\n',ebsd(ebsd.phase==GND(i).phase).mineral,GND(i).name,mean(GND(i).data(~isnan(GND(i).data))));
end

%% plot GND with grain boundaries

figure
for i = 1:length(GND)
[~,mP] = plot(ebsd, GND(i).data,'micronbar','on');
hold on
end
plot(grains.boundary,'linewidth',0.8,'linecolor','w')
hold off
mtexColorbar
setColorRange([0,300])

saveas(gcf,[fname(1:end-4) '_GND_map'],'fig');
print(gcf,[fname(1:end-4) '_GND_map'],'-dpng','-r400');

mP.micronBar.visible = 'off';
print(gcf,[fname(1:end-4) '_GND_map_NoScaleBar'],'-dpng','-r400');

%% GND plot ('.')

figure;hold on;
[~,mP] = plot(grains.boundary);

color = jet(length(GND)-1);
for i = 1:length(GND)-1
    ind{i} = find( GND(i).data > threshold);
    x = ebsd.x(ind{i}); y = ebsd.y(ind{i});
    hold on;
    plot(x,y,'.')
end

ii = 1;
for i = 1:length(GND)-1
    if ~isempty(ind{i})
        legendname(ii,1:length(systems{1}(i).name)) = systems{1}(i).name;
        ii = ii+1;
    end
end
legend(legendname,'color','none','edgecolor','none')
text(1.05,0.5,['Threshold = ' num2str(threshold)],'fontsize',12,'units','normalized','color','k','interpreter','latex')
hold off

saveas(gcf,[fname(1:end-4) '_GND_dotmap'],'fig');
print(gcf,[fname(1:end-4) '_GND_dotmap'],'-dpng','-r400');

mP.micronBar.visible = 'off';
print(gcf,[fname(1:end-4) '_GND_dotmap_NoScaleBar'],'-dpng','-r400');

%% Plot GND varies with y (GND calculated based on arithmetic mean)

% Set x range
x_min = 10;
x_max = max(ebsd.x)-10;

% Extract x and y of ebsd map
ind = find(ebsd.x > x_min & ebsd.x < x_max);
y = ebsd.y(ind);
y_unique = unique(y);

% Extract GND vs. y
for i = 1:length(GND)-1
    GND_all = GND(i).data(ind);
    for j = 1:length(y_unique)
        GND_temp = GND_all(find(y == y_unique(j)));
        GND_y{i}(j) = mean(GND_temp(~isnan(GND_temp)));
    end
end

% Plot all lines
figure
for i = 1:6
    hold on
    plot(y_unique,GND_y{i}); %,'o','markersize',4)
end
xlim([0,max(y_unique)*1.05])
xlabel('Distance along y ($\rm\mu m$)','interpreter','latex')
ylabel('Average GND density ($\rm\mu m^{-2}$)','interpreter','latex')
ax.FontName = 'Times New Roman';
hold off
boxaxes
legend(legendname,'color','none','edgecolor','none')
saveas(gcf,[fname(1:end-4) '_GND_along y'],'fig');
print(gcf,[fname(1:end-4) '_GND_along y'],'-dpng','-r400');

% % Compare <a> and <c+a>
% sum of <a>
GND_y{7} = GND_y{1} + GND_y{2} + GND_y{3};
% sum of <c+a>
GND_y{8} = GND_y{4} + GND_y{5} + GND_y{6};
figure
p = plot(y_unique,GND_y{7},y_unique,GND_y{8});
p(2).Color = '#77AC30';
xlim([0,max(y_unique)*1.05])
xlabel('Distance along y ($\rm\mu m$)','interpreter','latex');
ylabel('Average GND density ($\rm\mu m^{-2}$)','interpreter','latex');
ax.FontName = 'Times New Roman';
hold off
boxaxes;
legend({'<a>','<c+a>'},'color','none','edgecolor','none')
saveas(gcf,[fname(1:end-4) '_GND_along_compare y'],'fig');
print(gcf,[fname(1:end-4) '_GND_along_compare y'],'-dpng','-r400');

%% Plot average values with error bar

% Set average width of y
width_y = 3;
x_unique = unique(ebsd.x);
y_unique = unique(ebsd.y);

% Average GND_y and calculate errorbar
num = ceil(length(y_unique)/width_y);
step = y_unique(2) - y_unique(1);

for i = 1:6
%    GND_temp = reshape(GND(i).data,length(x_unique),length(y_unique));
    for j = 1:num
        if j < num
%            [lower, mid, upper] = bootStrapGND(GND_temp(1+(j-1)*width_y:j*width_y,:));
            [lower, mid, upper] = bootStrapGND(GND_y{i}(1+(j-1)*width_y:j*width_y));
            errorbar_low{i}(j) = lower;
            GND_y_ave{i}(j) = mid;
            errorbar_up{i}(j) = upper;
        else
%            [lower, mid, upper] = bootStrapGND(GND_temp(1+(j-1)*width_y:end,:));
            [lower, mid, upper] = bootStrapGND(GND_y{i}(1+(j-1)*width_y:end));
            errorbar_low{i}(j) = lower;
            GND_y_ave{i}(j) = mid;
            errorbar_up{i}(j) = upper;
        end
    end
end

% Plot results
figure
for i = 1:6
    if GND_y_ave{i} == 0
        continue;
    end
    hold on
    errorbar((1+width_y/2:width_y:width_y*num+width_y/2)*step,GND_y_ave{i},GND_y_ave{i}-errorbar_low{i},errorbar_up{i}-GND_y_ave{i},...
        'LineWidth',1.5,'Marker','o','MarkerSize',3);
end
xlim([0,max(y_unique)*1.05])
xlabel('Distance along y ($\rm\mu m$)','interpreter','latex');
ylabel('Average GND density ($\rm\mu m^{-2}$)','interpreter','latex');
ax.FontName = 'Times New Roman';
hold off
boxaxes;
legend(legendname,'color','none','edgecolor','none');
saveas(gcf,[fname(1:end-4) '_GND_ave y width' num2str(width_y)],'fig');
print(gcf,[fname(1:end-4) '_GND_ave y width' num2str(width_y)],'-dpng','-r400');

%% Compare <a> and <c+a> with errorbar

% Average <a>
for i = 7:8
    for j = 1:num
        if j < num
            [lower, mid, upper] = bootStrapGND(GND_y{i}(1+(j-1)*width_y:j*width_y));
            errorbar_low{i}(j) = lower;
            GND_y_ave{i}(j) = mid;
            errorbar_up{i}(j) = upper;
        else
            [lower, mid, upper] = bootStrapGND(GND_y{i}(1+(j-1)*width_y:end));
            errorbar_low{i}(j) = lower;
            GND_y_ave{i}(j) = mid;
            errorbar_up{i}(j) = upper;
        end
    end
end

% Plot results
figure
for i = 7:8
    hold on
    p(i) = errorbar((1+width_y/2:width_y:width_y*num+width_y/2)*step,GND_y_ave{i},GND_y_ave{i}-errorbar_low{i},errorbar_up{i}-GND_y_ave{i},...
        'LineWidth',1.5,'Marker','o','MarkerSize',3);
end
p(8).Color = '#77AC30';
xlim([0,max(y_unique)*1.05])
xlabel('Distance along y ($\rm\mu m$)','interpreter','latex');
ylabel('Average GND density ($\rm\mu m^{-2}$)','interpreter','latex');
ax.FontName = 'Times New Roman';
hold off
boxaxes;
legend({'<a>','<c+a>'},'color','none','edgecolor','none');
saveas(gcf,[fname(1:end-4) '_GND_ave_compare widthy ' num2str(width_y)],'fig');
print(gcf,[fname(1:end-4) '_GND_ave_compare widthy ' num2str(width_y)],'-dpng','-r400');

%% Save results

save([fname(1:end-4) '_GND_result'],'disArray','ebsd','ebsd_raw','fname','GND','systems','grains','width_y','legendname');

% Write results into text file
fid=fopen([fname(1:end-4) '_GND_result.txt'],'w');
fprintf(fid,' --- Summary of GND Results, geometrical mean--- \n');
for i = 1:length(GND)
    fprintf(fid,'%10s - %15s: %f\n',ebsd(ebsd.phase==GND(i).phase).mineral,GND(i).name,geomean(GND(i).data(GND(i).data>0)));
end
fprintf(fid,'\n');
fprintf(fid,' --- Summary of GND Results, arithmetic mean--- \n');
for i = 1:length(GND)
    fprintf(fid,'%10s - %15s: %f\n',ebsd(ebsd.phase==GND(i).phase).mineral,GND(i).name,mean(GND(i).data(~isnan(GND(i).data))));
end
fclose(fid);

%% GND plot (contour)

% color = jet(length(GND)-1);
% 
% figure;
% plot(grains.boundary)
% 
% x = unique(ebsd.x);
% y = unique(ebsd.y);
% [X,Y] = meshgrid(x,y);
% 
% for i = 1:length(GND)-1
%     if ~isnan(mean(GND(i).data(GND(i).data>0)))
%     GND_p = reshape(GND(i).data,length(y),length(x));
%     hold on;
%     contour(X,Y,GND_p,'linecolor',color(i,:),'showtext','on');
%     end
% end
% 
% ii = 1;
% for i = 1:length(GND)-1
%     if ~isnan(geomean(GND(i).data(GND(i).data>0)))
%         legendname(ii,1:length(systems{1}(i).name)) = systems{1}(i).name;
%         ii = ii+1;
%     end
% end
% legend(legendname)
% hold off

%% function 1

function [lower, mid, upper] = bootStrapGND(data)
% bootStrapGND Function
% Return a bootstrapped confidence interval on input data
% Adapted for use with the GND Code developed by Travis Skippon
% Written by Chris Cochrane (Mar. 20, 2018)
%
% Variables:
% alpha - Confidence interval is (100 - alpha)%

data=data(~isnan(data));
alpha=2;

nC = length(data);
if nC == 0
    upper = 0;
    mid = 0;
    lower = 0;
else
    nSelections = 1000;
    iterations = 1000;
    
    caseSampleIndices=randi([1 nC],nSelections,iterations);
    
    bootstrapVals = mean(data(caseSampleIndices),1);
    upper=prctile(bootstrapVals, 100 - alpha / 2);
    lower= prctile(bootstrapVals, alpha / 2);
    mid=mean(data);
end
end

%% function 2

function [disArray, systems]= GND_auto( ebsd, nthreads,poisson,cubicType ) 
  %Calculate GND densities for ebsd map using parallel procesing. 
  %ebsd= ebsd data with grain reconstruction complete (reccomended that the data is smoothed first) 
  %nthreads= number of threads you wish to use for parallel computing 
  %poisson= an array containing the poisson's ratios for each phase (if this 
  %is missing you will be prompted to enter the information) 
  %disArray= array containing GND data 
  %systems= structure containing info about how different dislocation types 
  %are stored in disArray 
  %cubicType = cell array setting either BCC or FCC for any cubic phases, in 
  %order.  Non-cubic phases should NOT be specified.  For example, if phase 1  
  %is BCC, phase 2 is HCP, and phase 3 is FCC, then cubicType={'BCC' 'FCC'} 
  
 
  
 
  %Start keeping track of elapsed time 
  tic; 
  
 
  %Set up parallel pool 
  if nthreads>1 
      pool=gcp; 
  end 
  
 
  if nargin==2 
      poissonDefined=0; 
  end 
  if nargin>2 
      poissonDefined=length(poisson); 
  end 
  if nargin>3 
      numCubics=0; 
      cubicTypeDefined=length(cubicType); 
      for i=unique(ebsd.phase(ebsd.phase>0))' 
          if strcmp(ebsd(ebsd.phase==i).CS.lattice,'cubic') 
              numCubics=numCubics+1; 
          end 
      end 
  end 
  
 
  %cubic counter used for checking that cubictype is defined for each cubic 
  %phase.  Should be set at 1 initially. If there are no cubic phases this 
  %will not be used. 
  cubicCounter=1; 
  
 
  %Allocate memory for temporary array to hold curvature data 
  tempcurve=zeros(length(ebsd),6); 
  tempcurve(:,:)=NaN; 
  
 
       
  tempPoisson=poisson; 
  %automated setup of dislocation types 
  
 
  %loop through all phases except phase 0, which is non-indexed data 
  for i=unique(ebsd.phase(ebsd.phase>0))' 
       
      %For each phase, gridify the ebsd data and get the x and y gradients 
      temp=ebsd(ebsd.phase==i); 
      [temp,newId]=gridify(temp); 
       
      gx=temp.gradientX; 
      gx=gx(temp.phase==i); 
       
      gy=temp.gradientY; 
      gy=gy(temp.phase==i); 
      
       
      %Put the curvature components into a temporary variable 
      tempcurve(ebsd.phase==i,:)=[gx.x gx.y gx.z gy.x gy.y gy.z]; 


       
       
       
      %If user didn't input enough Poissson's ratio values, then ask for them 
      if(length(unique(ebsd.phase(ebsd.phase>0)))>poissonDefined) 
          prompt=sprintf('Number of defined Poissons Ratio less than number of phases in ebsd data.  Enter Poisson Ratio for phase %i (%s) now.\n (To avoid this dialog send an array containing the values for all phases as an input when calling GND code)',i,ebsd(ebsd.phase==i).mineral); 
          name = 'Poissons ratio:'; 
          defaultans = {'0.30'}; 
          input = inputdlg(prompt,name,[1 40],defaultans); 
          poisson(i)=str2double(input{:}); 
      else 
          poisson(i)=tempPoisson(find(unique(ebsd.phase(ebsd.phase>0))'==i)); 
      end 
       
      %If current phase is hexagonal then run HCP Burgers vector setup (see 
      %doBurgersHCP.m) 
      if(strcmp(ebsd(ebsd.phase==i).CS.lattice,'hexagonal')) 
          [f{i}, A{i}, lb{i},systems{i}]=doBurgersHCP(ebsd,i,poisson(i)); 
      end 
       
     %If current phase is cubic, then check CubicTypes and run either BCC or 
      %FCC Burgers vector setup (see doBurgersBCC.m and doBurgersFCC.m) 
      if(strcmp(ebsd(ebsd.phase==i).CS.pointGroup,'m-3m')) 
          if(cubicTypeDefined<cubicCounter) 
          question=sprintf('Cubic phase detected for phase %i (%s).  Is this phase FCC or BCC?',i,ebsd(ebsd.phase==i).mineral); 
          default='BCC'; 
          cubicType{cubicCounter} = questdlg(question,'Specify cubic type','FCC','BCC',default); 
          end 
           
          if(strcmp(cubicType{cubicCounter},'BCC')) 
              fprintf('BCC structure selected for phase %i (%s)\n',i,ebsd(ebsd.phase==i).mineral); 
              [f{i}, A{i}, lb{i},systems{i}]=doBurgersBCC(ebsd,i,poisson(i)); 
          elseif(strcmp(cubicType{cubicCounter},'FCC')) 
              fprintf('FCC structure selected for phase %i (%s)\n',i,ebsd(ebsd.phase==i).mineral);
              [f{i}, A{i}, lb{i},systems{i}]=doBurgersFCC(ebsd,i,poisson(i)); 
          end 
          cubicCounter=cubicCounter+1;     
      end 
       
  end 
  
 
  
 
  %This makes sure that the curvature data calculated from the gridified 
  %ebsd map lines up with the original (non-gridified) map.  This  
  %allows the results to be easily plotted on the original ebsd map 
  [temp,newId]=gridify(ebsd); 
  indx=temp.id2ind(newId); 
  curve=tempcurve(indx,:); 
  
 
  
 
  %************************************************************************** 
  % Minimization 
  %************************************************************************** 
  disp('Minimizing dislocation energy...'); 
  
 
  %initialize variable for holding the dislocation densities of each type of 
  %dislocation for each point in the map. 
  disArray = zeros(size(curve,1),max(cellfun('length',f))); 
  
 
  %Define components of curvature tensor 
  kappa21=-curve(:,2); 
  kappa31=-curve(:,3); 
  
 
  kappa12=-curve(:,4); 
  kappa32=-curve(:,6); 
  
 
  kappa11=-curve(:,1); 
  kappa22=-curve(:,5); 
  
 
  %Set up options for linear programming solver 
  options=optimoptions('linprog','Algorithm','dual-simplex','Display','off'); 
  phase=ebsd.phase; 
  fmaxSize=max(cellfun('length',f)); 
  
 
  %Do short test run on the first 1000*nthreads points to estimate time 
  %required for the full dataset. 
  testrun=1:min(1000*nthreads,length(ebsd)); 
  tic 
  parfor (j=testrun,nthreads) 
      if(phase(j)~=0 && sum(isnan(curve(j,:)))==0) 
          x =linprog(f{phase(j)},A{phase(j)},[kappa11(j) kappa12(j) kappa21(j) kappa22(j) kappa31(j) kappa32(j)],[],[],lb{phase(j)},[],[],options); 
          disArray(j,:) = [x; zeros(fmaxSize-length(f{phase(j)}),1)]; 
      end 
  end 
  estimate=toc*(length(curve)-length(testrun))/length(testrun); 
  
 
  %Ptrint out time estimate in HH:MM:SS format. 
  %fprintf('Estimated time to completion: %s.\n',datestr(estimate/60/60/24,'HH:MM:SS')); 
  dt=datetime('now')+seconds(estimate); 
  fprintf('Analysis should be complete by: %s system clock time.\n',datestr(dt)); 
  

  
 
  
 
  %loop through all points 
  parfor (j=1:size(curve,1),nthreads) 
  %for j=1:size(curve,1) %Replace above line with this one for running 
  %wihtout parallel computing 
  
 
      %Only perform calculations on indexed phases, don't redo calculations 
      %that were done in the test run. 
      if(phase(j)~=0 && max(testrun==j)==0 && sum(isnan(curve(j,:)))==0)         
          x =linprog(f{phase(j)},A{phase(j)},[kappa11(j) kappa12(j) kappa21(j) kappa22(j) kappa31(j) kappa32(j)],[],[],lb{phase(j)},[],[],options); 
          %Place solved dislocation densities for various dislocation types 
          %in disArray 
          disArray(j,:) = [x; zeros(fmaxSize-length(f{phase(j)}),1)]; 
      end 
  end 
  
 
  %Print out system time that the analysis finished at, and total elapsed 
  %time. 
  disp('Minimization complete!'); 
  fprintf('Analysis completed at %s system clock time.\n',datestr(datetime('now'))) 
  fprintf('Total elapsed time was: %s \n',datestr(toc/60/60/24,'HH:MM:SS')) 
end

%% function 3

 function [ f, A, lb,systems] = doBurgersHCP(ebsd,phaseNum,poisson) 
 %% Preamble^M 
 %this code outputs array f (containing energies for all dislocation types) 
 %array A (containing Burgers and Line vector information),  
 %array lb (containing lower bounds for density of each dislocation type -i.e. zero),  
 %Structure systems (containing the possibleslip planes and directions and what family they belong to) 
 
 
 Ntypes=0; 
 
 
 %import the crystal symmetry from the loaded ebsd map.  iterated by calling 
 %function to cover all phases.  To test code type in i=1 first. 
 CS=ebsd(ebsd.phase==phaseNum).CS; 
 
 
 
 
 
 
 %%Get units of ebsd coordinates and set burgers vector unit conversion^M 
 %%appropriately. (units of Miller type objects are in Angstroms by default,^M 
 %%so need to convert from that).^M 
 
 
 if strcmp(ebsd.scanUnit,'nm') 
     unitConversion=1e-1; 
 elseif strcmp(ebsd.scanUnit,'um') 
     unitConversion=1e-4; 
 elseif strcmp(ebsd.scanUnit,'mm') 
     unitConversion=1e-7; 
 elseif strcmp(ebsd.scanUnit,'m') 
     unitConversion=1e-10; 
 else 
     disp('Warning! Units of EBSD scan coordinates not recognized! Assuming scan is in microns.') 
     unitConversion=1e-4; 
 end 
 
 
 
 
 
 
 %% Prism Slip System (Edge)^M 
 %b is burgers vector, n is slip plane normal 
 
 
 b=Miller(1,1,-2,0,CS,'uvw'); 
 n=Miller(1,0,-1,0,CS,'hkl'); 
 
 
 %get all equivilent vectors (b) or planes(n), and number of them (c)  
 %NOTE (c) is irrelevant and will be overwritten. 
 
 
 
 
 b = symmetrise(b); 
 n = symmetrise(n,'antipodal'); 
 
 
 %find the cases where these are perpendicular by taking the dot product. 
 
 
 [r,c] = find(isnull(dot_outer(vector3d(b),vector3d(n)))); 
 
 
 
 
 
 
 
 
 %reload (b) and (n) with only the valid systems  
 b = b(r); 
 n = n(c); 
 
 
 %save the valid slip systems and planes into a structure named 'systems' 
 %for prism, there should be three unique b. 
 systems(1).burgers=b; 
 systems(1).plane=n; 
 systems(1).name='Prism<a>'; 
 
 
 %convert the burgers vector index from Miller-Bravais to Miller (orthagonal) 
 %and store in the bt double 
 
 
 burgers(Ntypes+1:Ntypes+size(b,1))=b*unitConversion; 
 
 
 
 
 
 
 %calculate the line vector which is orthoganal to both the plane normal and the 
 %burgers vector. 
 
 
 for i=1:size(b,1)     
     t=Miller(cross(vector3d(b(i)),vector3d(n(i))),CS); 
     t.dispStyle='uvw'; 
     t=round(t); 
     line(i+Ntypes)=t.normalize; 
 end 
 
 
 
 
 % for i=1:size(b,1) 
 %     v1=[round(b(i).u) round(b(i).v) round(b(i).w)]; 
 %     v2=[round(n(i).h) round(n(i).k) round(n(i).l)]; 
 % %     v1=v1./norm(v1); 
 % %     v2=v2./norm(v2); 
 %     t=cross(v1,v2); 
 %     bt(i+Ntypes,2,:)=t/norm(t); 
 % end 
 
 
 
 
 Ntypes=size(burgers,2); 
 prismTypes=1:Ntypes; 
 
 
 %% Screw Dislocations <a>^M 
 %burgers vector of the screw dislocation, which is parallel to the line 
 %vector for screw 
 b=Miller(1,1,-2,0,CS,'uvw'); 
 
 
 %find all equivalents 
 b = symmetrise(b); 
 
 
 %save the valid slip systems and planes into a structure named 'systems' 
 %for <a> type screw, there should be three unique b. 
 
 
 systems(2).burgers=b; 
 systems(2).plane='screw'; 
 systems(2).name='screw<a>'; 
 
 
 %convert this to Miller indexes 
 burgers(Ntypes+1:Ntypes+size(b,1))=b*unitConversion; 
 line(Ntypes+1:Ntypes+size(b,1))=b.normalize; 
 
 
 
 
 
 
 Ntypes=size(burgers,2); 
 screwATypes=(prismTypes(end)+1):Ntypes; 
 
 
 %% Basal Slip System (Edge)^M 
 %b is burgers vector, n is slip plane normal 
 b=Miller(1,1,-2,0,CS,'uvw'); 
 n=Miller(0,0,0,1,CS,'hkl'); 
 
 
 %get all equivalent vectors (b) or planes(n), and number of them (c)  
 %NOTE (c) is irrelevant and will be overwritten. 
 
 
 b = symmetrise(b); 
 n = symmetrise(n,'antipodal'); 
 
 
 %find the cases where these are perpendicular by taking the dot product (=0). 
 [r,c] = find(isnull(dot_outer(vector3d(b),vector3d(n)))); 
 
 
 %reload (b) and (n) with only the valid systems  
 b = b(r); 
 n = n(c); 
 
 
 %save the valid slip systems and planes into a structure named 'systems' 
 %for basal slip, there should be three unique b. 
 systems(3).burgers=b; 
 systems(3).plane=n; 
 systems(3).name='basal<a>'; 
 
 
 %convert the burgers vector index from Miller-Bravais to Miller (orthagonal) 
 %and store in the bt double 
 burgers(Ntypes+1:Ntypes+size(b,1))=b*unitConversion; 
 
 
 
 
 %calculate the line vector which is orthogonal to both the plane normal and the 
 %burgers vector, and normalise it 
 
 
 for i=1:size(b,1)     
     t=Miller(cross(vector3d(b(i)),vector3d(n(i))),CS); 
     t.dispStyle='uvw'; 
     t=round(t); 
     line(i+Ntypes)=t.normalize; 
 end 
 
 
 
 
 
 
 Ntypes=size(burgers,2); 
 basalTypes=(screwATypes(end)+1):Ntypes; 
 
 
 
 
 
 
 
 
 
 
 
 
 %% Pyramidal Slip System (Edge)^M 
 %b is burgers vector, n is slip plane normal 
 b=Miller(1,1,-2,3,CS,'uvw'); 
 n=Miller(1,0,-1,1,CS,'hkl'); 
 
 
 %get all equivalent vectors (b) or planes(n), and number of them (c)  
 %NOTE (c) is irrelevant and will be overwritten. 
 
 
 
 
 b = symmetrise(b); 
 n = symmetrise(n,'antipodal'); 
 
 
 %find the cases where these are perpendicular by taking the dot product. 
 [r,c] = find(isnull(dot_outer(vector3d(b),vector3d(n)))); 
 
 
 %reload (b) and (n) with only the valid systems  
 b = b(r); 
 n = n(c); 
 
 
 %save the valid slip systems and planes into a structure named 'systems' 
 %for prismatic slip, there should be SIX unique b. 
 systems(4).burgers=b; 
 systems(4).plane=n; 
 systems(4).name='Pyramidal<c+a>'; 
 
 
 %convert the burgers vector index from Miller-Bravais to Miller (orthagonal) 
 %and store in the bt double 
 burgers(Ntypes+1:Ntypes+size(b,1))=b*unitConversion; 
 
 
 
 
 %calculate the line vector which is orthogonal to both the plane normal and the 
 %burgers vector, normalised to 1 
 
 
 for i=1:size(b,1)     
     t=Miller(cross(vector3d(b(i)),vector3d(n(i))),CS); 
     t.dispStyle='uvw'; 
     t=round(t); 
     line(i+Ntypes)=t.normalize; 
 end 
 
 
 %this is a counter for the # of systems 
 Ntypes=size(burgers,2); 
 pyramidalTypes=(basalTypes(end)+1):Ntypes; 
 
 
 
 
 %% Pyramidal Slip System 2 (Edge)^M 
 %b is burgers vector, n is slip plane normal 
 b=Miller(1,1,-2,3,CS,'uvw'); 
 n=Miller(1,1,-2,2,CS,'hkl'); 
 
 
 %get all equivalent vectors (b) or planes(n), and number of them (c)  
 %NOTE (c) is irrelevant and will be overwritten. 
 
 
 
 
 b = symmetrise(b); 
 n = symmetrise(n,'antipodal'); 
 
 
 %find the cases where these are perpendicular by taking the dot product. 
 [r,c] = find(isnull(dot_outer(vector3d(b),vector3d(n)))); 
 
 
 %reload (b) and (n) with only the valid systems  
 b = b(r); 
 n = n(c); 
 
 
 %save the valid slip systems and planes into a structure named 'systems' 
 %for prismatic slip, there should be SIX unique b. 
 systems(5).burgers=b; 
 systems(5).plane=n; 
 systems(5).name='Pyramidal2<c+a>'; 
 
 
 %convert the burgers vector index from Miller-Bravais to Miller (orthagonal) 
 %and store in the bt double 
 burgers(Ntypes+1:Ntypes+size(b,1))=b*unitConversion; 
 
 
 %calculate the line vector which is orthogonal to both the plane normal and the 
 %burgers vector, normalised to 1 
 
 
 
 
 for i=1:size(b,1)     
     t=Miller(cross(vector3d(b(i)),vector3d(n(i))),CS); 
     t.dispStyle='uvw'; 
     t=round(t); 
     line(i+Ntypes)=t.normalize; 
 end 
 
 
 % for i=1:size(b,1) 
 %     v1=[round(b(i).u) round(b(i).v) round(b(i).w)]; 
 %     v2=[round(n(i).h) round(n(i).k) round(n(i).l)]; 
 % %     v1=v1./norm(v1); 
 % %     v2=v2./norm(v2); 
 %     t=cross(v1,v2); 
 %     bt(i+Ntypes,2,:)=t/norm(t); 
 % end 
 
 
 %this is a counter for the # of systems 
 Ntypes=size(burgers,2); 
 pyramidalTypes2=(pyramidalTypes(end)+1):Ntypes; 
 
 
 
 
 
 
 
 
 
 
 %% Screw Dislocations <c+a>^M 
 %burgers vector of the screw dislocation, which is parallel to the line 
 %vector for screw 
 b=Miller(1,1,-2,3,CS,'uvw'); 
 
 
 %find all equivalents 
 b = symmetrise(b); 
 
 
 %save the valid slip systems and planes into a structure named 'systems' 
 %there should be 6 systems 
 systems(6).burgers=b; 
 systems(6).plane='screw'; 
 systems(6).name='screw<c+a>'; 
 
 
 %convert the burgers vector index from Miller-Bravais to Miller (orthagonal) 
 %and store in the bt double 
 burgers(Ntypes+1:Ntypes+size(b,1))=b*unitConversion; 
 line(Ntypes+1:Ntypes+size(b,1))=b.normalize; 
 
 
 
 
 Ntypes=size(burgers,2); 
 screwPyramidalTypes=(pyramidalTypes2(end)+1):Ntypes; 
 %% Calculate curvature matrix^M 
 %assign to temporary variables 
 bs=burgers; 
 ls=line; 
 %Calculate A matrix 
 A = -1*[bs.x.*ls.x-0.5*(bs.x.*ls.x + bs.y.*ls.y + bs.z.*ls.z); ls.x.*bs.y; ls.y.*bs.x; ls.y.*bs.y-0.5*(bs.x.*ls.x + bs.y.*ls.y + bs.z.*ls.z); ls.z.*bs.x; ls.z.*bs.y]; 
 
 
 
 
 
 
 
 
 
 
 %poisson for zr = 0.34 
 %for Mg poisson=0.29 
 %calculate the energy of each dislocation type, normalised to that of a basal dislocation   
 %If in the basal plane,this equals absolute value of a^2 
 %if screw=edge energy*(1-poisson's ratio) 
 f(prismTypes)=1; 
 f(screwATypes)=(1-poisson); 
 f(basalTypes)=1; 
 
 
 magBurgPrism=norm(bs(prismTypes(1))); 
 magBurgPyra=norm(bs(pyramidalTypes(1))); 
 
 
 f(pyramidalTypes)=(magBurgPyra^2)/(magBurgPrism^2); 
 f(pyramidalTypes2)=(magBurgPyra^2)/(magBurgPrism^2); 
 f(screwPyramidalTypes)=(magBurgPyra^2)/(magBurgPrism^2)*(1-poisson); 
 lb=zeros(Ntypes,1); 
 
 
 systems(1).indices=prismTypes; 
 systems(2).indices=screwATypes; 
 systems(3).indices=basalTypes; 
 systems(4).indices=pyramidalTypes; 
 systems(5).indices=pyramidalTypes2; 
 systems(6).indices=screwPyramidalTypes; 
 
 
 
 
 end 

 %% function 4
 
  function [GND] = sum_dislocations(disArray, systems,ebsd) 
 %sum_dislocations collects individual dislocation types and group them into 
 %slip systems for easier plotting/analysis 
 %   Detailed explanation goes here 
 
 
 %disArray = output from GND_auto.m containing a matrix of the dislocation 
 %densities for all dislcation types at each point in an ebsd map 
 
 
 %systems = output from GND_auto.m containing information about the burgers 
 %vectors, slip planes, phases, etc. of the dislocation types in disArray 
 
 
 %GND = a structure containing the information from disArray, arranged 
 %according to the slip systems in systems.  Properties of the structure 
 %include burgers, plane, phase, name, data.  Data contains all the GND 
 %densities for that system, name is the name of the slip system (for HCP 
 %systems), plane and burgers are the slip plane and Burgers vectors, and 
 %phase is the phase # 
 
 
 %For example, to get the slip plane of the first slip system, use 
 %GND(1).plane, and to get the dislocation density of dislocations on that 
 %slip system use GND(1).data 
 
 
 
 
 GND_counter=1; 
 
 
 for i=1:length(systems) 
 
 
 
 
 j=1; 
 while j<=length(systems{i}) 
     if(isfield(systems{i},'name')) 
         GND(GND_counter).name=systems{i}(j).name; 
     end 
     GND(GND_counter).data(ebsd.phase==i)=sum(disArray(ebsd.phase==i,systems{i}(j).indices),2); 
     GND(GND_counter).data(ebsd.phase==0 | ebsd.phase >i)=NaN; 
     GND(GND_counter).burgers=systems{i}(j).burgers; 
     GND(GND_counter).plane=systems{i}(j).plane;    
     GND(GND_counter).phase=i; 
     j=j+1; 
     GND_counter=GND_counter+1; 
 end 
 
 
 GND(GND_counter).name='total'; 
 GND(GND_counter).data(ebsd.phase==i)=sum(disArray(ebsd.phase==i,:),2); 
 GND(GND_counter).data(ebsd.phase~=i)=NaN; 
 GND(GND_counter).burgers=[]; 
 GND(GND_counter).plane=[]; 
 GND(GND_counter).phase=i; 
 GND_counter=GND_counter+1; 
 end 
 
 
  end 

 %% function 5
 function [a,b]=boxaxes
% get handle to current axes
% a = gca;
% % % get the original xlim, ylim
% % xlim_ori = get(a,'XLim');
% % ylim_ori = get(a,'YLim');
% % set box property to off and remove background color
% set(a,'box','off','color','none')
% % create new, empty axes with box but without ticks
% b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'xlim',get(a,'XLim'),'ylim',get(a,'YLim'));
% % set original axes as active
% axes(a)
% % link axes in case of zooming
% linkaxes([a b])
box off;
a = gca;
x = get(a,'XLim');
y = get(a,'YLim');
line([x(1) x(2)],[y(2) y(2)],'Color','k','linewidth',0.54)
line([x(2) x(2)],[y(1) y(2)],'Color','k','linewidth',0.54)
xlim(x)
ylim(y)
end
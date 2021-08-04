%DAAM version 1.0, optimized for taking in either apatite or zircon data 
%and running RDAAM or ZRDAAM, calls functions for annealing from 
%Ginster et al., 2019, Guenthner et al., 2013 and Flowers et al., 2009
%contains a modified input scheme for inputing a large time
%temperature matrix from excel (empty column between each individual t-T
%path) for forward modeling, space left for inverse modeling in subsequent
%versions

%decay constants in 1/yr
lambda235=9.8485e-10;
lambda238=1.55125e-10;
lambda232=4.9475e-11;
lambdaf=.0000000000000000846; %for Flowers et al. 2009 RDAAM
lambdaD=.000000000155125; %for Flowers et al. 2009 RDAAM

%functional form, amorphous fraction, mean distance to track parameters
El=165; %kJ/mol
D0l=193188; %cm2/s
D0N17=.0034; %cm2/s
EaN17=71; %kJ/mol
R=.008314472; %kJ/(K*mol)
Ba=5.48E-19; %g amorphized per alpha event
interconnection=3; %fudge factor in DI
SV=1.669; %nm^-1 track surface to volume ratio
lint_lattice=45920; %nm, extrapolated to a zircon with 1e14 alphas/g

%Flowers et al. 2009 RDAAM diffusion constants
psi=1e-13;
omega=1e-22;
D0L=0.6071; %cm2/s calculated from Flowers et al. 2009 Table 1 value and "typical" 60 micron radius
EsubL=122.3; %kJ/mol
EsubTrap=34; %kJ/mol
etaq=0.91; %Durango
L=0.000815; %cm, half the total etchable fission-track length

isGuenthner=false;
isGinster=false;
isnoAnneal=false;

%get file information, eventually I want to just make this into a dialog
%box
concen=input('Input concentration file>');
phase_input=input('Which phase are you running? Enter "ap" for apatite, "zirc" for zircon>');
if(phase_input=="ap")
    obs_data=input('Input measured data file (if none, enter 0)>');
    compType=2;
elseif(phase_input=="zirc")
    obs_data=input('Input measured data file (if none, enter 0)>');
    compType=input('Input 1 for running a model comparison test, 2 for grain size comparison>');
else
    disp('You must enter either "ap" or "zirc"!')
    return
end

NumGrains=0;

if(compType==1)
    meanGrain=input('What is the grain size (in microns)?>');
    if(meanGrain<=0)
        disp('You must pick a grain size greater than 0!')
        return
    end
    annealGuenthner=input('Do you want to run the Guenthner et al. 2013 annealing model? 1 for yes>');
    annealGinster=input('Do you want to run the Ginster et al., 2019 annealing model? 1 for yes>');
    noAnneal=input('Do you want to run a model with no annealing? 1 for yes>');
    anneal_model_counter=0;
    if(annealGuenthner==1) 
        NumGrains=NumGrains+size(concen,1);
        isGuenthner=true;
        anneal_model_counter=anneal_model_counter+1;
    end
    
    if(annealGinster==1)
        NumGrains=NumGrains+size(concen,1);
        isGinster=true;
        anneal_model_counter=anneal_model_counter+1;
    end
    
    if(noAnneal==1)
        NumGrains=NumGrains+size(concen,1);
        isnoAnneal=true;
        anneal_model_counter=anneal_model_counter+1;
    end
    
    if(anneal_model_counter<=1)
        disp('You must enter at least two annealing models!')
        return
    end
    
    model_parameters=cell(NumGrains,5);
    if(NumGrains==size(concen,1)*2)
        model_parameters(:,1)={"zirc"};
        model_parameters(:,2)=num2cell([concen(:,1); concen(:,1)]);
        model_parameters(:,3)=num2cell([concen(:,2); concen(:,2)]);
        model_parameters(:,4)={meanGrain}; 
    
        if(isGinster && isGuenthner)
            model_parameters(1:size(concen,1),5)={'Guenthner'};
            model_parameters(size(concen,1)+1:end,5)={'Ginster'};
        elseif(isGinster && isnoAnneal)
            model_parameters(1:size(concen,1),5)={'Ginster'};
            model_parameters(size(concen,1)+1:end,5)={'none'};
        else
            model_parameters(1:size(concen,1),5)={'Guenthner'};
            model_parameters(size(concen,1)+1:end,5)={'none'};
        end
    
    else
        model_parameters(:,1)={"zirc"};
        model_parameters(:,2)=num2cell([concen(:,1); concen(:,1); concen(:,1)]);
        model_parameters(:,3)=num2cell([concen(:,2); concen(:,2); concen(:,2)]);
    
        if(compType==1)
            model_parameters(1:end,4)={meanGrain};
            model_parameters(1:size(concen,1),5)={'Guenthner'};
            model_parameters(size(concen,1)+1:size(concen,1)*2,5)={'Ginster'};
            model_parameters(size(concen,1)*2+1:end,5)={'none'};
        else
            model_parameters(1:size(concen,1),4)={meanGrain};
            model_parameters(size(concen,1)+1:size(concen,1)*2,4)={meanGrain+stdGrain};
            model_parameters(size(concen,1)*2+1:end,4)={meanGrain-stdGrain};
            if(modelType==1)
                model_parameters(:,5)={'Guenthner'};
            elseif(modelType==2)
                model_parameters(:,5)={'Ginster'};
            else
                model_parameters(:,5)={'none'};
            end
        end
    
    end
    
    plotType='model_comp';
elseif(compType==2)
    if(obs_data~=0)
        meanGrain=mean(obs_data(:,4));
        stdGrain=std(obs_data(:,4));
    else
        meanGrain=input('What is the mean grain size (in microns)?>');
        stdGrain=input('What is the standard deviation in grain size (in microns)?>');
        if(meanGrain<=0 || meanGrain-stdGrain<=0)
            disp('You cannot have any grain sizes less than 0!')
            return
        end
    end
    NumGrains=size(concen,1);
    if(stdGrain>0)
        NumGrains=NumGrains*3;
    end
    
    model_parameters=cell(NumGrains,5);
    if(phase_input=="zirc")
        modelType=input('What annealing model do you want? 0=no annealing, 1=Guenthner et al. 2013, 2=Ginster et al. 2019>');
        model_parameters(:,1)={"zirc"};
        if(stdGrain>0)
            model_parameters(:,2)=num2cell([concen(:,1); concen(:,1); concen(:,1)]);
            model_parameters(:,3)=num2cell([concen(:,2); concen(:,2); concen(:,2)]);
            model_parameters(1:size(concen,1),4)={meanGrain};
            model_parameters(size(concen,1)+1:size(concen,1)*2,4)={meanGrain+stdGrain};
            model_parameters(size(concen,1)*2+1:end,4)={meanGrain-stdGrain};
        else
            model_parameters(:,2)=num2cell(concen(:,1));
            model_parameters(:,3)=num2cell(concen(:,2));
            model_parameters(:,4)={meanGrain};
        end
        
        if(modelType==2)
            model_parameters(:,5)={'Ginster'};
        elseif(modelType==1)
            model_parameters(:,5)={'Guenthner'};
        elseif(modelType==0)
            model_parameters(:,5)={'none'};
        else
            disp('You must pick a valid annealing model!')
            return
        end
    end
    
    if(phase_input=="ap")
        model_parameters(:,1)={"ap"};
        if(stdGrain>0)
            model_parameters(:,2)=num2cell([concen(:,1); concen(:,1); concen(:,1)]);
            model_parameters(:,3)=num2cell([concen(:,2); concen(:,2); concen(:,2)]);
            model_parameters(1:size(concen,1),4)={meanGrain};
            model_parameters(size(concen,1)+1:size(concen,1)*2,4)={meanGrain+stdGrain};
            model_parameters(size(concen,1)*2+1:end,4)={meanGrain-stdGrain};
        else
            model_parameters(:,2)=num2cell(concen(:,1));
            model_parameters(:,3)=num2cell(concen(:,2));
            model_parameters(:,4)={meanGrain};
        end
        model_parameters(:,5)={'Flowers'};
    end
    plotType='grain_size';
else
    disp('You must pick a valid comparison type!')
    return
end

temp_tT=input('Input time-temperature matrix>');
tic


tTsize=size(temp_tT,2);
howMany=(tTsize+1)/3;

%set some variables and allocate space for matrices used throughout
nodes=513;
dateeU=zeros(NumGrains,3);
output_date_eU=zeros(NumGrains,3*howMany);
main_index=1;

for each_grain=1:howMany
    if (each_grain==1)
        for index=2:size(temp_tT,1)
            if (temp_tT(index,1)~=0)
                location=index;
            else
                break
            end
        end
        tTInput=zeros(location,2);
        tTInput(:,1)=temp_tT(1:location,1);
        tTInput(:,2)=temp_tT(1:location,2);
    else
        for index=2:size(temp_tT,1)
            if (temp_tT(index,1+((each_grain-1)*3))~=0)
                location=index;
            else
                break
            end
        end
        tTInput=zeros(location,2);
        tTInput(:,1)=temp_tT(1:location,1+((each_grain-1)*3));
        tTInput(:,2)=temp_tT(1:location,2+((each_grain-1)*3));
    end
    
%create tT path with tT path interpolater
%allocate matrices for diffusion and annealing loops
tT=tTInterpolate(tTInput); %in degrees C
steps=size(tT,1);
alphai=zeros(steps,1);
rhov=zeros(steps,1);
diffusivities=zeros(steps,1);
compareDiff=zeros(steps,NumGrains);
compareDam=zeros(steps,NumGrains);

%now calculate damage array, only need to do this once per tT path IF YOU 
%SELECTED GUENTHNER OR FLOWERS and then calculate total damage within the 
%grain loop below... exception is the Ginster model which is grain
%chemistry specific. Search the input array for instances of "Guenthner"
%and "Flowers" or don't execute either. I kept this here to keep the
%ability to run multiple annealing models with the same file but slightly
%increase the efficiency of the annealing code (fewer calls if Guenthner or
%Flowers has more than one instance)
Guenthner=0;
Flowers=0;
modelName=model_parameters(:,5);
for i=1:size(modelName,1)
    if(strcmp(modelName{i,1},'Guenthner'))
        Guenthner=Guenthner+1;
    elseif(strcmp(modelName{i,1},'Flowers'))
        Flowers=Flowers+1;
    end
end
            
if(Guenthner>0)
    ZircDamage=Guenthner2013(steps,tT);
    annealedDam=zeros(steps,1);
elseif(Flowers>0)
    ApDamage=Flowers2009(steps,tT);
    annealedDam=zeros(steps,1);
end
   

for i=1:NumGrains
    phase=model_parameters{i,1};
    annealModel=model_parameters{i,5};
    Uppm=model_parameters{i,2};    
    Thppm=model_parameters{i,3};
    eU=Uppm+.235*Thppm;
    radius=model_parameters{i,4}; %microns
    U235atom=((Uppm/1000000)*.007204/238.02891)*(6.02214179E23); %atoms/g
	U238atom=((Uppm/1000000)*.992745/238.02891)*(6.02214179E23); %atoms/g
	Thatom=((Thppm/1000000)/232.03805)*(6.02214179E23); %atoms/g  
    alphai(:,1)=0;%reset
    rhov(:,1)=0; %reset
    
    %calculating diffusivities, note that temperature is averaged over the
    %duration of each time step so diffusivity is therefore averaged over
    %the same duration
    if(strcmp(phase,'zirc'))
        %unannealed damage array: damage created at each tT step, needs to be
        %created first to use in assessments of damage level below in
        %Ginster model
        for j=steps-1:-1:1
            t2=tT(j+1,1)*1000000; %in years
            t1=tT(j,1)*1000000; %in years
            alphai(j,1)=(8*U238atom*(exp(lambda238*t1)-exp(lambda238*t2)))+...
            (7*U235atom*(exp(lambda235*t1)-exp(lambda235*t2)))+...
            (6*Thatom*(exp(lambda232*t1)-exp(lambda232*t2)));
        end
        
        %create damage array, annealedDam vector that gets passed to
        %diffusion model, and alpha dose production vector
        %Ginster needs the alphai vector to be passed to it, the other
        %don't need that and annealing damage array is calculated up above
        %(i.e. outside of grain 
        if(strcmp(annealModel,'Ginster'))
            annealedDam=Ginster2019(steps,alphai,tT);
        elseif(strcmp(annealModel,'Guenthner'))
            annealedDam(:,1)=0; %reset
            for tsNode=steps-1:-1:1
                for node=tsNode:-1:1
                    annealedDam(tsNode,1)=annealedDam(tsNode,1)+alphai(node,1)*ZircDamage(tsNode,node);
                end
            end
        else  %aka no annealing
            annealedDam=zeros(steps,1);
            for tsNode=steps-1:-1:1
                for node=tsNode:-1:1
                    annealedDam(tsNode,1)=annealedDam(tsNode,1)+alphai(node,1);
                end
            end
        end
    compareDam(:,i)=annealedDam;
        for j=1:steps-1
            temp=(tT(j,2)+273.15+tT(j+1,2)+273.15)/2;
            damage_amount=annealedDam(j,1); %atoms/g
            fa=1-exp(-Ba*damage_amount);
            DI=1-exp(-Ba*damage_amount*interconnection);
            lint=4.2/(fa*SV)-2.5;
            tort=(lint_lattice/lint)^2;
            %the back and forth between cm2/s and 1/s for diffusivities is
            %necessary to have the functional form of equation 8 from
            %Guenthner et al. 2013
            Dtort=(1/tort)*D0l*exp(-El/(R*temp)); %cm2/s
            Dtorta2=Dtort/(radius*(1/10000)*(1-DI))^2; %1/s
            DN17=D0N17*exp(-EaN17/(R*temp)); %cm2/s
            DN17a2=DN17/(radius*(1/10000)*DI)^2; %1/s
            diffusivities(j,1)=(DI/DN17a2+(1-DI)/Dtorta2)^-1; %cm2/s
            diffusivities(j,1)=diffusivities(j,1)*radius^2; %microns2/s
        end
    elseif(strcmp(phase,'ap'))
        %create damage array, annealedDam vector that gets passed to
        %diffusion model
        
        for j=steps-1:-1:1
            t2=tT(j+1,1)*1000000; %in years
            t1=tT(j,1)*1000000; %in years
            %convert to atoms/volume with apatite density of 3.19 g/cm3 and
            %then make it proportational to each isotope ala Flowers et al.
            %2009
            rhov(j,1)=((8/8)*(U238atom*3.19)*(exp(lambda238*t1)-exp(lambda238*t2)))+...
            ((7/8)*(U235atom*3.19)*(exp(lambda235*t1)-exp(lambda235*t2)))+...
            ((6/8)*(Thatom*3.19)*(exp(lambda232*t1)-exp(lambda232*t2)));
        end
        
    
        annealedDam(:,1)=0; %reset
        for tsNode=steps-1:-1:1
            for node=tsNode:-1:1
                annealedDam(tsNode,1)=annealedDam(tsNode,1)+(lambdaf/lambdaD)*rhov(node,1)*etaq*L*ApDamage(tsNode,node);
            end
        end
       
        
        for j=1:steps-1
            temp=(tT(j,2)+273.15+tT(j+1,2)+273.15)/2;
            damage_amount=annealedDam(j,1); %tracks/cm3, aka erhos
            trapDiff=psi*damage_amount+omega*damage_amount^3;
            diffusivities(j,1)=(D0L*exp(-EsubL/(R*temp)))/(trapDiff*exp(EsubTrap/(R*temp))+1); %cm2/s
        end
        diffusivities(:,1)=diffusivities(:,1)*1e8; %microns2/s
    else
        diffusivities(:,1)=NaN; %so you'll know if you screwed up
    end
    compareDiff(:,i)=diffusivities;
    %pass diffusivities to diffusion model 
    dateeU(i,1)=HeDiff(phase,radius,nodes,diffusivities,tT,U235atom,U238atom,Thatom);
    dateeU(i,2)=eU;
    dateeU(i,3)=radius;
   
end

output_date_eU(:,main_index:main_index+2)=dateeU(:,:);
main_index=each_grain*3+1;

end

[dateeU_fig,top,bottom,plotMatrix]=plotDateeU(output_date_eU,plotType,temp_tT,obs_data);

toc
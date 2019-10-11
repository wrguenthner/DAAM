function [ tT_out ] = tTInterpolate( tTInput )
%tT PATH INTERPOLATER
%function takes in a simple tT path with a few points and returns more
%complicated tT path with many points via an interpolation algorithm, this
%particular interpolator was adapted from c++ code provided by Rich Ketcham 
%(UT Austin). Performs faster than interp1() function when output is fed
%into the diffusion code. Mainly because time slices are much bigger when
%almost no change, or little change occurs across time steps, so diffusion
%code doesn't have to work as hard for these particular time slices.

%initialize tT variables, greater precision is acheived by decreasing these
%variables at the cost of increased comp. time
maxRateAccel=3; %max step-to-step increase in duration
maxTempStep=1; %max allowed temperature step in degrees C

%initialize the tT path and perallocate the space so you don't have to add
%to it with each loop iteration... but make it really big so that you don't
%run out of space (extra zeros get cut off at the end of the funcion

tT_out=zeros(100000,2); %in degrees C
tT_out(1,1)=tTInput(end,1);
tT_out(1,2)=tTInput(end,2);
startTime=tTInput(end,1);
%initialize tT_position, which increments by one within the while loop in
%order to keep track of where you are within the tT_out array
tT_position=1;
 
nextPoint=zeros(1,2);

defTimeStep=startTime*.01; %default time step is 1% of total duration
prevTimeStep=defTimeStep;

%for loop steps through each span of tTInput
for dN=size(tTInput,1):-1:2
    currMaxTimeStep=(tTInput(dN,1)-tTInput(dN-1,1))*.2; %minimum of 5 steps
    %rate of temp change (K/m.y.)
    rate=(tTInput(dN,2)-tTInput(dN-1,2))/(tTInput(dN,1)-tTInput(dN-1,1));
    absRate=abs(rate);
    tempPerTimeStep=absRate*defTimeStep; %temp change per default time step (K)
    if(tempPerTimeStep<=maxTempStep)
        currDefTimeStep=defTimeStep; %default time step for current path segment (m.y.)
    else
        currDefTimeStep=maxTempStep/absRate;
    end
    
    if(currDefTimeStep>currMaxTimeStep)
        currDefTimeStep=currMaxTimeStep;
    end
    
    %make sure that time step is large enough to register
    if (currDefTimeStep<tTInput(dN,1)*1e-14)
        currDefTimeStep = tTInput(dN,1)*1e-14;
    end
    
    endTemp=tTInput(dN-1,2); %temperature at end of current tT segment (K)
    
    %main while loop that adds points to tTPath
    while(tT_out(tT_position,1)>tTInput(dN-1,1))
        timeStep=currDefTimeStep; %size of individual time step (m.y.)
        if(timeStep>prevTimeStep*maxRateAccel)
            timeStep=prevTimeStep*maxRateAccel;
        end
        
        %use nextPoint array to create the next tT point
        %check to see if this is final step for the segment
        %small factor is added in case of a roundoff
        if(timeStep*1.01>tT_out(tT_position,1)-tTInput(dN-1,1))
            nextPoint(1,1)=tTInput(dN-1,1);
            nextPoint(1,2)=endTemp;
        else
            nextPoint(1,1)=tT_out(tT_position,1)-timeStep;
            nextPoint(1,2)=tT_out(tT_position,2)-rate*timeStep;
        end
        
        %increment the current poisition variable by 1, put nextPoint array 
        %point into tT_out array, 
        tT_position=tT_position+1;
        tT_out(tT_position,1)=nextPoint(1,1);
        tT_out(tT_position,2)=nextPoint(1,2);
        
        prevTimeStep=timeStep;
    end
    
end

%remove the end of the array, which is all zeros because of pre-allocation
tT_out(tT_position+1:end,:)=[];

end


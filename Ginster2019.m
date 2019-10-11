function [annealedDam] = Ginster2019(steps, alphai, tT)
    %damage annealing model based on Ginster et al. 2019 kinetics
    
    %requires inputs of tT() matrix and alphai() or "unannealed" damage
    %matrix, returns annealedDam matrix
        
    %create damage array, annealedDam vector that gets passed to
    %diffusion model, and alpha dose production vector
    
    damage=zeros(steps,steps);
    annealer=zeros(steps,1); %doesn't really need to be a vector, but useful for error checking
    annealedDam=zeros(steps,1);
    memoryDam=0;
    
    %annealed damage array, arrangment is damage(1,1) = oldest, oldest tT
    %step; damage(step, steps) = youngest, youngest tT step. Annealing steps
    %forward in time (old -> young) in terms of rows (outer for loop)in
    %order to calculate damage level of previous time step, but backwards
    %in time (young -> old) in terms of columns within each row (inner for
    %loop) in order to maintain consitency with legacy diffusion code from
    %ZRDAAM_v1. As Duddy et al. (1988) demonstrated, the order in which 
    %equivalent time is calculated DOES NOT MATTER. Matrix is arranged such 
    %that the rows describe tT path for a given damage zone, and the 
    %columns give the time of damage creation.
    
    %NOTE: damage matrix convention is 1 = complete annealing, 0 = no
    %annealing as per Ginster et al. derivations, BUT need to convert to 
    %1 = no annealing, 0 = complete annealing in order to keep consistent 
    %with legacy code from ZFT annealing framework
    %This conversion is done by 1-damage(:,:) below
    
    for tsNode=1:steps-1
        teq=0; %reset for each new row
        temp=(tT(tsNode,2)+273.15+tT(tsNode+1,2)+273.15)/2; %temperature between current and next oldest node 
        
        %get damage ingrowth for each step from alphai vector and combine
        %with annealer vector and damage matrix from next oldest time step
        %(i.e.row), number is needed to decide on correct annealing 
        %equation to use, first look-up in damage loop will be a zero (1-0)
        %, okay because we want to add the damage created in the current
        
          if(tsNode>1) %can't look backwards for first step
              for j=tsNode:-1:1 
                  %sum the damage in a given row of damage matrix from 
                  %NEXT OLDEST time step row to the row you're currently in
                  %, i.e. tsNode-1, convert to 1 = no annealing framework 
                  %for if, else statements below
                  annealer(tsNode,1)=annealer(tsNode,1)+alphai(j,1)*(1-damage(tsNode-1,j)); 
              end
          else
              annealer(tsNode,1)=annealer(tsNode,1)+alphai(tsNode,1);
          end

        annealer(tsNode,1)=annealer(tsNode,1)/1e16; %set to a more palatable number
        
        %this is the threshold check, have we exceeded the old
        %amount of damage retained in memoryDam variable
        %comment out the first line after the if statement to turn this
        %feature off
        if(annealer(tsNode,1)<memoryDam)
            annealer(tsNode,1)=memoryDam; 
        end
         
        %inner for loop goes from young -> old. Note that order does not
        %matter for equivalent time calculations (see Duddy et al., 1988),
        %even though the referenced tT spot here is different than the
        %Guenthner et al 2013 annealing code.
        for node=tsNode:-1:1
           
           %now step through each possible annealing model using annealer
           %each case needs to have its own teq that can be added to the t
           %variable listed below. there will be slight mistmatch of teq 
           %when damage crosses a threshold to a new annealing model 
           %(but is using the previous steps annealing model teq) but
           %the error should be slight
           
           %intuitively, the inner selection statements should be looking
           %up the damage matix row above (tsNode-1) not the column to the
           %left (node+1), but this is a) not possible when tsNode=1, and
           %b) not really necessary because these cells are equivalent. But
           %you still need to consider the case when node==tsNode as
           %special because there will be a zero at the node+1 cell in that
           %case. Hence, the or statements for the first nested ifs.
           
           t=(tT(node,1)*1000000*365.25*24*60*60)-(tT(node+1,1)*1000000*365.25*24*60*60)+teq; %in seconds
           
           if(annealer(tsNode,1)<=46) 
               %model parameters from table 4 of Ginster et al. 2019
               if(damage(tsNode,node+1)>0.6 || node==tsNode)
                C0=0.3156;
                C1=1.06e-5;
                C2=-30;
                C3=0;
                alpha=0.3319;
                
                %in certain scenarios for some sets of parameters, the
                %inside part of the brackets (C0+...) can be negative,
                %which in turn leads to imaginary values as in (1/alpha) of
                %a negative. This translates into essentially r values that
                %are less than 0, which can therefore be treated as a lower
                %bound. the next if else statement handles this problem and
                %is repeated for all of the instances below (even where
                %it's not totally necessary) for consistency purposes
                
                r=(C0+C1*((log(t)-C2)/((1/temp)-C3))); 
                if(r<0)
                    r=0; 
                else
                    r=r^(1/alpha);
                end
                damage(tsNode,node)=r;
                if(node==1) %so that it doesn't try to calc. teq for non-existent step
                    break
                end
                
                if(damage(tsNode,node)>0.000001)
                    temp=(tT(node-1,2)+273.15+tT(node,2)+273.15)/2; %temperature at next (older) node
                    teq=exp(C2+((1/temp)-C3)*((r^alpha)-C0)/C1);
                end
               elseif(damage(tsNode-1,node)>0.3) 
                C0=0;
                C1=1.46e-5;
                C2=-11.4;
                C3=3.72e-4;
                alpha=0.6233;

                r=(C0+C1*((log(t)-C2)/((1/temp)-C3))); 
                if(r<0)
                    r=0; 
                else
                    r=r^(1/alpha);
                end
                damage(tsNode,node)=r;
                if(node==1) %so that it doesn't try to calc. teq for non-existent step
                    break
                end
                
                if(damage(tsNode,node)>0.000001)
                    temp=(tT(node-1,2)+273.15+tT(node,2)+273.15)/2; %temperature at next (older) node in row
                    teq=exp(C2+((1/temp)-C3)*((r^alpha)-C0)/C1);
                end
               else
                C0=-0.2886;
                C1=1.32e-5;
                C2=-20.2;
                C3=0.0001;
                alpha=1.6396;

                r=(C0+C1*((log(t)-C2)/((1/temp)-C3))); 
                if(r<0)
                    r=0; 
                else
                    r=r^(1/alpha);
                end
                damage(tsNode,node)=r;
                if(node==1) %so that it doesn't try to calc. teq for non-existent step
                    break
                end
                
                if(damage(tsNode,node)>0.000001)
                    temp=(tT(node-1,2)+273.15+tT(node,2)+273.15)/2; %temperature at next (older) node 
                    teq=exp(C2+((1/temp)-C3)*((r^alpha)-C0)/C1);
                end
               end
           elseif(annealer(tsNode,1)<=95)
               %model parameters from table 5 of Ginster et al. 2019
               if(damage(tsNode-1,node)>0.4 || node==tsNode) 
                C0=0.1675;
                C1=1.27e-5;
                C2=-29.8;
                C3=0;
                alpha=0.4539;

                r=(C0+C1*((log(t)-C2)/((1/temp)-C3))); 
                if(r<0)
                    r=0; 
                else
                    r=r^(1/alpha);
                end
                damage(tsNode,node)=r;
                if(node==1) %so that it doesn't try to calc. teq for non-existent step
                    break
                end
                
                if(damage(tsNode,node)>0.000001)
                    temp=(tT(node-1,2)+273.15+tT(node,2)+273.15)/2; %temperature at next (older) node 
                    teq=exp(C2+((1/temp)-C3)*((r^alpha)-C0)/C1);
                end
               elseif(damage(tsNode-1,node)>0.3) 
                C0=0;
                C1=1.64e-5;
                C2=-23.0;
                C3=1.75e-4;
                alpha=0.4922;

                r=(C0+C1*((log(t)-C2)/((1/temp)-C3))); 
                if(r<0)
                    r=0; 
                else
                    r=r^(1/alpha);
                end
                damage(tsNode,node)=r;
                if(node==1) %so that it doesn't try to calc. teq for non-existent step
                    break
                end
                
                if(damage(tsNode,node)>0.000001)
                    temp=(tT(node-1,2)+273.15+tT(node,2)+273.15)/2; %temperature at next (older) node 
                    teq=exp(C2+((1/temp)-C3)*((r^alpha)-C0)/C1);
                end
               else
                C0=-0.1604;
                C1=6.42e-6;
                C2=-10.0;
                C3=5.6e-4;
                alpha=2.4010;
             
                r=(C0+C1*((log(t)-C2)/((1/temp)-C3))); 
                if(r<0)
                    r=0; 
                else
                    r=r^(1/alpha);
                end
                damage(tsNode,node)=r;
                if(node==1) %so that it doesn't try to calc. teq for non-existent step
                    break
                end
                
                if(damage(tsNode,node)>0.000001)
                    temp=(tT(node-1,2)+273.15+tT(node,2)+273.15)/2; %temperature at next (older) node 
                    teq=exp(C2+((1/temp)-C3)*((r^alpha)-C0)/C1);
                end
                   
               end
           
           elseif(annealer(tsNode,1)<=220)
               %model parameters from Table 6A of Ginster et al. 2019
               if(damage(tsNode-1,node)>0.5 || node==tsNode) 
                C0=-0.2153;
                C1=1.75e-5;
                C2=-30.2;
                C3=1e-9;
                alpha=0.9165;

                r=(C0+C1*((log(t)-C2)/((1/temp)-C3))); 
                if(r<0)
                    r=0; 
                else
                    r=r^(1/alpha);
                end
                damage(tsNode,node)=r;
                if(node==1) %so that it doesn't try to calc. teq for non-existent step
                    break
                end
                
                if(damage(tsNode,node)>0.000001)
                    temp=(tT(node-1,2)+273.15+tT(node,2)+273.15)/2; %temperature at next (older) node 
                    teq=exp(C2+((1/temp)-C3)*((r^alpha)-C0)/C1);
                end
               elseif(damage(tsNode-1,node)>0.4) 
                C0=0;
                C1=1.84e-5;
                C2=-66.4;
                C3=-8.46e-4;
                alpha=0.3500;

                r=(C0+C1*((log(t)-C2)/((1/temp)-C3))); 
                if(r<0)
                    r=0; 
                else
                    r=r^(1/alpha);
                end
                damage(tsNode,node)=r;
                if(node==1) %so that it doesn't try to calc. teq for non-existent step
                    break
                end
                
                if(damage(tsNode,node)>0.000001)
                    temp=(tT(node-1,2)+273.15+tT(node,2)+273.15)/2; %temperature at next (older) node 
                    teq=exp(C2+((1/temp)-C3)*((r^alpha)-C0)/C1);
                end
               else
                C0=-0.0963;
                C1=6.66e-6;
                C2=0;
                C3=8.39e-4;
                alpha=1.9597;
                
                r=(C0+C1*((log(t)-C2)/((1/temp)-C3))); 
                if(r<0)
                    r=0; 
                else
                    r=r^(1/alpha);
                end
                damage(tsNode,node)=r;
                if(node==1) %so that it doesn't try to calc. teq for non-existent step
                    break
                end
                
                if(damage(tsNode,node)>0.000001)
                    temp=(tT(node-1,2)+273.15+tT(node,2)+273.15)/2; %temperature at next (older) node 
                    teq=exp(C2+((1/temp)-C3)*((r^alpha)-C0)/C1);
                end
               end
           else
               %Table 7, Low-T, high damage model from Ginster et al. 2019
               C0=-0.4376;
               C1=1.72e-5;
               C2=-18;
               C3=3.35e-4;
               alpha=0.89999;
               
               
                r=(C0+C1*((log(t)-C2)/((1/temp)-C3))); 
                if(r<0)
                    r=0; 
                else
                    r=r^(1/alpha);
                end
                damage(tsNode,node)=r;
               if(node==1) %so that it doesn't try to calc. teq for non-existent step
                  break
               end
               
                if(damage(tsNode,node)>0.000001)
                    temp=(tT(node-1,2)+273.15+tT(node,2)+273.15)/2; %temperature at next (older) node 
                    teq=exp(C2+((1/temp)-C3)*((r^alpha)-C0)/C1);
                end
           
           end
            
        end
        
        memoryDam=annealer(tsNode,1);
        
    end

    %create final annealedDam vector that gets passed to diffusion solver
    for tsNode=steps-1:-1:1
        for node=tsNode:-1:1
            annealedDam(tsNode,1)=annealedDam(tsNode,1)+alphai(node,1)*(1-damage(tsNode,node));
        end
    end
end


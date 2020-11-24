function [annealedDam] = Ginster2019(steps, alphai, tT)
    %damage annealing model based on Ginster et al. 2019 kinetics
    
    %requires inputs of tT() matrix and alphai() or "unannealed" damage
    %matrix, returns annealedDam matrix
        
    %create damage array, annealedDam vector that gets passed to
    %diffusion model, and alpha dose production vector
    
    %NOTE: damage matrix convention is 1 = complete annealing, 0 = no
    %annealing as per Ginster et al. derivations. But, need to convert to 1
    %= no annealing, 0 = complete annealing in order to keep consistent
    %with legacy code from ZFT annealing framework. This conversion is done
    %by 1-damage(:,:) for annealer and annealedDam
    damage=zeros(steps,steps);
    teq_mat=zeros(steps,steps); %for debugging purposes only
    annealer=zeros(steps,1); %doesn't really need to be a vector, but useful for error checking
    annealedDam=zeros(steps,1);
    memoryDam=0;
    
    %annealed damage array, arrangment is damage(1,1) = oldest, oldest tT
    %step; damage(steps, steps) = youngest, youngest tT step. Annealing steps
    %forward in time (old -> young) in terms of rows (outer for loop)in
    %order to calculate damage level of previous time step, and forward
    %in time (old -> young) in terms of columns within each row (inner for
    %loop). Matrix is arranged such that the rows describe tT path for 
    %a given damage zone, and the columns give the time of damage creation.
    
    %need to calculate damage for the first oldest tT step to satisfy
    %look-ups in for loops below, use the table 4 kinetics because damage
    %here is the "freshest"
    t=(tT(1,1)*1000000*365.25*24*60*60)-(tT(2,1)*1000000*365.25*24*60*60);
    temp=tT(1,2)+273.15;
    C0=-0.2886;
    C1=1.32e-5;
    C2=-20.2;
    C3=0.0001;
    alpha=1.6396;
    r=(C0+C1*((log(t)-C2)/((1/temp)-C3)));
    r=r^(1/alpha);
    damage(1,1)=r;
    
    for tsNode=2:steps-1 %has to be to steps-1 for tT matrix
        %get damage ingrowth for each step from alphai vector and combine
        %with annealer vector and damage matrix from next oldest time step
        %(i.e.row), number is needed to decide on correct annealing 
        %equation to use

              for j=tsNode-1:-1:1 
                  %sum the damage in a given row of damage matrix from 
                  %NEXT OLDEST time step row to the row you're currently in
                  %, i.e. tsNode-1, convert to 1 = no annealing framework 
                  %for if, else statements below
                  annealer(tsNode-1,1)=annealer(tsNode-1,1)+alphai(j,1)*(1-damage(tsNode-1,j)); 
              end
          
        annealer(tsNode-1,1)=annealer(tsNode-1,1)/1e16; %set to a more palatable number
        
        %this is the threshold check, have we exceeded the old
        %amount of damage retained in memoryDam variable
        if(annealer(tsNode-1,1)<memoryDam)
            annealer(tsNode-1,1)=memoryDam;
        else
            memoryDam=annealer(tsNode-1,1);
        end
        
        %inner for loop goes from old -> young
        for node=1:tsNode
           
           temp=(tT(tsNode,2)+273.15);
           %now step through each possible annealing model using annealer
           %each case needs to have its own teq that can be added to the t
           %variable listed below. there will be slight mistmatch of teq 
           %when damage crosses a threshold to a new annealing model 
           %(but is using the previous steps annealing model teq) but
           %the error should be slight
           
           %the inner selection statements looks up the damage matix row 
           %above (tsNode-1). 
           
           if(annealer(tsNode-1,1)<=46) 
               %model parameters from table 4 of Ginster et al. 2019
               if(damage(tsNode-1,node)>0.6)
                C0=0.3156;
                C1=1.06e-5;
                C2=-30;
                C3=0;
                alpha=0.3319;
                
                %calculate equivalent time using the previous tT step's r 
                %but the current kinetics (Note: this was the big change
                %that Kerry caught that necessitated the new version and
                %the redesign with respect to for loop setup
                
                teq=exp(C2+((1/temp)-C3)*((damage(tsNode-1,node)^alpha)-C0)/C1);
                t=(tT(node,1)*1000000*365.25*24*60*60)-(tT(node+1,1)*1000000*365.25*24*60*60)+teq; %in seconds
                
                %now get the new r for this time step
                r=(C0+C1*((log(t)-C2)/((1/temp)-C3))); 
                
                %in certain scenarios for some sets of parameters, the
                %inside part of the brackets (C0+...) can be negative,
                %which in turn leads to imaginary values as in (1/alpha) of
                %a negative. This translates into essentially r values that
                %are less than 0, which can therefore be treated as a lower
                %bound. the next if else statement handles this problem and
                %is repeated for all of the instances below (even where
                %it's not totally necessary) for consistency purposes
                
                if(r<0)
                    r=0; 
                else
                    r=r^(1/alpha);
                end
                damage(tsNode,node)=r;
                
               elseif(damage(tsNode-1,node)>0.3) 
                C0=0;
                C1=1.46e-5;
                C2=-11.4;
                C3=3.72e-4;
                alpha=0.6233;
                
                %calculate equivalent time using the previous tT step's r 
                %but the current kinetics
                
                teq=exp(C2+((1/temp)-C3)*(((damage(tsNode-1,node)^alpha)-C0)/C1));
                t=(tT(node,1)*1000000*365.25*24*60*60)-(tT(node+1,1)*1000000*365.25*24*60*60)+teq; %in seconds
                
                %now get the new r for this time step
                r=(C0+C1*((log(t)-C2)/((1/temp)-C3))); 
                if(r<0)
                    r=0; 
                else
                    r=r^(1/alpha);
                end
                damage(tsNode,node)=r;
               
               else
                C0=-0.2886;
                C1=1.32e-5;
                C2=-20.2;
                C3=0.0001;
                alpha=1.6396;
                
                %calculate equivalent time using the previous tT step's r 
                %but the current kinetics
                
                teq=exp(C2+((1/temp)-C3)*(((damage(tsNode-1,node)^alpha)-C0)/C1));
                t=(tT(node,1)*1000000*365.25*24*60*60)-(tT(node+1,1)*1000000*365.25*24*60*60)+teq; %in seconds
                
                %now get the new r for this time step

                r=(C0+C1*((log(t)-C2)/((1/temp)-C3))); 
                if(r<0)
                    r=0; 
                else
                    r=r^(1/alpha);
                end
                damage(tsNode,node)=r;
                
               end
           elseif(annealer(tsNode-1,1)<=95)
               %model parameters from table 5 of Ginster et al. 2019
               if(damage(tsNode-1,node)>0.4) 
                C0=0.1675;
                C1=1.27e-5;
                C2=-29.8;
                C3=0;
                alpha=0.4539;
                
                %calculate equivalent time using the previous tT step's r 
                %but the current kinetics
                
                teq=exp(C2+((1/temp)-C3)*(((damage(tsNode-1,node)^alpha)-C0)/C1));
                t=(tT(node,1)*1000000*365.25*24*60*60)-(tT(node+1,1)*1000000*365.25*24*60*60)+teq; %in seconds
                
                %now get the new r for this time step

                r=(C0+C1*((log(t)-C2)/((1/temp)-C3))); 
                if(r<0)
                    r=0; 
                else
                    r=r^(1/alpha);
                end
                damage(tsNode,node)=r;
                
               elseif(damage(tsNode-1,node)>0.3) 
                C0=0;
                C1=1.64e-5;
                C2=-23.0;
                C3=1.75e-4;
                alpha=0.4922;
                
                %calculate equivalent time using the previous tT step's r 
                %but the current kinetics
                
                teq=exp(C2+((1/temp)-C3)*(((damage(tsNode-1,node)^alpha)-C0)/C1));
                t=(tT(node,1)*1000000*365.25*24*60*60)-(tT(node+1,1)*1000000*365.25*24*60*60)+teq; %in seconds
                
                %now get the new r for this time step

                r=(C0+C1*((log(t)-C2)/((1/temp)-C3))); 
                if(r<0)
                    r=0; 
                else
                    r=r^(1/alpha);
                end
                damage(tsNode,node)=r;
                
               else
                C0=-0.1604;
                C1=6.42e-6;
                C2=-10.0;
                C3=5.6e-4;
                alpha=2.4010;
                
                %calculate equivalent time using the previous tT step's r 
                %but the current kinetics
                
                teq=exp(C2+((1/temp)-C3)*(((damage(tsNode-1,node)^alpha)-C0)/C1));
                t=(tT(node,1)*1000000*365.25*24*60*60)-(tT(node+1,1)*1000000*365.25*24*60*60)+teq; %in seconds
                
                %now get the new r for this time step
             
                r=(C0+C1*((log(t)-C2)/((1/temp)-C3))); 
                if(r<0)
                    r=0; 
                else
                    r=r^(1/alpha);
                end
                damage(tsNode,node)=r;
                
               end
           
           elseif(annealer(tsNode-1,1)<=220)
               %model parameters from Table 6A of Ginster et al. 2019
               if(damage(tsNode-1,node)>0.5) 
                C0=-0.2153;
                C1=1.75e-5;
                C2=-30.2;
                C3=1e-9;
                alpha=0.9165;
                
                %calculate equivalent time using the previous tT step's r 
                %but the current kinetics
                
                teq=exp(C2+((1/temp)-C3)*(((damage(tsNode-1,node)^alpha)-C0)/C1));
                t=(tT(node,1)*1000000*365.25*24*60*60)-(tT(node+1,1)*1000000*365.25*24*60*60)+teq; %in seconds
                
                %now get the new r for this time step

                r=(C0+C1*((log(t)-C2)/((1/temp)-C3))); 
                if(r<0)
                    r=0; 
                else
                    r=r^(1/alpha);
                end
                damage(tsNode,node)=r;
                
               elseif(damage(tsNode-1,node)>0.4) 
                C0=0;
                C1=1.84e-5;
                C2=-66.4;
                C3=-8.46e-4;
                alpha=0.3500;
                
                %calculate equivalent time using the previous tT step's r 
                %but the current kinetics
                
                teq=exp(C2+((1/temp)-C3)*(((damage(tsNode-1,node)^alpha)-C0)/C1));
                t=(tT(node,1)*1000000*365.25*24*60*60)-(tT(node+1,1)*1000000*365.25*24*60*60)+teq; %in seconds
                
                %now get the new r for this time step

                r=(C0+C1*((log(t)-C2)/((1/temp)-C3))); 
                if(r<0)
                    r=0; 
                else
                    r=r^(1/alpha);
                end
                damage(tsNode,node)=r;
               
               else
                C0=-0.0963;
                C1=6.66e-6;
                C2=0;
                C3=8.39e-4;
                alpha=1.9597;
                
                %calculate equivalent time using the previous tT step's r 
                %but the current kinetics
                
                teq=exp(C2+((1/temp)-C3)*(((damage(tsNode-1,node)^alpha)-C0)/C1));
                t=(tT(node,1)*1000000*365.25*24*60*60)-(tT(node+1,1)*1000000*365.25*24*60*60)+teq; %in seconds
                
                %now get the new r for this time step
                
                r=(C0+C1*((log(t)-C2)/((1/temp)-C3))); 
                if(r<0)
                    r=0; 
                else
                    r=r^(1/alpha);
                end
                damage(tsNode,node)=r;
                
               end
           else
               %Table 7, Low-T, high damage model from Ginster et al. 2019
               C0=-0.4376;
               C1=1.72e-5;
               C2=-18;
               C3=3.35e-4;
               alpha=0.89999;
               
               %calculate equivalent time using the previous tT step's r 
               %but the current kinetics
                
               teq=exp(C2+((1/temp)-C3)*(((damage(tsNode-1,node)^alpha)-C0)/C1));
               t=(tT(node,1)*1000000*365.25*24*60*60)-(tT(node+1,1)*1000000*365.25*24*60*60)+teq; %in seconds
                
               %now get the new r for this time step
               
               r=(C0+C1*((log(t)-C2)/((1/temp)-C3))); 
               if(r<0)
                   r=0; 
               else
                   r=r^(1/alpha);
               end
               damage(tsNode,node)=r;
               
           end
            teq_mat(tsNode,node)=teq;
        end
        
    end

    %create final annealedDam vector that gets passed to diffusion solver
    for tsNode=1:steps
        for node=1:tsNode
            annealedDam(tsNode,1)=annealedDam(tsNode,1)+alphai(node,1)*(1-damage(tsNode,node));
        end
    end
    
end


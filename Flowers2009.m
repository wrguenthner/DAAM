function [damage] = Flowers2009(steps,tT)

    %damage annealing model based on Flowers et al. 2009 AFT kinetics
    
    %requires inputs of tT() matrix and alphai() or "unannealed" damage
    %matrix, returns annealedDam matrix
    
    %AFT values from Flowers fit of Ketcham data
    C0=0.39528;
    C1=0.01073;
    C2=-65.12969;
    C3=-7.91715;
    alpha=0.04672;
    rmr0=0.83; %value used in HeFTy, but subject to change
    kappa=1.04-rmr0;

      
    %now create damage array to hold tracks
    damage=zeros(steps,steps);
      
      for tsNode=steps-1:-1:1
          teq=0;
          temp=(tT(tsNode,2)+273.15+tT(tsNode+1,2)+273.15)/2;
          for node=tsNode:-1:1
              t=(tT(node,1)*1000000*365.25*24*60*60)-(tT(node+1,1)*1000000*365.25*24*60*60)+teq;
              r=((C0+C1*((log(t)-C2)/(log(1/temp)-C3)))^(1/alpha)+1)^-1;
              damage(tsNode,node)=r;
              if(node==1)
                  break
              end
              temp=(tT(node-1,2)+273.15+tT(node,2)+273.15)/2;
              teq=exp(C2+(log(1/temp)-C3)*(((1/r)-1)^alpha-C0)/C1);
          end
      end
    
    %rmr0 and density conversion (rhor)
    for tsNode=steps-1:-1:1
        for node=tsNode:-1:1
            if(damage(tsNode,node)>=rmr0)
                damage(tsNode,node)=((damage(tsNode,node)-rmr0)/(1-rmr0)).^kappa;
            else
                damage(tsNode,node)=0;
            end
            
            if(damage(tsNode,node)>=0.765)
                damage(tsNode,node)=1.6*damage(tsNode,node)-0.6;
            else
                damage(tsNode,node)=9.205*damage(tsNode,node)*damage(tsNode,node)-9.157*damage(tsNode,node)+2.269;
            end
        end
    end
                  
end


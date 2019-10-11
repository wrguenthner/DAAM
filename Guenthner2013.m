function [damage] = Guenthner2013(steps,tT)

    %damage annealing model based on Guenthner et al. 2013 FCT kinetics
    
    %requires inputs of tT() matrix and alphai() or "unannealed" damage
    %matrix, returns annealedDam matrix
    
    %ZFT values from Guenthner fit of Yamada data
    C0=6.24534;
    C1=-0.11977;
    C2=-314.937;
    C3=-14.2868;
    alpha=-0.05721;
    
    equivTotAnnealLen=0.36/1.25+.2;
      
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
      
      
    %volume conversion
    damage(:,:)=1.25*(damage(:,:)-.2);
    
    %zero out the stuff less than equivalent total annealing length
    damage(damage<equivTotAnnealLen)=0;
                  
end


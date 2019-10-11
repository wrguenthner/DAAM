function [heModelAge] = ZircDif(radius, nodes, diffusivities, tT, U235atom, U238atom, Thatom)
%Alpha ejection and diffusion model for zircon (specific wrt alpha stopping
%distances, but otherwise universal). Requires grain radius, number of
%nodes for CN algorithm, diffusivities, tT history, and atoms of U and Th
%isotopes. Calculates a He age from diffusion profile and returns. 
    
    %decay constants in 1/yr
    lambda235=.00000000098485;
    lambda238=.000000000155125;
    lambda232=.000000000049475;

    %alpha ejection correction factors
    aEj238=zeros(nodes,1);
    aEj235=zeros(nodes,1);
    aEj232=zeros(nodes,1);
    
    %alpha stopping distances in microns, from Hourigan et al., 2005
    S238=16.97;
    S235=19.64;
    S232=19.32;
    
    rStep=radius/nodes;
    
    %fill in alpha ejection array for each isotope (no zonation)
    %in atoms/g
    for j=1:nodes
        X0=(j-.5)*rStep;
        if(X0>=radius-S238)
            Xstar=(X0^2+radius^2-S238^2)/(2*X0);
            Ft=.5+(Xstar-X0)/(2*S238);
            aEj238(j,1)=Ft*U238atom;
        else
            aEj238(j,1)=U238atom;
        end
        
        if(X0>=radius-S235)
            Xstar=(X0^2+radius^2-S235^2)/(2*X0);
            Ft=.5+(Xstar-X0)/(2*S235);
            aEj235(j,1)=Ft*U235atom;
        else
            aEj235(j,1)=U235atom;
        end
        
        if(X0>=radius-S232)
            Xstar=(X0^2+radius^2-S232^2)/(2*X0);
            Ft=.5+(Xstar-X0)/(2*S232);
            aEj232(j,1)=Ft*Thatom;
        else
            aEj232(j,1)=Thatom;
        end
    end
    
    %calculate the total amount of isotope in each Xtal (no zonation)
    %in atoms/g
    innerVol=0;
    position=0;
    total238=0;
    total235=0;
    total232=0;
    for j=1:nodes
        position=position+rStep;
        outerVol=position^3;
        total238=total238+U238atom*(outerVol-innerVol);
        total235=total235+U235atom*(outerVol-innerVol);
        total232=total232+Thatom*(outerVol-innerVol);
        innerVol=outerVol;
    end
    
    %scale by He production and volume with a base of /(1.33333*PI)
    total238=8*total238;
    total235=7*total235;
    total232=6*total232;
    
    
    %alpha ejection correction factors
    innerVol=0;
    rad=0;
    d232=0;
    d235=0;
    d238=0;
    nd232=0;
    nd235=0;
    nd238=0;
    for j=1:nodes
        rad=rad+rStep;
        outerVol=rad*rad*rad;
        vol=outerVol-innerVol;
        d232=d232+vol*aEj232(j,1);
        d235=d235+vol*aEj235(j,1);
        d238=d238+vol*aEj238(j,1);
        nd232=nd232+vol*Thatom;
        nd235=nd235+vol*U235atom;
        nd238=nd238+vol*U238atom;
        innerVol=outerVol;
    end
    
    if(Thatom==0)
        ft232=0;
    else
        ft232=d232/nd232;
    end
    ft235=d235/nd235;
    ft238=d238/nd238;
   
    %diffusion model
    oldSphere=zeros(nodes,1);
    sphere=zeros(nodes,1);
    a=zeros(nodes,1);
    b=zeros(nodes,1);
    c=zeros(nodes,1);
    d=zeros(nodes,1);

    oldSphere(:,1)=100;
    for j=1:size(tT,1)-1
        t1=tT(j,1)*1000000; %years
		t2=tT(j+1,1)*1000000; %years
        timeStep=(t1-t2)*365.25*24*60*60;%seconds
        beta=(2*rStep^2)/(diffusivities(j,1)*timeStep);
        a(:,1)=1;
        b(:,1)=-2-beta;
        c(:,1)=1;
        sphere(:,1)=0; %reset sphere
        
        %production term for first node
        %A is in atoms/g
        A=8*aEj238(1,1)*(exp(lambda238*t1)-exp(lambda238*t2))+7*aEj235(1,1)*(exp(...
            lambda235*t1)-exp(lambda235*t2))+6*aEj232(1,1)*(exp(lambda232*t1)-exp(lambda232*t2));      
        %neumann boundary condition
        b(1,1)=-3-beta;
        d(1,1)=(3-beta)*oldSphere(1,1)-oldSphere(2,1)-A*(1-.5)*rStep*beta; 
        %zero Dirichlet boundary condition
        A=8*aEj238(nodes,1)*(exp(lambda238*t1)-exp(lambda238*t2))+7*aEj235(nodes,1)*(exp(...
            lambda235*t1)-exp(lambda235*t2))+6*aEj232(nodes,1)*(exp(lambda232*t1)-exp(lambda232*t2));
        d(nodes,1)=-oldSphere(nodes-1,1)+(2-beta)*oldSphere(nodes,1)-A*...
           (nodes-.5)*rStep*beta;
        for k=2:nodes-1
            A=8*aEj238(k,1)*(exp(lambda238*t1)-exp(lambda238*t2))+7*aEj235(k,1)*(exp(...
               lambda235*t1)-exp(lambda235*t2))+6*aEj232(k,1)*(exp(lambda232*t1)-exp(lambda232*t2));
            d(k,1)=-oldSphere(k-1,1)+(2-beta)*oldSphere(k,1)-...
                oldSphere(k+1,1)-A*(k-.5)*rStep*beta;
        end
        
        %tridiagonal matrix algorithm
        c(1,1)=c(1,1)/b(1,1);
        d(1,1)=d(1,1)/b(1,1);
        for k=2:nodes-1
            c(k,1)=c(k,1)/(b(k,1)-c(k-1,1)*a(k,1));
            d(k,1)=(d(k,1)-d(k-1,1)*a(k,1))/(b(k,1)-c(k-1,1)*a(k,1));
        end
        d(nodes,1)=(d(nodes,1)-d(nodes-1,1))/(b(nodes-1,1)-c(nodes-1,1));
                
        %back substitute
        sphere(nodes,1)=d(nodes,1);
        for k=nodes-1:-1:1
            sphere(k,1)=d(k,1)-c(k,1)*sphere(k+1,1);
        end
        oldSphere(:,1)=sphere(:,1); %oldSphere becomes sphere and repeat
    end
   
    %convert sphere u array to He
    minHe=sphere(1,1)/(.5*rStep); %zero He check setup
    heProfile=zeros(nodes,1);
    integrand=zeros(nodes,1);
    for j=1:nodes
        position=(j-.5)*rStep;
        heProfile(j,1)=sphere(j,1)/position;
        integrand(j,1)=heProfile(j,1)*4*pi*position^2; %integration function
    end
    
    %Romberg integration algorithms
    %from NR
    decdigs=10;
    rom=zeros(2,decdigs);
    b=radius;
    a=0;
    romall=integrand;
    h=b-a;
    rom(1,1)=h*(romall(1)+romall(end))/2;
    for j=2:decdigs
        st=2^(decdigs-j+1);
        rom(2,1)=(rom(1,1)+h*sum(romall((st/2)+1:st:2^(decdigs-1))))/2;
        for k=1:j-1
           rom(2,k+1)=((4^k)*rom(2,k)-rom(1,k))/((4^k)-1);
        end
        rom(1,1:j)=rom(2,1:j);
        h=h/2;
    end
    res=rom(1,decdigs);
    %fill in crystal center (half-node before first node)
    totalHe=res; %in atoms/g
    totalHe=totalHe+.5*rStep*(integrand(1,1)*55/24-integrand(2,1)*59/24+...
        integrand(3,1)*37/24-integrand(4,1)*9/24);
    %convert to same basis as the other isotope totals
    totalHe=totalHe/(4*pi/3);
    %zero He check
    if minHe<-heProfile(1,1)*.05
        totalHe=0;
    end
    
    %iterative age calculation
    leftSum=total238*ft238+total235*ft235+total232*ft232+totalHe;
    loAge=0;
    hiAge=tT(1,1)*1000000; %years
    ageConv=100; %get it to within a hundred years
    while hiAge-loAge>ageConv
        midAge=(hiAge+loAge)/2;
        midVal=total238*ft238*exp(lambda238*midAge)+total235*ft235*exp(...
            lambda235*midAge)+total232*ft232*exp(lambda232*midAge);
        if midVal<leftSum
            loAge=midAge;
        else
            hiAge=midAge;
        end
    end
    heModelAge=(hiAge+loAge)/2;
    heModelAge=heModelAge/1000000; %convert to Ma

end


%% LINMA2720 - Project : Synchronization of biological pulse-coupled oscillators 
%===================================================================================================================
% Minh-Phuong Tran
%
% PARAMETERS : 
% N = # of oscillators
% T = Period
% dt = Increments in time
% t = Time vector
% binit = first value of parameter b
% epsiloninit = first value of parameter epsilon
% K = # of different values of epsilon (and b) that will be considered
% M = # of times we compute/ sample -> to make an average
%
% calcmode = Mode of calculation ["time","iterate"]
%           "time" : calculate with functions f and g, with time t and phi
%           "iterate" : calculate with functions h and R, only when 1
%                        oscillator fires
%
% displacemode = Displacement of the oscillators ["staticline",
%                "staticsquare", "brownian"]
%           "staticline" : no motion, in line
%           "staticsquare" : no motion, in a square grid /!\ N must be a
%                            square nb /!\
%           "brownian" : brownian motion
%
% animode = Mode of animation [0 , 1]
%           0 : no animation
%           1 : animation
%----------------------------------------------------------------------------
% READ ME :  
%
% To use this script, run 1 set of parameters in the PARAMETERS section.
% Then, run the MAIN section.
% Finally, run any plot or animation (in the GRAPHS and ANIMATIONS section) of interest (corresponding to
% the selected parameters).





%% PARAMETERS
%====================================================================================================================
%% Parameters : 2 OSCILLATORS graph of time of sync
%--------------------------------------------------------------------------------------------------------------------

clc
clear all

N = 2; 
T = 5; 
dt = 0.01; 
t = 0:dt:200;
tfinanim = 0;

binit = 0.05; 
epsiloninit = 0.05;

K = floor(1/epsiloninit);
M = 100;

displacemode = "staticline";
calcmode = "iterate"; % "time" or "iterate" -> different results
animode = 0; 

%% Parameters : 2 OSCILLATORS animation
%---------------------------------------------------------------------------------------------------------------------

clc
clear all

N = 2; 
T = 1; 
dt = 0.01; 
t = 0:dt:5;
tfinanim = 5;

binit = 2; 
epsiloninit = rand()/2; 

K = 1;
M = 1;

displacemode = "staticline";
calcmode = "time"; 
animode = 1; 

%% Parameters : N OSCILLATORS animation
%---------------------------------------------------------------------------------------------------------------------

clc
clear all

N = 25; 
T = 3; 
dt = 0.05; 
t = 0:dt:10;
tfinanim = 6;

binit = 2; 
epsiloninit = rand([1,N])/10; 

K = 1;
M = 1;

displacemode = "staticsquare";
calcmode = "time"; 
animode = 1; 

%% Parameters : N OSCILLATORS column graph of firing times
%---------------------------------------------------------------------------------------------------------------------

clc
clear all

N = 100; 
T = 1; 
dt = 0.01; 
t = 0:dt:10;
tfinanim = 10;

binit = 2; 
epsiloninit = rand([1,N])/150; 

K = 1;
M = 1;

displacemode = "staticsquare";
calcmode = "time"; 
animode = 1; 

%% MAIN 
%===================================================================================================================
% Calculation
%---------------------------------------------------------------------------------------------------------------------

Tsync = zeros(1,K);
Tsyncall = zeros(1,M);
for k = 1:K
        epsilon = epsiloninit*k;
        b = binit*k;
        delta = 0;
        if(N== 2)
            delta = 1 - g1(1 - epsilon,b); %% Pour eviter la synchro des le premier tir
        end
        for m = 1:M
            phi0 = delta + (1-delta).*rand(1,N-1);
            [Tsyncall(m), y, pos] = tsync(phi0,b,epsilon,t,dt,T,N,calcmode,animode,displacemode);
        end
        condiTsync = (Tsyncall>0);
        Tsync(k) = sum(Tsyncall.*condiTsync)/max(sum(condiTsync),1);
end

%% GRAPHS and ANIMATIONS
%===================================================================================================================
%% Plot of f1 and g1
%--------------------------------------------------------------------------------------------------------------------

figure
plot(t,f1(b,t)); hold on; plot(t,g1(f1(b,t),b));

%% Plot of time of sync
%--------------------------------------------------------------------------------------------------------------------
B = binit.*(1:1:K);
Eps = epsiloninit.*(1:1:K);
Prodbeps = B.*Eps;
figure
plot(Prodbeps, Tsync); hold on; plot(Prodbeps, 1./Prodbeps);
xlabel("Produit (\epsilon b)");
title("Nombre d'iterations necessaires pour arriver a la synchronisation");
legend("Nombre d'iterations","1/(\epsilon b)");

%% Column graph of firing times
%--------------------------------------------------------------------------------------------------------------------
lights = not(y);
sumlights = sum(lights,2);
gap = floor(1/(dt*20));
nt = floor(length(t)/gap);
t_new = t(1:gap:gap*nt);
sumlightsbis = reshape(sumlights(1:gap*nt), [gap,nt])';
sumlights_new = sum(sumlightsbis,2);

figure
bar(t_new,sumlights_new);

xlabel("temps [période]");
ylabel("nombre d'oscillateurs qui tirent");


%% Animation
%---------------------------------------------------------------------------------------------------------------------

animation(y,t, N, pos, dt, tfinanim);

%% USEFUL FUNCTIONS
%====================================================================================================================
function animation(y,t, N, pos, dt, tfin)
M = max(max(max(abs(pos))))+1;
nfin = ceil(tfin/dt);
lights = not(y);
figure
for i = 1:nfin
    for k = 1:N
    if(lights(i,k) == 1)
        plot(pos(1,k,i),pos(2,k,i),'o','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','green');
        axis([-M,M,-M,M]); 
        hold on;
    else
        plot(pos(1,k,i),pos(2,k,i),'o','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','white');
        axis([-M,M,-M,M]); 
        hold on;
    end
    end
    title(sprintf('time: %3.2f s',t(i)));
    pause(.0001);
    hold off;
end

end


function pos = displacement(displacemode,timet,prev,dt)
N = size(prev,2);
if displacemode == "staticsquare"
    if timet== 0
        Nsq = sqrt(N);
        [posx,posy] = meshgrid(linspace(-Nsq,Nsq,Nsq));
        pos = [posx(:)';posy(:)'];
    else
        pos = prev;
    end
elseif displacemode=="staticline"
    if timet==0
        x = linspace(-N,N,N);
        pos = [x;zeros(1,N)];
    else
        pos = prev;
    end
elseif displacemode == "brownian"
    if timet== 0
        pos = normrnd(0,3,size(prev));
    else
        pos = prev + dt.*normrnd(0,1,size(prev));
    end
end

end

function y = f1(b, phi2)
phi = round(phi2,4);
y = log(1+phi.*(exp(b) - 1))./b;
y = round(y,4);
end
function phi = g1(y2, b)
y = round(y2,4);

phi = (exp(b.*y)-1)./(exp(b)-1);
phi=round(phi,4);
end

function y2 = fire(phi, b, epsilon)
y = f1(b,phi);
condy = (y<1);
y2 = (y + sum(epsilon.*(1-condy))).* condy;
y2 = y2.*(y2<1);

end
function h = firemap(phi, epsilon,b)
lambda = exp(b*epsilon);
phieq = (exp(b*(1+epsilon))-1)/((exp(b)-1)*(lambda+1));
h = -lambda*(phi-phieq)+phieq;
end
function R = returnmap(phi,epsilon,b)
h1 = firemap(phi,epsilon,b);
R = firemap(h1,epsilon,b);
end
function phi = fphi(t,T)
phi = mod(t./T,1);
end
function [Tsync, y, pos] = tsync(phi0,b,epsilon,t,dt,T,N,calcmode,animode,displacemode)
    
    Tsync = 0;
    y = 0;
    pos = 0;
    lengtht = length(t);
    if(calcmode == "time" || animode == 1)
        phi = zeros(lengtht,N);
        phi(1,:) = [0, phi0];
        
        y = zeros(lengtht,N);
        y(1,:) = f1(b,phi(1,:));
        
        pos = zeros(2,N,lengtht);
        pos(:,:,1) = displacement(displacemode,0,zeros(2,N),dt);

    end
    R = phi0;
    for i = 2 : lengtht

        if(calcmode == "time" || animode == 1)
            pos(:,:,i) = displacement(displacemode,t(i),pos(:,:,i-1),dt);
            phi(i,:) = min(g1(y(i-1,:),b)+dt./T,1);
            y(i,:) = fire(phi(i,:), b, epsilon);            
            if(sum(y(i,:))==0 && Tsync ==0)
                Tsync= t(i);
                if(animode == 0)
                    break
                end
            end
        end
        if(calcmode=="iterate")
            R = returnmap(R,epsilon,b);
            if(R>1|| R<0)
                Tsync = i;
                %disp("out of bound");
                break
            end
        end
    end
      
   
end
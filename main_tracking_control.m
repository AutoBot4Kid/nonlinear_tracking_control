clear all ; clc ; close all
data = load('energyDataSave.mat');
tvec = data.ocpsol(:,1);
xref = data.ocpsol(:,2); yref = data.ocpsol(:,3);
xdref = data.ocpsol(:,4); ydref = data.ocpsol(:,5);
xddref = data.ocpsol(:,6); yddref = data.ocpsol(:,7);
thref = data.ocpsol(:,8); vref = data.ocpsol(:,9);
wref = data.ocpsol(:,10); gammaref =data.ocpsol(:,11);



m= [xref(1)-0.02,yref(1)-0.02,thref(1),vref(1)];
Xgain = [1,1,1,1,1,1];
Ygain = [1,1,1,1,1,1];
sig = sqrt(0.05);
dt = (tvec(2)-tvec(1));
t = tvec;

for k = 1:length(t)-1
    x_current = m(k,1);
    y_current = m(k,2);
    th_current = m(k,3);
    v_current = m(k,4);


    x_ref = interp1(tvec,xref,t(k));
    x_refd = interp1(tvec,xdref,t(k));
    x_refdd = interp1(tvec,xddref,t(k));

    y_ref =interp1(tvec,yref,t(k));
    y_refd =interp1(tvec,ydref,t(k));
    y_refdd =interp1(tvec,yddref,t(k));

    v_refd = interp1(tvec,gammaref,t(k));
    v_ref = interp1(tvec,vref,t(k));
    w_ref  = interp1(tvec,wref,t(k));

    ex = x_ref-x_current;
    ey= y_ref-y_current;
    exdot = x_refd-v_current*cos(th_current);
    eydot = y_refd-v_current*sin(th_current);

    k1 =  Xgain(1);
    k2 =  Xgain(2);
    c1 =  Xgain(3);
    c2 =  Xgain(4);
    eps1 = Xgain(5);
    eps2 = Xgain(6);

    k3 = Ygain(1);
    k4 = Ygain(2);
    c3 = Ygain(3);
    c4 = Ygain(4);
    eps3 = Ygain(5);
    eps4 = Ygain(6);


    zeta1 = x_refdd + (k1 +(k2/(ex^2+eps1)^(1/2)))*ex + (c1 +(c2/(ex^2+eps2)^(1/2)))*exdot;
    zeta2 = y_refdd + (k3 +(k4/(ey^2+eps3)^(1/2)))*ey + (c3 +(c4/(ey^2+eps4)^(1/2)))*eydot;

    lam = exp(-1*(v_current^2/sig^2));

    Q = lam*eye(2) ;
    p =eye(2);
    Amat = [cos(th_current),-v_current*sin(th_current);
        sin(th_current),v_current*cos(th_current)];
    cont(k,:) = pinv(Amat'*p*Amat+Q)*(Amat'*p*[zeta1;zeta2]+Q*[v_refd;w_ref]);

    m(k+1,:) = m(k,:)+[v_current*cos(th_current)*dt,v_current*sin(th_current)*dt,cont(k,2)*dt, cont(k,1)*dt;]+ randn(1,4)*0;
end

figure(1)
subplot(2,2,1)
plot(m(:,1),m(:,2),'LineWidth',2)
hold on
plot(xref,yref,'--','LineWidth',2)
grid on
subplot(2,2,2)
plot(t(:),m(:,3),'LineWidth',2)
hold on
plot(tvec,thref,'--','LineWidth',2)
grid on
subplot(2,2,3)
plot(t(:),m(:,4),'LineWidth',2)
hold on
plot(tvec,vref,'--','LineWidth',2)
grid on
subplot(2,2,4)
plot(t(1:end-1),cont(:,2),'LineWidth',2)
hold on
plot(tvec,wref,'--','LineWidth',2)
grid on


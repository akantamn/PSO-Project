set N /1*5/;
alias (N,M);

table Y(M,N)
1 2 3 4 5
1 19.7642 15.8114 3.9528 0 0
2 15.8114 34.2580 5.2705 5.2705 7.9057
3 3.9528 5.2705 40.8461 31.6228 0
4 0 5.2705 31.6228 40.8461 3.9528
5 0 7.9057 0 3.9528 11.8585;

table theta(N,M)
1 2 3 4 5
1 -1.2490 1.8925 1.8925 0 0
2 1.8925 -1.2490 1.8925 1.8925 1.8925
3 1.8925 1.8925 -1.2490 1.8925 0
4 0 1.8925 1.8925 -1.2490 1.8925
5 0 1.8925 0 1.8925 -1.2490;

variables voltage(N);
voltage.l(N)=1;
voltage.lo(N)=0.95;
voltage.up(N)=1.05;
voltage.fx('1')=1.05;
voltage.fx('2')=1.05;
voltage.fx('4')=1.05;


variables delta(N);
delta.l(N)=0;
delta.lo(N)=-pi*20/180;
delta.up(N)=pi*20/180;
delta.fx('1')=0;


variables Pgen(N);
Pgen.fx('3') = 0;
Pgen.fx('5') = 0;

variables Qgen(N);
Qgen.fx('3') = 0;
Qgen.fx('5') = 0;


variables Pload(N);
Pload.fx('1') = 1.25;
Pload.fx('2') = 1.5;
Pload.fx('3') = 1.75;
Pload.fx('4') = 1.80;
Pload.fx('5') = 1.45;


variables Qload(N);
Qload.fx('1') = 0.45;
Qload.fx('2') = 0.5;
Qload.fx('3') = 0.65;
Qload.fx('4') = 0.6;
Qload.fx('5') = 0.55;


variables cost;


equations Costf;
equations Pbal(N);
equations Qbal(N);



$ontext

Pgen.up('1')=4;
Pgen.lo('1')=.50;

Pgen.fx('2')=3.25;
Pgen.fx('4')=3.25;



*Qgen.up('1')=1;
*Qgen.up('2')=1;
*Qgen.up('4')=1;
*Qgen.lo('1')=-0.50;
*Qgen.lo('2')=-0.50;
*Qgen.lo('4')=-0.50;
*$offtext
variables pbr(N,M);
equation pbreq(N,M);
*******************************************************************************

Costf.. cost=e=1;
Pbal(N).. Pgen(N)-Pload(N)=e=sum(M,(voltage(N)*voltage(M)*Y(N,M)*cos(theta(N,M)+delta(M)-delta(N))));
Qbal(N).. Qgen(N)-Qload(N)=e=-sum(M,(voltage(N)*voltage(M)*Y(N,M)*sin(theta(N,M)+delta(M)-delta(N))));
pbreq(N,M).. pbr(N,M)=e=Voltage(N)*Voltage(M)*Y(N,M)*cos(theta(N,M)+delta(M)-delta(N))-Voltage(N)*Voltage(N)*Y(N,M)*cos(theta(N,M));



model ELD /all/;

solve ELD using NLP minimizing cost;

option decimals = 5;


display Pgen.l;
display Qgen.l;
display delta.l;
display voltage.l;
display Pload.l;
display Qload.l;
display pbr.l;
close all;
clear;

% �����N1�̉�]������d�S�܂ł̒��� (m)
lg1 = 1.2;
% �����N1�̉�]�������[�܂ł̒��� (m)
l1 = 2.0;
% �����N2�̉�]������d�S�܂ł̒��� (m)
lg2 = 0.5;
% �����N2�̉�]�������[�܂ł̒��� (m)
l2 = 1.0;

% �{�f�B0�̎��� (kg)
m0 = 5.0;
% �����N1�̎��� (kg)
m1 = 1.0;
% �����N2�̎��� (kg)
m2 = 1.0;
% �����N1�̊������[�����g (kg.m^2)
J1 = 1.0;
% �����N2�̊������[�����g (kg.m^2)
J2 = 1.0;

% �d�͉����x (N/m^2)
g = 9.8;

%% ��ԕ�����
% 'diff(x_0,t,2)=-(49*theta_2+54929*theta_1-3525*F)/19070
% 'diff(theta_1,t,2)=-(1323*theta_2-12005*theta_1+175*F)/1907
% 'diff(theta_2,t,2)=(49294*theta_2-59682*theta_1+870*F)/9535
% x = [ x0, th1, th2, dx0, dth1, dth2 ]';

% ����W��
den = (((lg1^2*lg2^2*m0+J1*lg2^2+J2*lg1^2-2*J2*l1*lg1+J2*l1^2)*m1+(J1*lg2^2+J2*l1^2)*m0+J1*J2)*m2+(J2*lg1^2*m0+J1*J2)*m1+J1*J2*m0);

% ddx�ɂ�����x0, th1, th2, F�̌W��
ddx0_x0 = 0;
ddx0_th1 = ((-g*lg1^2*lg2^2*m1-J1*g*lg2^2-J2*g*l1^2)*m2^2+(-g*lg1^2*lg2^2*m1^2-2*J2*g*l1*lg1*m1)*m2-J2*g*lg1^2*m1^2);
ddx0_th2 = ((g*l1*lg1-g*lg1^2)*lg2^2*m1-J1*g*lg2^2)*m2^2;
ddx0_F = ((lg1^2*lg2^2*m1+J1*lg2^2+J2*l1^2)*m2+J2*lg1^2*m1+J1*J2);

% ddth1�ɂ�����x0, th1, th2, F�̌W��
ddth1_x0 = 0;
ddth1_th1 = -((-g*lg1*lg2^2*m1-J2*g*l1)*m2^2+(-g*lg1*lg2^2*m1^2+(-g*lg1*lg2^2*m0-J2*g*lg1-J2*g*l1)*m1-J2*g*l1*m0)*m2-J2*g*lg1*m1^2-J2*g*lg1*m0*m1);
ddth1_th2 = -((g*l1-g*lg1)*lg2^2*m1+g*l1*lg2^2*m0)*m2^2;
ddth1_F = -((lg1*lg2^2*m1+J2*l1)*m2+J2*lg1*m1);

% ddth2�ɂ�����x0, th1, th2, F�̌W��
ddth2_x0 = 0;
ddth2_th1 = ((((g*lg1^2-g*l1*lg1)*lg2-g*lg1*lg2^2)*m1+J1*g*lg2-J2*g*l1)*m2^2+(((g*lg1^2-g*l1*lg1)*lg2-g*lg1*lg2^2)*m1^2+(((g*lg1^2-g*l1*lg1)*lg2-g*lg1*lg2^2)*m0+J1*g*lg2-J2*g*lg1-J2*g*l1)*m1+(J1*g*lg2-J2*g*l1)*m0)*m2-J2*g*lg1*m1^2-J2*g*lg1*m0*m1);
ddth2_th2 = ((((g*l1-g*lg1)*lg2^2+(g*lg1^2-2*g*l1*lg1+g*l1^2)*lg2)*m1+(g*l1*lg2^2+g*l1^2*lg2)*m0+J1*g*lg2)*m2^2+((g*lg1^2*lg2*m0+J1*g*lg2)*m1+J1*g*lg2*m0)*m2);
ddth2_F = (((lg1*lg2^2+(l1*lg1-lg1^2)*lg2)*m1-J1*lg2+J2*l1)*m2+J2*lg1*m1);

A = [
    0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 1;
    (ddx0_x0/den), (ddx0_th1/den), (ddx0_th2/den), 0, 0, 0;
    (ddth1_x0/den), (ddth1_th1/den), (ddth1_th2/den), 0, 0, 0;
    (ddth2_x0/den), (ddth2_th1/den), (ddth2_th2/den), 0, 0, 0;
];

B = [
    0;
    0;
    0;
    (ddx0_F/den);
    (ddth1_F/den);
    (ddth2_F/den);
];

C = [
    1, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0;
];

D = [
    0;
    0;
    0;
];

Q = diag( [ 1, 1, 1, 1, 1, 1 ] );
R = diag( [ 1 ] );
gain = lqr( A ,B, Q, R )

% return;
sim( 'DoubleInvertedPendulum' );

%% �`��

figure(1); hold on;
x_pre = l1 * sin( Output(1,3) ) + l2 * sin( Output(1,3) + Output(1,4) );
y_pre = l1 * cos( Output(1,3) ) + l2 * cos( Output(1,3) + Output(1,4) );

for cnt = 1:length( Output )
    clf;
    x0 = Output(cnt,2);
    th1 = Output(cnt,3);
    th2 = Output(cnt,4);
    
    x = [
        x0;    
        x0 + l1 * sin( th1 );
        x0 + l1 * sin( th1 ) + l2 * sin( th1 + th2 );
    ];

    y = [
        0;    
        l1 * cos( th1 );
        l1 * cos( th1 ) + l2 * cos( th1 + th2 );
    ];

    plot( x, y, 'Color', [ 0, 0, 0 ] ); hold on;
    scatter( x, y, 'MarkerFaceColor', [ 0, 0, 0 ], 'MarkerEdgeColor', [ 0, 0, 0 ] );
    
    plot( [ x_pre;x(3); ], [ y_pre;y(3); ], 'Color', [ 1, 0, 0 ] );
    x_pre = x(3);
    y_pre = y(3);
    
    axis equal;
%     axis( [ -3, 3, -3, 3 ] );
    axis( [ -10, 10, -10, 10 ] );
    pause( 0.1 );
end
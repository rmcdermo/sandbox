% McDermott
% 3-7-2016
% curlchem.m
%
% Apply Curl's model to fast chemistry in a hypothetical cell
%
% Consider the following simple reaction scheme:
%
% F1 + A1 --> 2*P1
% F2 + A2 --> 2*P2
% F1 + A2 --> 2*P3

close all
clear all

% initial volume fractions

X_F1_0 = 0.1;
X_F2_0 = 0.4;
X_A2_0 = 0.2;
X_P1_0 = 0.0;
X_P2_0 = 0.0;
X_P3_0 = 0.0;
X_A1_0 = 1 - X_F1_0 - X_F2_0 - X_A2_0 - X_P1_0 - X_P2_0 - X_P3_0;

N = 100000;
N_F1_0 = X_F1_0*N; N1 = N_F1_0;
N_F2_0 = X_F2_0*N; N2 = N_F1_0+N_F2_0;
N_A1_0 = X_A1_0*N; N3 = N_F1_0+N_F2_0+N_A1_0;
N_A2_0 = X_A2_0*N;
N_P1_0 = X_P1_0*N;
N_P2_0 = X_P2_0*N;
N_P3_0 = X_P3_0*N;

N_F1 = N_F1_0;
N_F2 = N_F2_0;
N_A1 = N_A1_0;
N_A2 = N_A2_0;
N_P1 = N_P1_0;
N_P2 = N_P2_0;
N_P3 = N_P3_0;

% add concentration field for ODE integration

Z_F1_0 = X_F1_0;
Z_F2_0 = X_F2_0;
Z_A1_0 = X_A1_0;
Z_A2_0 = X_A2_0;
Z_P1_0 = X_P1_0;
Z_P2_0 = X_P2_0;
Z_P3_0 = X_P3_0;

Z_F1 = Z_F1_0;
Z_F2 = Z_F2_0;
Z_A1 = Z_A1_0;
Z_A2 = Z_A2_0;
Z_P1 = Z_P1_0;
Z_P2 = Z_P2_0;
Z_P3 = Z_P3_0;

% testing new scheme

Y_F1_0 = Z_F1_0;
Y_F2_0 = Z_F2_0;
Y_A1_0 = Z_A1_0;
Y_A2_0 = Z_A2_0;
Y_P1_0 = Z_P1_0;
Y_P2_0 = Z_P2_0;
Y_P3_0 = Z_P3_0;

% assign species to particles, add random position for visualization

%figure(1)
for i = 1:N
    p(i).x = rand(1);
    p(i).y = rand(1);
    if i<=N1
        p(i).spec = 'F1';
        %plot(p(i).x,p(i).y,'bo'); hold on
    elseif i>N1 & i<=N2
        p(i).spec = 'F2';
        %plot(p(i).x,p(i).y,'ko')
    elseif i>N2 & i<=N3
        p(i).spec = 'A1';
        %plot(p(i).x,p(i).y,'g.')
    elseif i>N3
        p(i).spec = 'A2';
        %plot(p(i).x,p(i).y,'b.')
    end
end
%title('initial condition')

iter=0;
n_iter=1000;


dt = .1;

while iter<n_iter

    iter = iter + 1;

    % create a random permutation vector

    % I = randperm(N);

    % % now react pairs of particles

    % for i=2:2:N
    %     % p(i-1) reacts (if possible) with p(i)
    %     if ( p(I(i-1)).spec=='F1' & p(I(i)).spec=='A1' ) | ( p(I(i-1)).spec=='A1' & p(I(i)).spec=='F1' )
    %         p(I(i-1)).spec='P1';
    %         p(I(i)).spec='P1';
    %         N_F1 = N_F1-1;
    %         N_A1 = N_A1-1;
    %         N_P1 = N_P1+2;
    %     elseif ( p(I(i-1)).spec=='F2' & p(I(i)).spec=='A2' ) | ( p(I(i-1)).spec=='A2' & p(I(i)).spec=='F2' )
    %         p(I(i-1)).spec='P2';
    %         p(I(i)).spec='P2';
    %         N_F2 = N_F2-1;
    %         N_A2 = N_A2-1;
    %         N_P2 = N_P2+2;
    %     elseif ( p(I(i-1)).spec=='F1' & p(I(i)).spec=='A2' ) | ( p(I(i-1)).spec=='A2' & p(I(i)).spec=='F1' )
    %         p(I(i-1)).spec='P3';
    %         p(I(i)).spec='P3';
    %         N_F1 = N_F1-1;
    %         N_A2 = N_A2-1;
    %         N_P3 = N_P3+2;
    %     end
    % end

    % % compute new mass fractions

    % X_F1 = N_F1/N;
    % X_F2 = N_F2/N;
    % X_A1 = N_A1/N;
    % X_A2 = N_A2/N;
    % X_P1 = N_P1/N;
    % X_P2 = N_P2/N;
    % X_P3 = N_P3/N;

    % figure(2)
    % plot(iter,X_F1,'bo'); hold on
    % plot(iter,X_F2,'ko')
    % plot(iter,X_A1,'go')
    % plot(iter,X_A2,'mo')
    % plot(iter,X_P1,'co')
    % plot(iter,X_P2,'yo')
    % plot(iter,X_P3,'mo')
    %legend('F1','F2','A1','A2','P1','P2','P3','location','eastoutside')

    % figure(3)
    % for i = 1:N
    %     if p(i).spec=='F1'
    %         plot(p(i).x,p(i).y,'bo'); hold on

    %     elseif p(i).spec=='F2'
    %         plot(p(i).x,p(i).y,'ko'); hold on

    %     elseif p(i).spec=='A1'
    %         plot(p(i).x,p(i).y,'g.'); hold on

    %     elseif p(i).spec=='A2'
    %         plot(p(i).x,p(i).y,'b.'); hold on

    %     elseif p(i).spec=='P1'
    %         plot(p(i).x,p(i).y,'c+'); hold on

    %     elseif p(i).spec=='P2'
    %         plot(p(i).x,p(i).y,'y+'); hold on

    %     elseif p(i).spec=='P3'
    %         plot(p(i).x,p(i).y,'m+'); hold on

    %     end
    % end

    % analogous ODE integration

    Z_F1_0 = Z_F1;
    Z_F2_0 = Z_F2;
    Z_A1_0 = Z_A1;
    Z_A2_0 = Z_A2;
    Z_P1_0 = Z_P1;
    Z_P2_0 = Z_P2;
    Z_P3_0 = Z_P3;

    % EDC-type reaction rate
    % R1 = max(0,min(Z_F1_0,Z_A1_0));
    % R2 = max(0,min(Z_F2_0,Z_A2_0));
    % R3 = max(0,min(Z_F1_0,Z_A2_0));

    % First-order reactions
    R1 = Z_F1_0*Z_A1_0;
    R2 = Z_F2_0*Z_A2_0;
    R3 = Z_F1_0*Z_A2_0;

    Z_F1 = Z_F1_0 + dt*( -R1 - R3 );
    Z_F2 = Z_F2_0 + dt*( -R2 );
    Z_A1 = Z_A1_0 + dt*( -R1 );
    Z_A2 = Z_A2_0 + dt*( -R2 - R3 );
    Z_P1 = Z_P1_0 + dt*( 2*R1 );
    Z_P2 = Z_P2_0 + dt*( 2*R2 );
    Z_P3 = Z_P3_0 + dt*( 2*R3 );

    plot(iter,Z_F1,'b.'); hold on
    plot(iter,Z_F2,'k.')
    plot(iter,Z_A1,'g.')
    plot(iter,Z_A2,'m.')
    % plot(iter,Z_P1,'c.')
    % plot(iter,Z_P2,'y.')
    % plot(iter,Z_P3,'m.')

    pause(0.001)

    % if (N_F1==0 | N_A1==0) & (N_F2==0 | N_A2==0) & (N_F1==0 | N_A2==0)
    %     disp('all reactions complete!')
    %     X_F1 = N_F1/N;
    %     X_F2 = N_F2/N;
    %     X_A1 = N_A1/N;
    %     X_A2 = N_A2/N;
    %     X_P1 = N_P1/N;
    %     X_P2 = N_P2/N;
    %     X_P3 = N_P3/N;
    %     X = [X_F1,X_F2,X_A1,X_A2,X_P1,X_P2,X_P3]
    %     X_Sum = sum(X)
    %     Z = [Z_F1,Z_F2,Z_A1,Z_A2,Z_P1,Z_P2,Z_P3]
    %     break
    % end

end

Z = [Z_F1,Z_F2,Z_A1,Z_A2,Z_P1,Z_P2,Z_P3]

% new scheme test

dt = 10;

R1 = Y_F1_0*Y_A1_0;
R2 = Y_F2_0*Y_A2_0;
R3 = Y_F1_0*Y_A2_0;

Y_F1 = Y_F1_0 + dt*( -R1 - R3 );
Y_F2 = Y_F2_0 + dt*( -R2 );
Y_A1 = Y_A1_0 + dt*( -R1 );
Y_A2 = Y_A2_0 + dt*( -R2 - R3 );
Y_P1 = Y_P1_0 + dt*( 2*R1 );
Y_P2 = Y_P2_0 + dt*( 2*R2 );
Y_P3 = Y_P3_0 + dt*( 2*R3 );

sumY_0 = sum(Y_F1_0+Y_F2_0+Y_A1_0+Y_A2_0+Y_P1_0+Y_P2_0+Y_P3_0);
sumY = sum(Y_F1+Y_F2+Y_A1+Y_A2+Y_P1+Y_P2+Y_P3);

figure
tvec = [0,dt];
plot(tvec,[Y_F1_0,Y_F1],'ko-'); hold on
plot(tvec,[Y_F2_0,Y_F2],'bo-')
plot(tvec,[Y_A1_0,Y_A1],'ro-')
plot(tvec,[Y_A2_0,Y_A2],'mo-')
plot(tvec,[Y_P1_0,Y_P1],'k+-')
plot(tvec,[Y_P2_0,Y_P2],'r+-')
plot(tvec,[Y_P3_0,Y_P3],'b+-')
plot(tvec,[0,0],'k--')
plot(tvec,[1,1],'k--')
plot(tvec,[sumY_0,sumY],'bsq-')

% new scheme

R1 = Y_F1_0*Y_A1_0;
R2 = Y_F2_0*Y_A2_0;
R3 = Y_F1_0*Y_A2_0;

% now, solve for dt_i that sets Y_F1 to zero, the min(dt_i) is our first step

DT(1) = -Y_F1_0/( -R1 - R3 );
DT(2) = -Y_F2_0/( -R2 );
DT(3) = -Y_A1_0/( -R1 );
DT(4) = -Y_A2_0/( -R2 - R3 );

dt = min(DT)
Y_F1 = Y_F1_0 + dt*( -R1 - R3 );
Y_F2 = Y_F2_0 + dt*( -R2 );
Y_A1 = Y_A1_0 + dt*( -R1 );
Y_A2 = Y_A2_0 + dt*( -R2 - R3 );
Y_P1 = Y_P1_0 + dt*( 2*R1 );
Y_P2 = Y_P2_0 + dt*( 2*R2 );
Y_P3 = Y_P3_0 + dt*( 2*R3 );

sumY_0 = sum(Y_F1_0+Y_F2_0+Y_A1_0+Y_A2_0+Y_P1_0+Y_P2_0+Y_P3_0);
sumY = sum(Y_F1+Y_F2+Y_A1+Y_A2+Y_P1+Y_P2+Y_P3);

figure
tvec = [0,dt];
plot(tvec,[Y_F1_0,Y_F1],'bo-'); hold on
plot(tvec,[Y_F2_0,Y_F2],'ko-')
plot(tvec,[Y_A1_0,Y_A1],'go-')
plot(tvec,[Y_A2_0,Y_A2],'mo-')
% plot(tvec,[Y_P1_0,Y_P1],'k+-')
% plot(tvec,[Y_P2_0,Y_P2],'r+-')
% plot(tvec,[Y_P3_0,Y_P3],'b+-')
plot(tvec,[0,0],'k--')
plot(tvec,[1,1],'k--')
plot(tvec,[sumY_0,sumY],'bsq-')

% step 2

Y_F1_0 = Y_F1;
Y_F2_0 = Y_F2;
Y_A1_0 = Y_A1;
Y_A2_0 = Y_A2;
Y_P1_0 = Y_P1;
Y_P2_0 = Y_P2;
Y_P3_0 = Y_P3;

R1 = Y_F1_0*Y_A1_0
R2 = Y_F2_0*Y_A2_0
R3 = Y_F1_0*Y_A2_0

Y = [Y_F1,Y_F2,Y_A1,Y_A2,Y_P1,Y_P2,Y_P3]
























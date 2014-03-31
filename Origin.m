%==============================Impulse Response============================
%             The impulse response is calculated by Monte Carlo method
%                       with a collimated laser source
%                                   2014.3.26
%==========================================================================
% harbor laser 10m 500;12m 500;14m 514;16m 1000 
% harbor LED  div=20 10m 500;12m 500;14m 1000;16 1000
% coastal LED div=40 30m 1000;40m 1000
clear; clc;
load scattering_phase_table_coastal.mat;

% Discard the historical data
fid = fopen('impulse_response_z_40_coastal_t_LED_div_40.bin','w');
fclose(fid);
fid = fopen('impulse_response_z_40_coastal_W_LED_div_40.bin','w');
fclose(fid);
fid = fopen('impulse_response_z_40_coastal_u_z_LED_div_40.bin','w');
fclose(fid);
fid = fopen('impulse_response_z_40_coastal_t_LED_div_40.idx','w');
fclose(fid);


% Parameters settings
N = 1e6;           % the number of the photons
z_rec =40;          % the location of the receiver m 30 40 50 60
W_th = 1e-6;
a = 0.179;  % absorbing coefficient m^-1
b = 0.219;  % scattering coefficient m^-1
c = 0.398;  % attenuation coefficient m^-1
v = 2.25257e8;  % the speed of light in the sea m/s
aperture = 0.5;     % the diameter of the receiver m
t_0 = z_rec/v;
LoopNum = 1000;

for n = 1:LoopNum
    tic;
    x = zeros(1,N);     % initialize the coordinates m
    y = zeros(1,N);
    z = zeros(1,N);
    t = zeros(1,N);     % initialize the time s
    W = ones(1,N);      % initialize the power W

%===============================Source=====================================
    theta_s = 10*pi/180;    %coastal 
    random_theta = unifrnd(zeros(1,N),ones(1,N));
    random_phi = unifrnd(zeros(1,N),ones(1,N));
    theta_incident = acos(1-random_theta.*(1-cos(theta_s/2)));   % transmitting angle rad
    phi_incident = 2*pi*random_phi;
    u_x = sin(theta_incident).*cos(phi_incident);
    u_y = sin(theta_incident).*sin(phi_incident);
    u_z = cos(theta_incident);

    t_i = [];
    W_i = [];
    u_z_i = [];
    
    while ~isempty(z)
        number = length(x);
        random_z = unifrnd(zeros(1,number), ones(1,number));
        delta_z = -log(random_z)/c;                              % something happens between photons and water
        x = x+delta_z.*u_x;     % update coordinates
        y = y+delta_z.*u_y;
        z = z+delta_z.*u_z;
        t = t+delta_z/v;        % update time
        W = W.*(1-a/c);
        
        z_receive = z(z>=z_rec);                           % photons reach the receiving plane 
        x_receive = x(z>=z_rec) - (z_receive-z_rec)./u_z(z>=z_rec).*u_x(z>=z_rec); % fix the reaching coordinates
        y_receive = y(z>=z_rec) - (z_receive-z_rec)./u_z(z>=z_rec).*u_y(z>=z_rec);
        t_receive = t(z>=z_rec) - (z_receive-z_rec)./u_z(z>=z_rec)/v;
        u_z_receive = u_z(z>=z_rec);
        
        W_receive = W(z>=z_rec)/(1-a/c);
        
        z_receive = z_receive(W_receive>=W_th);     % the photons whose weights are higher than the threshold
        x_receive = x_receive(W_receive>=W_th);
        y_receive = y_receive(W_receive>=W_th);
        t_receive = t_receive(W_receive>=W_th);
        u_z_receive = u_z_receive(W_receive>=W_th);
        W_receive = W_receive(W_receive>=W_th);
        
        position = x_receive.^2+y_receive.^2;               % the parameter of the detected photons
        order_detect = position <= (aperture/2).^2;   % decide whether the photon have be detected
        t_detect = t_receive(order_detect);                  % the parameter of the detected photons
        W_detect = W_receive(order_detect);
        u_z_detect = u_z_receive(order_detect);
        
        t_i = [t_i, t_detect];                                  % the last output
        W_i = [W_i, W_detect];
        u_z_i = [u_z_i, u_z_detect];
        
        label_1 = W>=W_th;
        x = x(label_1);           % the photons whose weights aren't below the threshold
        y = y(label_1);
        z = z(label_1);
        t = t(label_1);
        W = W(label_1);
        u_x = u_x(label_1);
        u_y = u_y(label_1);
        u_z = u_z(label_1);
        label_2 = z<z_rec;
        x = x(label_2);           % the photons that haven't reached the receiving plane
        y = y(label_2);
        z = z(label_2);
        t = t(label_2);
        W = W(label_2);
        u_x = u_x(label_2);
        u_y = u_y(label_2);
        u_z = u_z(label_2);
        
        if ~isempty(z)
            random_theta_scatter = unifrnd(zeros(1,length(u_x)), ones(1,length(u_x)));                    % the process of scattering
            random_phi_scatter = unifrnd(zeros(1,length(u_x)), ones(1,length(u_x)));
            u_scatter = cos(interp1(scattering_table, theta_table, random_theta_scatter));
            phi_scatter = 2*pi*random_phi_scatter;

            u = (sqrt(1-u_scatter.^2))./(sqrt((1-u_z.^2)));
            u_x_1 = u.*(u_x.*u_z.*cos(phi_scatter)-u_y.*sin(phi_scatter))+u_x.*u_scatter;
            u_y_1 = u.*(u_y.*u_z.*cos(phi_scatter)+u_x.*sin(phi_scatter))+u_y.*u_scatter;
            u_z_1 = -(sqrt(1-u_scatter.^2)).*cos(phi_scatter).*(sqrt((1-u_z.^2)))+u_z.*u_scatter;
            
            u_x = u_x_1;
            u_y = u_y_1;
            u_z = u_z_1;
            clear u_x_1 u_y_1 u_z_1 u;
        else
            break;
        end
    end
    
    
    % low-level file I/O
    
    % open the file to append and read
    fid1 = fopen('impulse_response_z_40_coastal_t_LED_div_40.bin','a');
    fid2 = fopen('impulse_response_z_40_coastal_W_LED_div_40.bin','a');
    fid3 = fopen('impulse_response_z_40_coastal_u_z_LED_div_40.bin','a');
    
    fid4 = fopen('impulse_response_z_40_coastal_t_LED_div_40.idx','a');


    % write values at end of file
    fwrite(fid1, t_i, 'double');
    fwrite(fid2, W_i, 'double');
    fwrite(fid3, u_z_i, 'double');
    fwrite(fid4, length(t_i), 'int32');

    % close the file
    fclose(fid4);
    fclose(fid3);
    fclose(fid2);
    fclose(fid1);
    
    
    clearvars -except LoopNum n N z_rec W_th v t_0 a b c h aperture theta_table scattering_table;
    
    display(sprintf('%d/%d loops.', n, LoopNum));
    toc;
end


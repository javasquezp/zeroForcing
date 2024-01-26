%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Written by Juan A. Vasquez-Peralvo
%Luxembourg-Luxembourg
%5/12/2023
%Modified 7/28/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%Description%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code allows tto obtain the radiation pattern of 7 beams considering
%constrains of power, beamwidth, SLL and activation elements

%%%%%%%%%%%%%%%%%%%%%%%%Initialize Matlab%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Add Required Paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------------------
%                            Main Script Logic
% ------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
fprintf( 1, '------------------------------------------------------------\n' );
fprintf( 1, '------ Zero Forcing  --------\n' );
fprintf( 1, '------------------------------------------------------------\n' );
fprintf( 1, '------ Author(s):                                   --------\n' );
fprintf( 1, '------          Juan Andrés Vásquez Peralvo         --------\n' );
fprintf( 1, '------------------------------------------------------------\n' );
fprintf( 1, '------ Date:   January 2024                         --------\n' );
fprintf( 1, '------------------------------------------------------------\n' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start Timer
tini=clock;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Initiate constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Nx = 10;
    Ny = 10;
    dx = 0.5;
    dy = 0.5;
    desiredAzimuth = 10;
    desiredElevation = -10;
    nullAzimuth = [-16 0 10];
    nullElevation = [-7 0 6];
% Define the range of azimuth and elevation angles
    azRange = -180:0.5:180; % Azimuth range in degrees
    elRange = -90:0.5:90;   % Elevation range in degrees
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Initial Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert angles to radians
    desiredAzimuthRad = deg2rad(desiredAzimuth);
    desiredElevationRad = deg2rad(desiredElevation);
    nullAzimuthRad = deg2rad(nullAzimuth);
    nullElevationRad = deg2rad(nullElevation);

% Wave number
    k = 2 * pi;

% Generate the URA matrix
    [x, y] = meshgrid(0:Nx-1, 0:Ny-1);
    x = x * dx;
    y = y * dy;

% Steering vector for the desired direction
    aDesired = exp(-1i * k * (x * sin(desiredAzimuthRad) * cos(desiredElevationRad) + y * sin(desiredElevationRad))).';

% Steering vector for the null direction
for m = 1 : length(nullAzimuthRad)
    aNull(:, :, m) = exp(-1i * k * (x * sin(nullAzimuthRad(m)) * cos(nullElevationRad(m)) + y * sin(nullElevationRad(m)))).';
end
aNull = reshape(aNull, Nx*Ny, length(nullAzimuthRad));                              
% Channel matrix
    H = [aDesired(:) aNull];

% Zero forcing beamforming weights calculation
    weights = pinv(H)' * [1; zeros(length(nullAzimuthRad), 1)];

% Calculate the beam pattern
    for az_idx = 1:length(azRange)
        for el_idx = 1:length(elRange)
            az_rad = deg2rad(azRange(az_idx));
            el_rad = deg2rad(elRange(el_idx));
            steering_vector = exp(-1i * k * (x * sin(az_rad) * cos(el_rad) + y * sin(el_rad))).';
            beam_pattern(el_idx, az_idx) = abs(weights' * steering_vector(:))^2;
        end
    end

% Normalize the beam pattern
beam_pattern = beam_pattern / max(beam_pattern(:));
% dB the pattern
beam_patterndB = pow2db(beam_pattern);
% Plot
[AZ, EL] = meshgrid(azRange, elRange);
surf(AZ, EL, beam_patterndB, 'LineStyle', 'none');
xlabel('Azimuth (degrees)');
ylabel('Elevation (degrees)');
zlabel('Normalized Power');
title('Beam Pattern of URA');
caxis([ min(pow2db(abs(beam_patterndB(:))+0.0001)) max(beam_patterndB(:))])
xlim([-90 90])
ylim([-90 90])
view(2); % 2D view
hold on
scatter3(nullAzimuth, nullElevation, max(beam_patterndB(:)).*ones(1, length(nullAzimuthRad)), 'r*')
legend('Pattern', 'Nulls')
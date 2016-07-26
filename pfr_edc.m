% McDermott
% 11-30-2015
% pfr_edc.m
%
% This script solves a set of ODEs for a plug flow reactor model with EDC as the combustion model.
%
%  --------------------------------
%
%      m_dot, cp, T, Y_fuel =>
%
%  --------------------------------

close all
clear all

tau_mix        = 0.01; % mixing time (s)
mass_flux_fuel = 0.04; % mass flux of fuel (kg/m2/s)
N_A            = 2/0.21; % molar stoich coefficient of air
W_A            = 29; % mole wt of air
W_F            = 16; % mole wt of methane
s              = N_A*W_A/W_F; % mass stoich coef of air
cp             = 1000; % specific heat, (J/kg/K)
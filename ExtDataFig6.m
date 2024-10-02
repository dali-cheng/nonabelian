clearvars; clc; close all;
% Run Fig4b_experiment.m first to generate the data files
% 'Data_ExtDataFig6_1.mat' and 'Data_ExtDataFig6_2.mat'
load('Data_ExtDataFig6_2.mat'); % Change as needed
equirect_project([S(:, 3), S(:, 1), S(:, 2)]);
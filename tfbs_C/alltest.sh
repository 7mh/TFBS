#!/bin/bash
######################################

# Filtering and Brute force Mutation Alignment for all data files.

######################################
./tf8merFilter 0.3 < ../data_files/Jumeau_8mers_11111111.txt | ./tf8merAlign ../data_files/CG5080-WT.asc

./tf8merFilter 0.3 < ../data_files/ Twi_8mers_11111111.txt | ./tf8merAlign ../data_files/CG5080-WT.asc
./tf8merFilter 0.3 < ../data_files/Tup_8mers_11111111.txt | ./tf8merAlign ../data_files/CG5080-WT.asc
./tf8merFilter 0.3 < ../data_files/Tin_8mers_11111111.txt | ./tf8merAlign ../data_files/CG5080-WT.asc
./tf8merFilter 0.3 < ../data_files/ Pnt_8mers_11111111.txt | ./tf8merAlign ../data_files/CG5080-WT.asc
./tf8merFilter 0.3 < ../data_files/CHES-1-like_8mers_11111111.txt | ./tf8merAlign ../data_files/CG5080-WT.asc

# File containing the appropriate options for
# the gfortran compiler

FC = gfortran
CC = gcc
FCFLAGS = -O2 -Wall
DEBUG_FCFLAGS = -g -O0 -fbounds-check -Wall -Wextra -Wconversion -pedantic -ffpe-trap=zero,overflow,underflow


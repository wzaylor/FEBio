# Make include file for FEBio on Linux

include $(FEBDIR)build/lnx64.mk

CC = g++

FLG = -g -O0 -fPIC -fopenmp -std=c++11

INC = -I$(FEBDIR) -I$(FEBDIR)build/include

#Makefile for acc-lanczos

CXX=		nvc++

#only use for lanczos_std, lanczos_cpu and lanczos_gpu
#CXX=		g++

CPPFLAGS=	-MMD -MP
OPT_FLAG=	-O3
ACC_FLAGS=	-acc -Minfo=accel
CXXERR_FLAGS=	-W -Wall

#Possible combination of flags

#CPU version
CPU_CPLEX=	-DCOMPLEX3D

#mvmul accelerated
GPU_CPLEX=	-DCOMPLEX3D -DPARALLEL_1

#Lanczos algorithm accelerated
ACC_CPLEX=	-DCOMPLEX3D -DPARALLEL_2 -DACC_ALG

SRCDIR=		src/
OBJDIR=		obj/
INCDIRS=	-I inc/custom_lqcd -I inc/lambda_lanczos

TARGETS= 	lanczos_std lanczos_cpu lanczos_gpu lanczos_acc

.PHONY: all full clean .gitignore

all:	$(TARGETS) .gitignore

full:	clean all


lanczos_std:		$(OBJDIR)grid_std.o $(OBJDIR)main_std.o
	$(CXX) $(LINK_NVTX) -o $@ $^

lanczos_cpu:		$(OBJDIR)grid_cpu.o $(OBJDIR)main_cpu.o
	$(CXX) $(LINK_NVTX) -o $@ $^

lanczos_gpu:		$(OBJDIR)grid_gpu.o $(OBJDIR)main_gpu.o
	$(CXX) $(LINK_NVTX) $(ACC_FLAGS) -o $@ $^

lanczos_acc:		$(OBJDIR)grid_acc.o $(OBJDIR)main_acc.o
	$(CXX) $(LINK_NVTX) $(ACC_FLAGS) -o $@ $^



$(OBJDIR)grid_std.o:		$(SRCDIR)grid.cpp
	g++ -E -MMD -MP -MF $(OBJDIR)$(*F).d -MT $(OBJDIR)$(*F).o $(INCDIRS) $< >/dev/null 2>&1
	$(CXX) $(OPT_FLAG) $(INCDIRS) $(CXXERR_FLAGS) -o $@ -c $<

$(OBJDIR)main_std.o:		$(SRCDIR)old_main.cpp
	g++ -E -MMD -MP -MF $(OBJDIR)$(*F).d -MT $(OBJDIR)$(*F).o $(INCDIRS) $< >/dev/null 2>&1
	$(CXX) $(OPT_FLAG) $(INCDIRS) $(CXXERR_FLAGS) -o $@ -c $<

$(OBJDIR)grid_cpu.o:		$(SRCDIR)grid.cpp
	g++ -E -MMD -MP $(CPU_CPLEX) -MF $(OBJDIR)$(*F).d -MT $(OBJDIR)$(*F).o $(INCDIRS) $< >/dev/null 2>&1
	$(CXX) $(OPT_FLAG) $(INCDIRS) $(CXXERR_FLAGS) $(CPU_CPLEX) -o $@ -c $<

$(OBJDIR)main_cpu.o:		$(SRCDIR)main.cpp
	g++ -E -MMD -MP $(CPU_CPLEX) -MF $(OBJDIR)$(*F).d -MT $(OBJDIR)$(*F).o $(INCDIRS) $< >/dev/null 2>&1
	$(CXX) $(OPT_FLAG) $(INCDIRS) $(CXXERR_FLAGS) $(CPU_CPLEX) -o $@ -c $<

$(OBJDIR)grid_gpu.o:		$(SRCDIR)grid.cpp
	g++ -E -MMD -MP $(GPU_CPLEX) -MF $(OBJDIR)$(*F).d -MT $(OBJDIR)$(*F).o $(INCDIRS) $< >/dev/null 2>&1
	$(CXX) $(OPT_FLAG) $(ACC_FLAGS) $(INCDIRS) $(CXXERR_FLAGS) $(GPU_CPLEX) -o $@ -c $<

$(OBJDIR)main_gpu.o:		$(SRCDIR)main.cpp
	g++ -E -MMD -MP $(GPU_CPLEX) -MF $(OBJDIR)$(*F).d -MT $(OBJDIR)$(*F).o $(INCDIRS) $< >/dev/null 2>&1
	$(CXX) $(OPT_FLAG) $(ACC_FLAGS) $(INCDIRS) $(CXXERR_FLAGS) $(GPU_CPLEX) -o $@ -c $<

$(OBJDIR)grid_acc.o:		$(SRCDIR)grid.cpp
	g++ -E -MMD -MP $(ACC_CPLEX) -MF $(OBJDIR)$(*F).d -MT $(OBJDIR)$(*F).o $(INCDIRS) $< >/dev/null 2>&1
	$(CXX) $(OPT_FLAG) $(ACC_FLAGS) $(INCDIRS) $(CXXERR_FLAGS) $(ACC_CPLEX) -o $@ -c $<

$(OBJDIR)main_acc.o:		$(SRCDIR)main.cpp
	g++ -E -MMD -MP $(ACC_CPLEX) -MF $(OBJDIR)$(*F).d -MT $(OBJDIR)$(*F).o $(INCDIRS) $< >/dev/null 2>&1
	$(CXX) $(OPT_FLAG) $(ACC_FLAGS) $(INCDIRS) $(CXXERR_FLAGS) $(ACC_CPLEX) -o $@ -c $<



.gitignore:
	echo "obj/*.o" "obj/*.d" ${TARGETS} | xargs -n1 > .gitignore

-include $(OBJDIR)*.d

clean:
	rm -f $(OBJDIR)*.o $(OBJDIR)*.d $(TARGETS) 

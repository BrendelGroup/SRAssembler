
RM := rm -rf
DFLAG:=
CPPFLAG:=
LDFlag:=
PROGDIR:=../bin
PROG_NO_MPI:=$(PROGDIR)/SRAssembler
PROG_MPI:=$(PROGDIR)/SRAssembler_MPI
PROG:=$(PROG_NO_MPI)
INSTALLDIR:=/usr/local/bin
MPICXX:=mpicxx
MPIFLAG:=-DMPIMODE
ifdef with-boost
  CPPFLAG:=-I$(with-boost)
endif

# All of the sources participating in the build are defined here
CPP_SRCS += \
AbyssAssembler.cpp \
Aligner.cpp \
Assembler.cpp \
GSQAligner.cpp \
GTHAligner.cpp \
ExonerateAligner.cpp \
GeneFinder.cpp \
SnapGeneFinder.cpp \
Logger.cpp \
MPIWrapper.cpp \
SOAPDenovoAssembler.cpp \
SRAssembler.cpp \
SRAssemblerMaster.cpp \
SRAssemblerSlave.cpp \
SplicedAligner.cpp \
Utility.cpp \
VmatchBWAfastmapAligner.cpp \
Library.cpp 

OBJS += \
AbyssAssembler.o \
Aligner.o \
Assembler.o \
GSQAligner.o \
GTHAligner.o \
ExonerateAligner.o \
GeneFinder.o \
SnapGeneFinder.o \
Logger.o \
MPIWrapper.o \
SOAPDenovoAssembler.o \
SRAssembler.o \
SRAssemblerMaster.o \
SRAssemblerSlave.o \
SplicedAligner.o \
Utility.o \
VmatchBWAfastmapAligner.o \
Library.o 

CPP_DEPS += \
AbyssAssembler.d \
Aligner.d \
Assembler.d \
GSQAligner.d \
GTHAligner.d \
ExonerateAligner.d \
GeneFinder.d \
SnapGeneFinder.d \
Logger.d \
MPIWrapper.d \
SOAPDenovoAssembler.d \
SRAssembler.d \
SRAssemblerMaster.d \
SRAssemblerSlave.d \
SplicedAligner.d \
Utility.d \
VmatchBWAfastmapAligner.d \
Library.d \


# Each subdirectory must supply rules for building sources it contributes
%.o: %.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(DFLAG) $(CPPFLAG) -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
endif

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: SRAssembler

mpi: CXX=$(MPICXX)
mpi: DFLAG+=$(MPIFLAG)
mpi: PROG=$(PROG_MPI)
mpi: SRAssembler
# Tool invocations
SRAssembler: $(OBJS) $(USER_OBJS)
	mkdir -p $(PROGDIR)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	$(CXX)  -o"$(PROG)" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
all-clean:
	-$(RM) $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) $(PROG_NO_MPI) $(PROG_MPI)
	-@echo ' '
clean:
	-$(RM) $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS)
	-@echo ' '
install:
	cp ../bin/* $(INSTALLDIR)/


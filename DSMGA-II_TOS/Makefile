# ---------------------------------------------------------------------------- #
#define TRIMMING
#define MEDIAN
#define MEAN
#define MAX
# ---------------------------------------------------------------------------- #
#define ORDERING

#define EPOCH_RESTART
#define GENERATION_RESTART

#define COUNTING
#define MOVING_ONE
#define MOVING_FRONT
# ---------------------------------------------------------------------------- #
#define SIMILARITY_CHECK
# ---------------------------------------------------------------------------- #
CXX = g++
#CXXFLAGS = -O0 -g -std=c++11
CXXFLAGS = -O2 -Wall -march=native -std=c++11
INCLUDE = 
TLIB = -lm

#-----Suffix Rules---------------------------
# set up C++ suffixes and relationship between .cc and .o files

.SUFFIXES: .cpp .otttmeanmoes .otttmeanmogs .otttmeances .otttmeancgs

.o :
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

.cpp :
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< -o $@ -lm $(TLIB) 

# .cpp.ottt : 
# 	$(CXX) -DTRIMMING -DORDERING -DSIMILARITY_CHECK $(CXXFLAGS) $(INCLUDE) -c $< -o $@

.cpp.otttmeanmoes:
	$(CXX) -DTRIMMING -DMEAN -DORDERING -DMOVING_ONE -DEPOCH_RESTART -DSIMILARITY_CHECK $(CXXFLAGS) $(INCLUDE) -c $< -o $@
.cpp.otttmeanmogs:
	$(CXX) -DTRIMMING -DMEAN -DORDERING -DMOVING_ONE -DGENERATION_RESTART -DSIMILARITY_CHECK $(CXXFLAGS) $(INCLUDE) -c $< -o $@
.cpp.otttmeances:
	$(CXX) -DTRIMMING -DMEAN -DORDERING -DCOUNTING -DEPOCH_RESTART -DSIMILARITY_CHECK $(CXXFLAGS) $(INCLUDE) -c $< -o $@
.cpp.otttmeancgs:
	$(CXX) -DTRIMMING -DMEAN -DORDERING -DCOUNTING -DGENERATION_RESTART -DSIMILARITY_CHECK $(CXXFLAGS) $(INCLUDE) -c $< -o $@
	
#-----File Dependencies----------------------

SRC = $(OBJ_NOFLAG) $(SRC_FLAG) $(SRC_Zobrist)

SRC_FLAG = chromosome.cpp dsmga2.cpp global.cpp
SRC_NOFLAG = fastcounting.cpp main.cpp mt19937ar.cpp myrand.cpp spin.cpp nk-wa.cpp sat.cpp maxcut.cpp
SRC_Zobrist = genZobrist.cpp

OBJ = $(addsuffix .o, $(basename $(SRC))) $(OBJ_TTT_MEAN_MOE_S) $(OBJ_TTT_MEAN_MOG_S) $(OBJ_TTT_MEAN_CE_S) $(OBJ_TTT_MEAN_CG_S)
OBJ_NOFLAG = $(addsuffix .o, $(basename $(SRC_NOFLAG)))

OBJ_Z = $(addsuffix .o, $(basename $(SRC_Zobrist)))

OBJ_O = $(addsuffix .o, $(basename $(SRC_FLAG)))

# all: DSMGA2 dsmga2-tos_ttt genZobrist  # sweep
all: DSMGA2 TTT_MEAN_MOE_S TTT_MEAN_MOG_S TTT_MEAN_CE_S TTT_MEAN_CG_S genZobrist  # sweep

DSMGA2: $(OBJ_NOFLAG) $(OBJ_O)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_NOFLAG) $(OBJ_O)

genZobrist: $(OBJ_Z)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_Z)

# ---------------------------------------------------------------------------- #
#                                      FFF                                     #
# ---------------------------------------------------------------------------- #

# - dsmga2-tos-fff

# ---------------------------------------------------------------------------- #
#                                      TFF                                     #
# ---------------------------------------------------------------------------- #

# - dsmga2-tos-tft
# - dsmga2-tos-tff_max
# - dsmga2-tos-tff_median

# ---------------------------------------------------------------------------- #
#                                      FTF                                     #
# ---------------------------------------------------------------------------- #

# - dsmga2-tos-ftf_moE
# - dsmga2-tos-ftf_moG
# - dsmga2-tos-ftf_moN
# - dsmga2-tos-tff_cntE
# - dsmga2-tos-ftf_cntG
# - dsmga2-tos-ftf_cntN
# - dsmga2-tos-ftf_mfE
# - dsmga2-tos-ftf_mfG
# - dsmga2-tos-ftf_mfN
# - dsmga2-tos-ftf_cntE

# ---------------------------------------------------------------------------- #
#                                      FFT                                     #
# ---------------------------------------------------------------------------- #

# - dsmga2-tos-fft

# ---------------------------------------------------------------------------- #
#                                      TTF                                     #
# ---------------------------------------------------------------------------- #

# - dsmga2-tos-ttf_moG
# - dsmga2-tos-ttf_cntE
# - dsmga2-tos-ttf_cntG
# - dsmga2-tos-ttf_moE

# ---------------------------------------------------------------------------- #
#                                      TFT                                     #
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
#                                      FTT                                     #
# ---------------------------------------------------------------------------- #

# - dsmga2-tos-ftt_moG
# - dsmga2-tos-ftt_moE
# - dsmga2-tos-ftt_count

# ---------------------------------------------------------------------------- #
#                                      TTT                                     #
# ---------------------------------------------------------------------------- #

# --------------------------- - dsmga2-tos-ttt_moE --------------------------- #
OBJ_TTT_MEAN_MOE_S = $(addsuffix .otttmeanmoes, $(basename $(SRC_FLAG)))
TTT_MEAN_MOE_S: $(OBJ_NOFLAG) $(OBJ_TTT_MEAN_MOE_S)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_NOFLAG) $(OBJ_TTT_MEAN_MOE_S)
# --------------------------- - dsmga2-tos-ttt_moG --------------------------- #
OBJ_TTT_MEAN_MOG_S = $(addsuffix .otttmeanmogs, $(basename $(SRC_FLAG)))
TTT_MEAN_MOG_S: $(OBJ_NOFLAG) $(OBJ_TTT_MEAN_MOG_S)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_NOFLAG) $(OBJ_TTT_MEAN_MOG_S)
# --------------------------- - dsmga2-tos-ttt_cntE -------------------------- #
OBJ_TTT_MEAN_CE_S = $(addsuffix .otttmeances, $(basename $(SRC_FLAG)))
TTT_MEAN_CE_S: $(OBJ_NOFLAG) $(OBJ_TTT_MEAN_CE_S)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_NOFLAG) $(OBJ_TTT_MEAN_CE_S)
# --------------------------- - dsmga2-tos-ttt_cntG -------------------------- #
OBJ_TTT_MEAN_CG_S = $(addsuffix .otttmeancgs, $(basename $(SRC_FLAG)))
TTT_MEAN_CG_S: $(OBJ_NOFLAG) $(OBJ_TTT_MEAN_CG_S)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_NOFLAG) $(OBJ_TTT_MEAN_CG_S)

#-----Other stuff----------------------------
depend:
	makedepend -Y. $(SRC)

clean:
	rm -f $(OBJ)

# DO NOT DELETE

chromosome.o: spin.h chromosome.h global.h myrand.h mt19937ar.h
chromosome.o: bitwisedistance.h nk-wa.h doublelinkedlistarray.h zkey.h sat.h maxcut.h
dsmga2.o: chromosome.h global.h myrand.h mt19937ar.h bitwisedistance.h
dsmga2.o: spin.h nk-wa.h doublelinkedlistarray.h zkey.h sat.h dsmga2.h maxcut.h
dsmga2.o: statistics.h trimatrix.h fastcounting.h
fastcounting.o: global.h myrand.h mt19937ar.h bitwisedistance.h spin.h
fastcounting.o: nk-wa.h doublelinkedlistarray.h zkey.h sat.h fastcounting.h maxcut.h
global.o: myrand.h mt19937ar.h statistics.h doublelinkedlistarray.h zkey.h
global.o: chromosome.h global.h bitwisedistance.h spin.h nk-wa.h sat.h maxcut.h
main.o: statistics.h dsmga2.h chromosome.h global.h myrand.h mt19937ar.h
main.o: bitwisedistance.h spin.h nk-wa.h doublelinkedlistarray.h zkey.h
main.o: sat.h trimatrix.h fastcounting.h maxcut.h
myrand.o: myrand.h mt19937ar.h
spin.o: global.h myrand.h mt19937ar.h bitwisedistance.h spin.h
spin.o: nk-wa.h doublelinkedlistarray.h zkey.h sat.h maxcut.h
nk-wa.o: nk-wa.h
sat.o: sat.h 
maxcut.o: maxcut.h
chromosome.o: spin.h chromosome.h global.h myrand.h mt19937ar.h
chromosome.o: bitwisedistance.h nk-wa.h doublelinkedlistarray.h zkey.h sat.h maxcut.h
dsmga2.o: chromosome.h global.h myrand.h mt19937ar.h bitwisedistance.h
dsmga2.o: spin.h nk-wa.h doublelinkedlistarray.h zkey.h sat.h dsmga2.h maxcut.h
dsmga2.o: statistics.h trimatrix.h fastcounting.h
fastcounting.o: global.h myrand.h mt19937ar.h bitwisedistance.h spin.h
fastcounting.o: nk-wa.h doublelinkedlistarray.h zkey.h sat.h fastcounting.h maxcut.h
global.o: myrand.h mt19937ar.h statistics.h doublelinkedlistarray.h zkey.h
global.o: chromosome.h global.h bitwisedistance.h spin.h nk-wa.h sat.h maxcut.h
# sweep.o: statistics.h dsmga2.h chromosome.h global.h myrand.h mt19937ar.h
# sweep.o: bitwisedistance.h spin.h nk-wa.h doublelinkedlistarray.h
# sweep.o: zkey.h sat.h trimatrix.h fastcounting.h maxcut.h
myrand.o: myrand.h mt19937ar.h
spin.o: global.h myrand.h mt19937ar.h bitwisedistance.h spin.h
spin.o: nk-wa.h doublelinkedlistarray.h zkey.h sat.h maxcut.h
nk-wa.o: nk-wa.h
sat.o: sat.h
maxcut.o: maxcut.h
nono.o: nono.h

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
# CXXFLAGS = -O2 -Wall -march=native -std=c++11
CXXFLAGS = -O2 -Wall -std=c++11
INCLUDE = 
TLIB = -lm

#-----Suffix Rules---------------------------
# set up C++ suffixes and relationship between .cc and .o files

# ---------------------------------------------------------------------------- #
#                                STEP 1 SUFFFIX                                #
# ---------------------------------------------------------------------------- #

.SUFFIXES: .cpp .otttmeancgs .otffmean .oftfcg
# .otttmeancgs

# -------------------------------- STEP 1 ENDS ------------------------------- #

.o :
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

.cpp :
	$(CXX) $(CXXFLAGS) $(INCLUDE) $< -o $@ -lm $(TLIB) 

# ---------------------------------------------------------------------------- #
#                                 STEP 2 .CPP.O                                #
# ---------------------------------------------------------------------------- #

.cpp.oftfce:
	$(CXX) -DORDERING -DCOUNTING -DEPOCH_RESTART $(CXXFLAGS) $(INCLUDE) -c $< -o $@

.cpp.oftfcg:
	$(CXX) -DORDERING -DCOUNTING -DGENERATION_RESTART $(CXXFLAGS) $(INCLUDE) -c $< -o $@

.cpp.otffmean:
	$(CXX) -DTRIMMING -DMEAN $(CXXFLAGS) $(INCLUDE) -c $< -o $@

.cpp.otffmedian:
	$(CXX) -DTRIMMING -DMEDIAN $(CXXFLAGS) $(INCLUDE) -c $< -o $@

.cpp.otffmax:
	$(CXX) -DTRIMMING -DMAX $(CXXFLAGS) $(INCLUDE) -c $< -o $@

.cpp.ofttmfe : 
	$(CXX) -DORDERING -DMOVING_FRONT -DEPOCH_RESTART $(CXXFLAGS) $(INCLUDE) -c $< -o $@

.cpp.offts : 
	$(CXX) -DSIMILARITY_CHECK $(CXXFLAGS) $(INCLUDE) -c $< -o $@

.cpp.otftmedians : 
	$(CXX) -DTRIMMING -DMEDIAN -DSIMILARITY_CHECK $(CXXFLAGS) $(INCLUDE) -c $< -o $@

.cpp.otftmeans : 
	$(CXX) -DTRIMMING -DMEAN -DSIMILARITY_CHECK $(CXXFLAGS) $(INCLUDE) -c $< -o $@

.cpp.ofttcgs : 
	$(CXX) -DORDERING -DCOUNTING -DGENERATION_RESTART -DSIMILARITY_CHECK $(CXXFLAGS) $(INCLUDE) -c $< -o $@

.cpp.ottfmeancg :
	$(CXX) -DTRIMMING -DMEAN -DORDERING -DCOUNTING -DGENERATION_RESTART $(CXXFLAGS) $(INCLUDE) -c $< -o $@

.cpp.otttmeancgs : 
	$(CXX) -DTRIMMING -DMEAN -DORDERING -DCOUNTING -DGENERATION_RESTART -DSIMILARITY_CHECK $(CXXFLAGS) $(INCLUDE) -c $< -o $@

# -------------------------------- STEP 2 ENDS ------------------------------- #
	
#-----File Dependencies----------------------

SRC = $(OBJ_NOFLAG) $(SRC_FLAG) $(SRC_Zobrist)

SRC_FLAG = chromosome.cpp dsmga2.cpp global.cpp
SRC_NOFLAG = fastcounting.cpp main.cpp mt19937ar.cpp myrand.cpp spin.cpp nk-wa.cpp sat.cpp maxcut.cpp traps.cpp
SRC_Zobrist = genZobrist.cpp

# ---------------------------------------------------------------------------- #
#                                  STEP 3 OBJ                                  #
# ---------------------------------------------------------------------------- #

OBJ = $(addsuffix .o, $(basename $(SRC))) $(OBJ_TTT_MEAN_CG_S) $(OBJ_FTF_CG) $(OBJ_TFF_MEAN)

# -------------------------------- STEP 3 ENDS ------------------------------- #

OBJ_NOFLAG = $(addsuffix .o, $(basename $(SRC_NOFLAG)))

OBJ_Z = $(addsuffix .o, $(basename $(SRC_Zobrist)))

OBJ_O = $(addsuffix .o, $(basename $(SRC_FLAG)))

# ---------------------------------------------------------------------------- #
#                                  STEP 4 ALL                                  #
# ---------------------------------------------------------------------------- #

# all: DSMGA2 TTT_MEAN_CG_S genZobrist  # sweep
all: DSMGA2 TFF_MEAN FTF_CG genZobrist  # sweep
# all: DSMGA2 TFF_MEAN genZobrist  # sweep
# all: FTF_CG genZobrist  # sweep


# -------------------------------- STEP 4 ENDS ------------------------------- #

DSMGA2: $(OBJ_NOFLAG) $(OBJ_O)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_NOFLAG) $(OBJ_O)

genZobrist: $(OBJ_Z)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_Z)

# ---------------------------------------------------------------------------- #
#                                   STEP 5 GA                                  #
# ---------------------------------------------------------------------------- #


# ------------------------------------ FFF ----------------------------------- #
# - dsmga2-tos-fff
# ---------------------------------------------------------------------------- #


# ------------------------------------ TFF ----------------------------------- #
OBJ_TFF_MEAN = $(addsuffix .otffmean, $(basename $(SRC_FLAG)))
TFF_MEAN: $(OBJ_NOFLAG) $(OBJ_TFF_MEAN)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_NOFLAG) $(OBJ_TFF_MEAN)

OBJ_TFF_MEDIAN = $(addsuffix .otffmedian, $(basename $(SRC_FLAG)))
TFF_MEDIAN: $(OBJ_NOFLAG) $(OBJ_TFF_MEDIAN)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_NOFLAG) $(OBJ_TFF_MEDIAN)

OBJ_TFF_MAX = $(addsuffix .otffmax, $(basename $(SRC_FLAG)))
TFF_MAX: $(OBJ_NOFLAG) $(OBJ_TFF_MAX)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_NOFLAG) $(OBJ_TFF_MAX)
# ---------------------------------------------------------------------------- #


# ------------------------------------ FTF ----------------------------------- #
OBJ_FTF_CE  = $(addsuffix .oftfce, $(basename $(SRC_FLAG)))
FTF_CE: $(OBJ_NOFLAG) $(OBJ_FTF_CE)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_NOFLAG) $(OBJ_FTF_CE)

OBJ_FTF_CG  = $(addsuffix .oftfcg, $(basename $(SRC_FLAG)))
FTF_CG: $(OBJ_NOFLAG) $(OBJ_FTF_CG)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_NOFLAG) $(OBJ_FTF_CG)
# ---------------------------------------------------------------------------- #


# ------------------------------------ FFT ----------------------------------- #
OBJ_FFT_S = $(addsuffix .offts, $(basename $(SRC_FLAG)))
FFT_S: $(OBJ_NOFLAG) $(OBJ_FFT_S)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_NOFLAG) $(OBJ_FFT_S)
# ---------------------------------------------------------------------------- #


# ------------------------------------ TTF ----------------------------------- #
OBJ_TTF_MEAN_CG = $(addsuffix .ottfmeancg, $(basename $(SRC_FLAG)))
TTF_MEAN_CG: $(OBJ_NOFLAG) $(OBJ_TTF_MEAN_CG)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_NOFLAG) $(OBJ_TTF_MEAN_CG)
# ---------------------------------------------------------------------------- #


# ------------------------------------ TFT ----------------------------------- #
OBJ_TFT_MEDIAN_S = $(addsuffix .otftmedians, $(basename $(SRC_FLAG)))
TFT_MEDIAN_S: $(OBJ_NOFLAG) $(OBJ_TFT_MEDIAN_S)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_NOFLAG) $(OBJ_TFT_MEDIAN_S)

OBJ_TFT_MEAN_S = $(addsuffix .otftmeans, $(basename $(SRC_FLAG)))
TFT_MEAN_S: $(OBJ_NOFLAG) $(OBJ_TFT_MEAN_S)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_NOFLAG) $(OBJ_TFT_MEAN_S)

# ------------------------------------ FTT ----------------------------------- #S
OBJ_FTT_CG_S = $(addsuffix .ofttcgs, $(basename $(SRC_FLAG)))
FTT_CG_S: $(OBJ_NOFLAG) $(OBJ_FTT_CG_S)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_NOFLAG) $(OBJ_FTT_CG_S)

# ------------------------------------ TTT ----------------------------------- #

# dsmga2-tos-ttt_moE
OBJ_TTT_MEAN_MOE_S = $(addsuffix .otttmeanmoes, $(basename $(SRC_FLAG)))
TTT_MEAN_MOE_S: $(OBJ_NOFLAG) $(OBJ_TTT_MEAN_MOE_S)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_NOFLAG) $(OBJ_TTT_MEAN_MOE_S)

# dsmga2-tos-ttt_moG
OBJ_TTT_MEAN_MOG_S = $(addsuffix .otttmeanmogs, $(basename $(SRC_FLAG)))
TTT_MEAN_MOG_S: $(OBJ_NOFLAG) $(OBJ_TTT_MEAN_MOG_S)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_NOFLAG) $(OBJ_TTT_MEAN_MOG_S)

# dsmga2-tos-ttt_cntE
OBJ_TTT_MEAN_CE_S = $(addsuffix .otttmeances, $(basename $(SRC_FLAG)))
TTT_MEAN_CE_S: $(OBJ_NOFLAG) $(OBJ_TTT_MEAN_CE_S)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_NOFLAG) $(OBJ_TTT_MEAN_CE_S)

# dsmga2-tos-ttt_cntG
OBJ_TTT_MEAN_CG_S = $(addsuffix .otttmeancgs, $(basename $(SRC_FLAG)))
TTT_MEAN_CG_S: $(OBJ_NOFLAG) $(OBJ_TTT_MEAN_CG_S)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(TLIB) -o $@ $(OBJ_NOFLAG) $(OBJ_TTT_MEAN_CG_S)

# -------------------------------- STEP 5 ENDS ------------------------------- #

#-----Other stuff----------------------------
depend:
	makedepend -Y. $(SRC)

clean:
	rm -f $(OBJ)

# DO NOT DELETE

chromosome.o: spin.h chromosome.h global.h myrand.h mt19937ar.h
chromosome.o: bitwisedistance.h nk-wa.h doublelinkedlistarray.h zkey.h sat.h maxcut.h traps.h
dsmga2.o: chromosome.h global.h myrand.h mt19937ar.h bitwisedistance.h
dsmga2.o: spin.h nk-wa.h doublelinkedlistarray.h zkey.h sat.h dsmga2.h maxcut.h traps.h
dsmga2.o: statistics.h trimatrix.h fastcounting.h
fastcounting.o: global.h myrand.h mt19937ar.h bitwisedistance.h spin.h
fastcounting.o: nk-wa.h doublelinkedlistarray.h zkey.h sat.h fastcounting.h maxcut.h traps.h
global.o: myrand.h mt19937ar.h statistics.h doublelinkedlistarray.h zkey.h
global.o: chromosome.h global.h bitwisedistance.h spin.h nk-wa.h sat.h maxcut.h traps.h
main.o: statistics.h dsmga2.h chromosome.h global.h myrand.h mt19937ar.h
main.o: bitwisedistance.h spin.h nk-wa.h doublelinkedlistarray.h zkey.h
main.o: sat.h trimatrix.h fastcounting.h maxcut.h traps.h
myrand.o: myrand.h mt19937ar.h
spin.o: global.h myrand.h mt19937ar.h bitwisedistance.h spin.h
spin.o: nk-wa.h doublelinkedlistarray.h zkey.h sat.h maxcut.h traps.h
nk-wa.o: nk-wa.h
sat.o: sat.h 
maxcut.o: maxcut.h
traps.o: traps.h
chromosome.o: spin.h chromosome.h global.h myrand.h mt19937ar.h
chromosome.o: bitwisedistance.h nk-wa.h doublelinkedlistarray.h zkey.h sat.h maxcut.h traps.h
dsmga2.o: chromosome.h global.h myrand.h mt19937ar.h bitwisedistance.h
dsmga2.o: spin.h nk-wa.h doublelinkedlistarray.h zkey.h sat.h dsmga2.h maxcut.h traps.h
dsmga2.o: statistics.h trimatrix.h fastcounting.h
fastcounting.o: global.h myrand.h mt19937ar.h bitwisedistance.h spin.h
fastcounting.o: nk-wa.h doublelinkedlistarray.h zkey.h sat.h fastcounting.h maxcut.h traps.h
global.o: myrand.h mt19937ar.h statistics.h doublelinkedlistarray.h zkey.h
global.o: chromosome.h global.h bitwisedistance.h spin.h nk-wa.h sat.h maxcut.h traps.h
# sweep.o: statistics.h dsmga2.h chromosome.h global.h myrand.h mt19937ar.h
# sweep.o: bitwisedistance.h spin.h nk-wa.h doublelinkedlistarray.h
# sweep.o: zkey.h sat.h trimatrix.h fastcounting.h maxcut.h
myrand.o: myrand.h mt19937ar.h
spin.o: global.h myrand.h mt19937ar.h bitwisedistance.h spin.h
spin.o: nk-wa.h doublelinkedlistarray.h zkey.h sat.h maxcut.h traps.h
nk-wa.o: nk-wa.h
sat.o: sat.h
maxcut.o: maxcut.h
traps.o: traps.h
nono.o: nono.h
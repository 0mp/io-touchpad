CC=gcc
CFLAGS=-c -fPIC
LDFLAGS=-shared

APP_DIR=app
APP_CLASSIFIER_DATA_DIR=app/classifier/data
APP_CLASSIFIER_DATA_USERDEFINED=app/classifier/data/user-defined
APP_CLASSIFIER_DATA_32=app/classifier/data/32
APP_CLASSIFIER_DATA_64=app/classifier/data/64
APP_TOOLS_DATA=app/tools/data
APP_TOOLS_DATA_MATRIXANALYSER=app/tools/data/matrixanalyser
APP_TOOLS_DATA_MATRIXANALYSER_FIGURES_DIR=app/tools/data/matrixanalyser/figures
PYCACHE_DIR=__pycache__

LIB_DIR=lib
OBJ_DIR=obj
SRC_DIR=src

DATA_DIRS=$(APP_CLASSIFIER_DATA_USERDEFINED) \
		  $(APP_TOOLS_DATA) \
		  $(APP_TOOLS_DATA_MATRIXANALYSER) \
		  $(APP_TOOLS_DATA_MATRIXANALYSER_FIGURES_DIR)

PERSISTANT_DATA_DIRS=$(APP_CLASSIFIER_DATA_DIR) \
					 $(APP_CLASSIFIER_DATA_32) \
					 $(APP_CLASSIFIER_DATA_64)

SILENCE_ERROR_MESSAGES=2>/dev/null || true

.PHONY: all check clean cleanpychache directories touchpadlib

all: directories touchpadlib


directories: $(DATA_DIRS) $(PERSISTANT_DATA_DIRS)

$(DATA_DIRS):
	-@mkdir $@ 2>/dev/null || true

$(PERSISTANT_DATA_DIRS):
	-@mkdir $@ 2>/dev/null || true

touchpadlib: $(LIB_DIR)/touchpadlib.so

$(LIB_DIR)/touchpadlib.so: $(OBJ_DIR)/touchpadlib.o
	-@mkdir $(LIB_DIR) 2>/dev/null || true
	$(CC) $(LDFLAGS) -Wl,-soname,touchpadlib.so -o $(LIB_DIR)/touchpadlib.so $(OBJ_DIR)/touchpadlib.o

$(OBJ_DIR)/touchpadlib.o: $(SRC_DIR)/touchpadlib.c $(SRC_DIR)/touchpadlib.h
	-@mkdir $(OBJ_DIR) 2>/dev/null || true
	$(CC) $(CFLAGS) -o $(OBJ_DIR)/touchpadlib.o $(SRC_DIR)/touchpadlib.c


clean: cleanpychache clean_classifier_data_userdefined
	-@rm $(OBJ_DIR)/*.o $(LIB_DIR)/touchpadlib.so $(SILENCE_ERROR_MESSAGES)
	-@rmdir $(LIB_DIR) $(SILENCE_ERROR_MESSAGES)
	-@rmdir $(OBJ_DIR) $(SILENCE_ERROR_MESSAGES)
	-@rm -rf $(DATA_DIRS) $(SILENCE_ERROR_MESSAGES)

cleanpychache: $(wildcard $(APP_DIR)/*/$(PYCACHE_DIR))
	@rm -rf $^ || true

clean_classifier_data_userdefined: $(APP_CLASSIFIER_DATA_USERDEFINED)
	@rm -f $(APP_CLASSIFIER_DATA_USERDEFINED)/* || true



check:
	@pylint $(APP_DIR)/*
	@pep8 $(APP_DIR)/*
	@pep257 $(APP_DIR)/*

#---------------------------------------------------------------------#
# Makefile for mgos_prog1.cpp
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
# 0 Prepare the environment for making the executable
#
#    0.1 Set directory information
#---------------------------------------------------------------------#

MY_INC = ./src/include
MY_LIB = ./src/lib

#CGAL_INCDIR = /home/wiktor/CGAL-5.1/include
#CGAL_LIBDIR = 



#---------------------------------------------------------------------#
#    0.2 Set compiler flags
#---------------------------------------------------------------------#

CXX       = g++ -std=c++14 
CXXFLAGS  = -Wall -g -m64 -w
INCLUDE_DIR = -I$(MY_INC) #-I$(CGAL_INCDIR)

#---------------------------------------------------------------------#
#    0.3 Set linker flags
#---------------------------------------------------------------------#

LDFLAGS = \
	-Xlinker --start-group \
		-L$(MY_LIB) \
		#-L$(CGAL_INCDIR) \
		#-lgmp -lCGAL \
		#-L$(BULL_LIBDIR) -lbull \
		#-L$(BOOST_INCDIR) \
	-Xlinker --end-group
	
#---------------------------------------------------------------------#
#    0.4. Define source files 
#---------------------------------------------------------------------#

APP_SRC  = main.cpp
TARGET   = main
POINT_3D = ./src/lib/Point_3d.cpp
NODE_VV = ./src/lib/Node_VV.cpp
NODE_QTV = ./src/lib/Node_QTV.cpp
NODE_QTC = ./src/lib/Node_QTC.cpp
PROGRAM_INSTANCE = ./src/lib/Program_instance.cpp
UTILS = ./src/lib/Utils.cpp
VORONOID_DIAGRAM = ./src/lib/Voronoi_diagram.cpp
QUASI_TRIANGULATION = ./src/lib/Quasi_triangulation.cpp
#NODE_QTC = ./my_lib/Node_QTC.cpp
#FUNCTIONS = ./my_lib/Functions.cpp
#GATE = ./my_lib/Gate.cpp
#ALPHA_FACE = ./my_lib/Alpha_Face.cpp
#FACE = ./my_lib/Face.cpp

#---------------------------------------------------------------------#
# 1. Create  
#---------------------------------------------------------------------#

all: $(TARGET) 

$(TARGET): $(TARGET).o $(POINT_3D).o $(PROGRAM_INSTANCE).o $(NODE_VV).o $(NODE_QTV).o $(NODE_QTC).o $(UTILS).o $(VORONOID_DIAGRAM).o $(QUASI_TRIANGULATION).o
	${CXX} $(CXXFLAGS) -static $(TARGET).o $(POINT_3D) $(PROGRAM_INSTANCE) $(NODE_VV) $(NODE_QTV) $(NODE_QTC) $(UTILS) $(VORONOID_DIAGRAM) $(QUASI_TRIANGULATION) $(LDFLAGS) -o $(TARGET)

$(TARGET).o: $(APP_SRC) $(POINT_3D) $(PROGRAM_INSTANCE) $(NODE_VV) $(NODE_QTV) $(NODE_QTC) $(UTILS) $(VORONOID_DIAGRAM) $(QUASI_TRIANGULATION) $(GLOBALS)
	${CXX} -c $(CXXFLAGS) $(APP_SRC) $(INCLUDE_DIR) -o $(TARGET).o

	

$(POINT_3D).o: $(POINT_3D)
	${CXX} -c $(CXXFLAGS) $(POINT_3D) $(INCLUDE_DIR) -o $(POINT_3D).o

$(PROGRAM_INSTANCE).o: $(PROGRAM_INSTANCE)
	${CXX} -c $(CXXFLAGS) $(PROGRAM_INSTANCE) $(INCLUDE_DIR) -o $(PROGRAM_INSTANCE).o

$(NODE_VV).o: $(NODE_VV)
	${CXX} -c $(CXXFLAGS) $(NODE_VV) $(INCLUDE_DIR) -o $(NODE_VV).o

$(NODE_QTV).o: $(NODE_QTV)
	${CXX} -c $(CXXFLAGS) $(NODE_QTV) $(INCLUDE_DIR) -o $(NODE_QTV).o

$(NODE_QTC).o: $(NODE_QTC)
	${CXX} -c $(CXXFLAGS) $(NODE_QTC) $(INCLUDE_DIR) -o $(NODE_QTC).o

$(UTILS).o: $(UTILS)
	${CXX} -c $(CXXFLAGS) $(UTILS) $(INCLUDE_DIR) -o $(UTILS).o

$(VORONOID_DIAGRAM).o: $(VORONOID_DIAGRAM)
	${CXX} -c $(CXXFLAGS) $(VORONOID_DIAGRAM) $(INCLUDE_DIR) -o $(VORONOID_DIAGRAM).o

$(QUASI_TRIANGULATION).o: $(QUASI_TRIANGULATION)
	${CXX} -c $(CXXFLAGS) $(QUASI_TRIANGULATION) $(INCLUDE_DIR) -o $(QUASI_TRIANGULATION).o


#---------------------------------------------------------------------#
# 2. Clean 
#---------------------------------------------------------------------#
clean:
	rm -f $(TARGET).o $(TARGET) $(POINT_3D).o $(PROGRAM_INSTANCE).o $(NODE_VV).o $(NODE_QTV).o $(NODE_QTC).o $(UTILS).o $(VORONOID_DIAGRAM).o $(QUASI_TRIANGULATION).o


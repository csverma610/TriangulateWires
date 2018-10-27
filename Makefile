OBJS = main.o TriangulateWires.o ./ConvexHull/QuickHull.o
#CXX = clang++-6.0
CXX  = g++

CPPFLAGS = -O3 -fPIC -fopenmp -DBOOST_TIMER -DNDEBUG -fPIC
CPPFLAGS += -I$(QGLVIEWER_DIR)/include
CPPFLAGS += -I$(QTDIR)/include -I$(QTDIR)/include/QtCore -I$(QTDIR)/include/QtXml -I$(QTDIR)/include/QtWidgets -I$(QTDIR)/include/QtGui -I$(QTDIR)/include/QtOpenGL
CPPFLAGS += -I$(BOOST_DIR)/include

LIBS     += -L$(QGLVIEWER_DIR)/lib -lQGLViewer
LIBS     += -lGL -lGLU
LIBS     += -L$(QTDIR)/lib -lQt5Widgets -lQt5Core -lQt5Gui -lQt5OpenGL
LIBS     += -L$(BOOST_DIR)/lib -lboost_timer
LIBS     += -lgomp

all = twires lview

twires:$(OBJS)
	$(CXX) -o twires $(OBJS) $(LIBS)

lview:LatticeViewer.o 
	$(CXX) -o lview LatticeViewer.o $(LIBS)

.o:.cpp
	$(CXX) $(CPPFLAGS) $<

clean:
	\rm -rf *.o twires lview

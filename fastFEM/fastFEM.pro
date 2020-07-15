#-------------------------------------------------
#
# Project created by QtCreator 2017-03-27T09:40:41
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport#necessary for qcustomplot

#output excution file name
TARGET = fastFEM

#project type:app lib vcapp vclib
TEMPLATE = app

QMAKE_CXXFLAGS += /wd"4819"

DEFINES += _CRT_SECURE_NO_WARNINGS

DESTDIR = $$PWD/../bin


!debug_and_release|build_pass {
    CONFIG(debug, debug|release) {
        TARGET = $$member(TARGET, 0)d
    }
}

unix {

CONFIG(debug,debug|release){
#output directory
#DESTDIR = ../x64/Linux/debug
}else{
#output directory
#DESTDIR = ../x64/Linux/release
}
#librarys
LIBS += \
    -L../SuperLU_MT_3.1 -lsuperlu_mt_OPENMP \
    -lgomp\
    -lpthread \
    ../armadillo/lib/libarmadillo.so \
    ../openblas/lib/libopenblas.so

#additional include path
INCLUDEPATH += \
    ../armadillo/include \
    ../qcustomplot \
    ../SuperLU_MT_3.1/SRC \
}

win32 {
CONFIG(debug,debug|release){
#output directory
#DESTDIR = ../x64/Debug
}else{
#output directory
#DESTDIR = ../x64/Release
}
#librarys
LIBS += \
    ..\SuperLU_MT_3.1\SuperLU_MT_3.1.lib \
    ..\OpenBLAS-v0.2.15-Win64-int32\lib\libopenblas.dll.a

#additional include path
INCLUDEPATH += \
    ..\armadillo\include \
    ..\qcustomplot \
    ..\SuperLU_MT_3.1\SRC \
        ..\OpenBLAS-v0.2.15-Win64-int32\include \

}


include($$PWD/../gmsh/gmsh.pri)

#include($$PWD/../qcustomplot/qcustomplot.pri)

QMAKE_CXXFLAGS += /openmp

CONFIG += debug_and_release


#source file lists
SOURCES += \
    coil3dtest.cpp \
    datatype.cpp \
    libstructure.cpp \
    mag2dtime.cpp \
    mag3dstatic.cpp \
    mag3dtime.cpp \
    main.cpp \
    material.cpp \
    SuperLU_MT.cpp \
    mesh/meshtype.cpp \
    relay.cpp \
    relay1250.cpp \
    testremesh.cpp \
    tetralheral.cpp \
    to5relay.cpp \
    to5time.cpp \
    valverelay.cpp \
    valvetime.cpp


#.ui file lists
FORMS    +=

#header file lists
HEADERS += \
    coil3dtest.h \
    datatype.h \
    libstructure.h \
    mag2dtime.h \
    mag3dstatic.h \
    mag3dtime.h \
    relay.h \
    relay1250.h \
    spline.h \
    SuperLU_MT.h \
    mesh/meshtype.h \
    testremesh.h \
    tetralheral.h \
    to5relay.h \
    to5time.h \
    valverelay.h \
    valvetime.h


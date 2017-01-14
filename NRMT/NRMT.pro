#-------------------------------------------------
#
# Project created by QtCreator 2016-06-26T14:16:48
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = NRMT
TEMPLATE = app

SOURCES += main.cpp\
        mainwindow.cpp \
    datatype.cpp\
    qcustomplot.cpp

HEADERS  += mainwindow.h \
    datatype.h\
    qcustomplot.h

FORMS    += mainwindow.ui
INCLUDEPATH += .\
    E:\Projects\cplusplus\SuperLU_MT_3.1\SRC \
	E:\Projects\cplusplus\armadillo-6.500.4\include \
	E:\Projects\cplusplus\OpenBLAS-v0.2.15-Win64-int32\include
	

LIBS += SuperLU_MT_3.1.lib \
    E:\Projects\cplusplus\OpenBLAS-v0.2.15-Win64-int32\lib\libopenblas.dll.a \
    

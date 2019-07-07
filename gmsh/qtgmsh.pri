
win32:CONFIG(release, debug|release): LIBS += LIBS += -L$$PWD/../bin/ -lgmsh
     else:LIBS += LIBS += -L$$PWD/../bin/ -lgmsh

INCLUDEPATH += $$PWD/
DEPENDPATH += $$PWD/



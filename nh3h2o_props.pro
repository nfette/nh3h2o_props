TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    aqua.cpp

include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    aqua.h


include(../defaults.pri)

CONFIG   -= app_bundle
CONFIG += c++11

LIBS += -L$$OUT_PWD/../Model -lmodel
LIBS += -larmadillo
LIBS += -L/usr/local/lib/RpolyPlusPlus -lrpoly_plus_plus

TEMPLATE = lib

TARGET = wellindexcalculator

HEADERS += \
    geometry_functions/geometryfunctions.h \
    geometry_functions/geometryfunctions_exceptions.h \
    geometry_classes/line.h \
    geometry_classes/squareplane.h \
    well_constraint_projections/well_constraint_projections.h

SOURCES += \
    geometry_functions/geometryfunctions.cpp \
    geometry_classes/line.cpp \
    geometry_classes/squareplane.cpp \
    well_constraint_projections/well_constraint_projections.cpp

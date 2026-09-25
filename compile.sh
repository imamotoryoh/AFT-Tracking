g++ -Iinclude \
    src/AFTGeometryDisplay.cc src/AFTGeometry.cc \
    $(root-config --cflags --libs) \
    -o bin/AFTGeometryDisplay

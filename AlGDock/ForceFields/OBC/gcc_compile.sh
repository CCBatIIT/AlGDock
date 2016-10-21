
g++ -c ObcParameters.cpp -o ObcParameters.o
g++ -c ReferenceForce.cpp -o ReferenceForce.o
g++ -c ReferenceObc.cpp -o ReferenceObc.o
g++ -c ObcWrapper.cpp -o ObcWrapper.o

# c++  
g++ -c test.cpp -o test_cpp.o
g++ test_cpp.o ObcParameters.o ReferenceForce.o ReferenceObc.o -o test_cpp

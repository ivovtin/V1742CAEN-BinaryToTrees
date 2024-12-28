# V1742CAEN-BinaryToTrees

To get the code, you need to run the command: <br />
```
git clone https://github.com/ivovtin/V1742CAEN-BinaryToTrees.git
```

To assemble use:
```
gcc -o procn24dec24 -lstdc++ `root-config --libs --glibs --cflags` procn24dec24.cc
```

For launch convertation: <br />
```
procn24dec24 7 wave_0.dat wave_2.dat wave_4.dat wave_6.dat wave_9.dat wave_11.dat wave_13.dat test.root -100 2
```

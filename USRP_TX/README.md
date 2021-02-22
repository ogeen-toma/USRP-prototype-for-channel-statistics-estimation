1) Once this folder is downloaded open terminal window on this directory, which should shows something like .../USRP_TX$.
2) Type the following commands to build the transmitter program:

.../USRP_TX$ mkdir build

.../USRP_TX$ cd build

.../USRP_TX/build$ cmake ../

.../USRP_TX/build$ make

3) Finally, run the transmitter program either from GUI or from terminal window:

a) GUI:

b) Terminal window

5) Then copy config.txt file into build folder, also can be done from terminal window as:

.../USRP_TX/build$ cp ../config.txt ../build/

Note that config.txt file is used by the source code USRP_TX.c to configure the USRP as a PU, which generates a traffic activity in a frequency channel.
This configuration can be modified from config.txt file, where detailed description for each parameter is provided.

4) Finally use the following command to run the transmitter program:

.../USRP_TX/build$ sudo chrt --rr 99 ./USRP_TX config.txt

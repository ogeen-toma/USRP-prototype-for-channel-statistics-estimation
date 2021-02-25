1) Once USRP_TX folder is downloaded open the terminal window on this directory, which should show something like .../USRP_TX$.
2) Type the following commands to build the transmitter program:

        .../USRP_TX$ mkdir build

        .../USRP_TX$ cd build

        .../USRP_TX/build$ cmake ../

        .../USRP_TX/build$ make

3) Finally, run the transmitter program either from GUI or terminal window:

     a) GUI: simply open the executable file TX in the main directory and change the configurations therein as desired and then start transmitting.
     Note: if TX could not display, make sure it is given execution permission from its properties.
     
     b) Terminal window: transmitter program can also be run manually from terminal window as follows:
     
     First, copy config.txt file into build folder or by using the following command in terminal window:
     
        .../USRP_TX/build$ cp ../config.txt ../build/
     
     Note: config.txt file is used by the source code USRP_TX.c to configure the USRP as a PU, which generates a traffic activity in a frequency channel. The settings of the PU can be configured manually from config.txt file, where detailed description for each parameter is provided.

     Finally, use the following command to run the transmitter program:

        .../USRP_TX/build$ sudo chrt --rr 99 ./USRP_TX config.txt



1) Once USRP_RX folder is downloaded open the terminal window on this directory, which should show something like .../USRP_RX$.
2) Type the following commands to build the receiver program:

        .../USRP_RX$ mkdir build

        .../USRP_RX$ cd build

        .../USRP_RX/build$ cmake ../

        .../USRP_RX/build$ make

3) Finally, run the receiver program either from GUI or terminal window:

     a) GUI: simply open the executable file RX in the main directory and change the configurations therein as desired and then start receiving.
     Note: if RX could not display, make sure it is given execution permission from its properties.
     
     b) Terminal window: receiver program can also be run manually from terminal window as follows:
     
     First, copy config.txt file into build folder or by using the following command in terminal window:
     
        .../USRP_RX/build$ cp ../config.txt ../build/
     
     Note: config.txt file is used by the source code USRP_RX.c to configure the USRP as a SU, which performs spectrum sensing in a frequency channel. This configuration can be modified manually from config.txt file, where detailed description for each parameter is provided.

     Finally, use the following command to run the receiver program:

        .../USRP_RX/build$ sudo chrt --rr 99 ./USRP_RX config.txt




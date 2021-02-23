1) Once this folder is downloaded open terminal window on this directory, which should show something like .../USRP_TX$.
2) Type the following commands to build the transmitter program:

.../USRP_TX$ mkdir build

.../USRP_TX$ cd build

.../USRP_TX/build$ cmake ../

.../USRP_TX/build$ make

3) Finally, run the transmitter program either from GUI or terminal window:

     a) GUI: simply open the executable file TX in the main directory and change the configurations therein as desired and then start transmitting.
     (Note: make sure the TX file is given the permission of execution from its properties)

     b) Terminal window: to run the the transmitter program manually from terminal window. First, copy config.txt file into build folder or by using this command:
     
     .../USRP_TX/build$ cp ../config.txt ../build/

     Second, use the following command to run the transmitter program:

     .../USRP_TX/build$ sudo chrt --rr 99 ./USRP_TX config.txt


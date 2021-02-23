
/***********************************************************************
 *
 * USRP-Based Prototype for Real-Time Estimation of Channel Activity Statistics 
 *				in Smart Spectrum Sharing
 *
 ***********************************************************************
 *
 *                    		RECEIVER PROGRAM
 *
 ***********************************************************************
 *
 * Mr Ogeen H. Toma
 * Dr Miguel López-Benítez
 * Department of Electrical Engineering and Electronics
 * University of Liverpool, United Kingdom
 *
 * Email: Ogeen.toma@liverpool.ac.uk
 * Email: M.Lopez-Benitez@liverpool.ac.uk
 *
 ***********************************************************************
 *
 * This program was developed and tested on USRP B200mini connected to 
 * a host computer with Linux distribution Ubuntu 18.04.3 LTS.
 *
 ***********************************************************************
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ***********************************************************************
 *
 * Version history:
 *
 * 		v1.00 	01 Feb 2021
 *
 **********************************************************************/

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <termios.h>
#include <time.h>
#include <unistd.h>

#include <uhd.h>

#define EXECUTE_OR_GOTO(label, ...) \
if(__VA_ARGS__){ \
*return_code = EXIT_FAILURE; \
goto label; \
}

// Current version
#define VERSION 1.00

// Number of points in a plot
#define NOF_PLOT_POINTS 300

// Number of points to calculate CDF and percentiles
#define NOF_CDF_POINTS 1000

// Maximum number of characters per line in the configuration file
#define MAX_LINE_LENGTH_IN_CONFIG_FILE	100

// Structure with all configuration parameters
struct config_struct
{
    // Program operation parameters
    double ReceptionTime;
    double Ts;
    double tau;
    double DecisionThreshold;
    int TestMode;
    int DisplayEnergyStats;
    int SaveEnergyValues;
    int ShowDetectedPeriods;

    // USRP parameters
    double CenterFrequency;
    double Bandwidth;
    double SampleRate;
    double TunerGain;
    char USRPdevice[16];

    // Other calculated parameters
    int N;					// This is the exact number of (complex I/Q) samples per sensing event required
    // to meet the desired duration of tau at the specified sampling rate. This is
    // the number of (complex I/Q) samples that will be used for sensing decisions.
};

// Structure with channel sensing results
struct channel_sensing_results_struct
{
    double *EnergyValues;
    double *detected_IdlePeriods;
    double *detected_BusyPeriods;
    unsigned long long int nof_EnergyValues;
    unsigned long long int nof_detected_Periods;
};

// Function prototypes
void display_header(void);
void set_default_config_params(struct config_struct *config);
void load_config_file(int argc, char *argv[], struct config_struct *config);
void check_config_params(struct config_struct config);
void show_config_params(struct config_struct config);
void open_plotter(FILE **gnuPlotter);
void close_plotter(FILE *gnuPlotter);
void sense_channel(struct config_struct config, struct channel_sensing_results_struct *channel_sensing_results, FILE *gnuPlotter, int *return_code);
void display_energy_statistics(struct channel_sensing_results_struct channel_sensing_results);
void save_energy_values(struct channel_sensing_results_struct channel_sensing_results);
void save_detected_periods(struct channel_sensing_results_struct channel_sensing_results);
int continue_program(void);
double elapsed_time_ms(struct timespec start_time, struct timespec end_time);
void plot_energy_history(FILE *gnuPlotter, int DisplayEnergyStats, double *x, double *y, int nof_points);
void plot_mean_history(FILE *gnuPlotter_mean, double *x, double *y, int nof_points);
void plot_DC_history(FILE *gnuPlotter_DC, double *x, double *y, int nof_points);
void plot_distribution_history(FILE *gnuPlotter_distribution, double *y, double ts, int nof_points);
void calculate_cdf(double *vector, unsigned long long int vector_size, double *cdf_x_vector, double *cdf_y_vector, int cdf_size);
double calculate_percentile(double *cdf_x_vector, double *cdf_y_vector, int cdf_size, double percentile);


int main(int argc, char *argv[])
{
    // Configuration structure, channel sensing results, device handle and plotter
    struct config_struct config;
    struct channel_sensing_results_struct channel_sensing_results;
    FILE *gnuPlotter;
    //FILE *gnuPlotter_mean, *gnuPlotter_DC, *gnuPlotter_distribution;

    // Display header
    display_header();

    // Set default configuration parameters
    set_default_config_params(&config);

    // Load configuration parameters from the config file
    load_config_file(argc, argv, &config);

    // Check configuration parameters
    check_config_params(config);

    // Show configuration parameters
    show_config_params(config);

    // Open plotter
    if(config.TestMode)
        open_plotter(&gnuPlotter);

//    if(config.DisplayEnergyStats)
//    {
//        open_plotter(&gnuPlotter_mean);
//        open_plotter(&gnuPlotter_DC);
//        open_plotter(&gnuPlotter_distribution);
//    }

    // Final confirmation before starting experiment
    printf("\n\nReady to sense the channel for %.0f seconds (%.2d:%.2d:%.2d)\n", config.ReceptionTime, (int) floor(config.ReceptionTime/3600.0), (int) floor(((unsigned long int) config.ReceptionTime % 3600)/60.0), (int) (((unsigned long int) config.ReceptionTime % 3600) % 60));
    printf("\n  Sensing period (Ts)         : %g ms", config.Ts*1000);
    printf("\n  Sensing time (tau)          : %g ms", config.tau*1000);
    printf("\n  Decision threshold          : %g", config.DecisionThreshold);
    printf("\n  Samples per sensing event   : %d", config.N);

    int return_code = EXIT_SUCCESS;
    //if(continue_program())
	if(1)
    {
        // Sense channel
        sense_channel(config, &channel_sensing_results, gnuPlotter, &return_code);

        // Close plotter
        if(config.TestMode)
            close_plotter(gnuPlotter);

        // Display energy statistics
//        if(config.DisplayEnergyStats){
//            display_energy_statistics(channel_sensing_results);
//            close_plotter(gnuPlotter_mean);
//            close_plotter(gnuPlotter_DC);
//            close_plotter(gnuPlotter_distribution);
//        }

        // Save energy values
        if(config.SaveEnergyValues)
            save_energy_values(channel_sensing_results);

        // Save detected periods
        save_detected_periods(channel_sensing_results);

        if (return_code == EXIT_FAILURE){
            printf("\nProgram execution failed.\n\n");
            return 0;
        }

        printf("\nProgram execution finished successfully\n\n");
    }
    else
    {

        // Close plotter
        if(config.TestMode)
            close_plotter(gnuPlotter);

//        if(config.DisplayEnergyStats)
//        {
//            close_plotter(gnuPlotter_mean);
//            close_plotter(gnuPlotter_DC);
//            close_plotter(gnuPlotter_distribution);
//        }

        printf("\nProgram execution finished by the user\n\n");
    }

    return 0;
}

void display_header(void)
{
    printf("\n########################################################");
    printf("\n#                                                      #");
    printf("\n#               USRP_RX (B200mini) v%.2f               #", VERSION);
    printf("\n#                                                      #");
    printf("\n########################################################");
    printf("\n#                                                      #");
    printf("\n# Mr Ogeen H. Toma                                     #");
    printf("\n# Dr Miguel López-Benítez                              #");
    printf("\n# Department of Electrical Engineering and Electronics #");
    printf("\n# University of Liverpool, United Kingdom              #");
    printf("\n#                                                      #");
    printf("\n# Email : Ogeen.toma@liverpool.ac.uk                   #");
    printf("\n# Email : M.Lopez-Benitez@liverpool.ac.uk              #");
    printf("\n#                                                      #");
    printf("\n########################################################\n");
}

void set_default_config_params(struct config_struct *config)
{
    // Sets default values for the configuration parameters
    // (will be changed if other values are specified in the config file)

    // Program operation parameters
    config->ReceptionTime = 60.0;			// 60 seconds
    config->Ts = 0.01;						// 10 ms
    config->tau = 0.001;					// 1 ms
    config->DecisionThreshold = -1;			// Not set
    config->TestMode = 0;					// OFF
    config->DisplayEnergyStats = 0;			// OFF
    config->SaveEnergyValues = 0;			// OFF
    config->ShowDetectedPeriods = 0;		// OFF

    // USRP parameters
    config->CenterFrequency = 433920000;	// 433.92 MHz
    config->Bandwidth = 1.74e6;             // 1.74 MHz
    config->SampleRate = 1000000;			// 1 MHz
    config->TunerGain = 0;					// 0 dB
    strncpy(config->USRPdevice,"serial=319474A",16);// "type=b200", 10

    // Other calculated parameters
    config->N = (int) ceil(config->tau*((double) config->SampleRate));
//    config->N =2040;
}

void load_config_file(int argc, char *argv[], struct config_struct *config)
{
    FILE *config_file;
    char line[MAX_LINE_LENGTH_IN_CONFIG_FILE];
    char *parameter_name, *parameter_value;

    // Check that the program has the right number of arguments (just one, the name of the config file)
    if(argc != 2)
    {
        printf("\nError: Invalid arguments. Use: %s config_filename\n\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    // Try to open the config file
    if((config_file = fopen(argv[1], "r")) == NULL)
    {
        printf("\nError: Could not open file %s\n\n", argv[1]);
        exit(EXIT_FAILURE);
    }
    else
    {
        printf("\nLoading configuration parameters from file %s ... ", argv[1]);
    }

    // Load configuration parameters from the config file
    while(fgets(line, MAX_LINE_LENGTH_IN_CONFIG_FILE, config_file) != NULL)
    {
        // Process the line only if it isn't a blank line (doesn't start with '\n') and it isn't a comment (doesn't start with #)
        if(line[0] != '\n' && line[0] != '#')
        {
            parameter_name = strtok (line, " =");
            parameter_value = strtok (NULL, " =");

            // printf("\nDebug: %s = %s", parameter_name, parameter_value);

            if(strcasecmp(parameter_name, "ReceptionTime") == 0)
            {
                config->ReceptionTime = atof(parameter_value);
            }

            if(strcasecmp(parameter_name, "Ts") == 0)
            {
                config->Ts = atof(parameter_value);
            }

            if(strcasecmp(parameter_name, "tau") == 0)
            {
                config->tau = atof(parameter_value);
            }

            if(strcasecmp(parameter_name, "DecisionThreshold") == 0)
            {
                config->DecisionThreshold = atof(parameter_value);
            }

            if(strcasecmp(parameter_name, "TestMode") == 0)
            {
                config->TestMode = atoi(parameter_value);
            }

            if(strcasecmp(parameter_name, "DisplayEnergyStats") == 0)
            {
                config->DisplayEnergyStats = atoi(parameter_value);
            }

            if(strcasecmp(parameter_name, "SaveEnergyValues") == 0)
            {
                config->SaveEnergyValues = atoi(parameter_value);
            }

            if(strcasecmp(parameter_name, "ShowDetectedPeriods") == 0)
            {
                config->ShowDetectedPeriods = atoi(parameter_value);
            }

            if(strcasecmp(parameter_name, "CenterFrequency") == 0)
            {
                config->CenterFrequency = atof(parameter_value);
            }

            if(strcasecmp(parameter_name, "Bandwidth") == 0)
            {
                config->Bandwidth = atof(parameter_value);
            }

            if(strcasecmp(parameter_name, "SampleRate") == 0)
            {
                config->SampleRate = atof(parameter_value);
            }

            if(strcasecmp(parameter_name, "TunerGain") == 0)
            {
                config->TunerGain = atof(parameter_value);
            }

            if(strcasecmp(parameter_name, "USRPdevice") == 0)
            {
                strncpy(config->USRPdevice,"serial=",16);
                strcat(config->USRPdevice,parameter_value);
            }
        }
    }

    // Recalculate these parameters
    config->N = (int) ceil(config->tau*((double) config->SampleRate));
//    config->N =2040;

    // Close the config file
    if(fclose(config_file))
    {
        printf("\nError: Could not close file %s\n\n", argv[1]);
        exit(EXIT_FAILURE);
    }
    else
    {
        printf("Done\n");
    }
}

void check_config_params(struct config_struct config)
{
    printf("\nChecking configuration parameters ... ");

    if(config.ReceptionTime <= 0.0)
    {
        printf("\nError: The reception time should be > 0\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.Ts <= 0.0)
    {
        printf("\nError: The sensing period (Ts) should be > 0\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.tau <= 0.0)
    {
        printf("\nError: The sensing time (tau) should be > 0\n\n");
        exit(EXIT_FAILURE);
    }
    else if(config.tau >= config.Ts)
    {
        printf("\nError: The sensing time (tau) should be shorter than the sensing period (Ts)\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.DecisionThreshold == -1)
    {
        printf("\nWarning: The decision threshold has not been set\n");
    }
    else if(config.DecisionThreshold <= 0)
    {
        printf("\nError: The decision threshold should be > 0 or -1\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.TestMode != 0 &&
       config.TestMode != 1)
    {
        printf("\nError: TestMode must be either 0 or 1\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.DisplayEnergyStats != 0 &&
       config.DisplayEnergyStats != 1)
    {
        printf("\nError: DisplayEnergyStats must be either 0 or 1\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.SaveEnergyValues != 0 &&
       config.SaveEnergyValues != 1)
    {
        printf("\nError: SaveEnergyValues must be either 0 or 1\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.ShowDetectedPeriods != 0 &&
       config.ShowDetectedPeriods != 1)
    {
        printf("\nError: ShowDetectedPeriods must be either 0 or 1\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.CenterFrequency < 70e6)
    {
        printf("\nError: The center frequency is too low\n\n");
        exit(EXIT_FAILURE);
    }
    else if(config.CenterFrequency > 6e9)
    {
        printf("\nError: The center frequency is too high\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.Bandwidth <= 200000 || config.Bandwidth > 56000000)
    {
        printf("\nError: The analog frontend filter bandwidth should be within the interval [200 KHz, 56 MHz]\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.SampleRate < 200000 || config.SampleRate > 56000000)
    {
        printf("\nError: The sample rate should be within the interval [200 KHz, 56 MHz]\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.SampleRate > 56e6)
    {
        printf("\nWarning: The sample rate is > 56 MHz, some samples might be lost\n");
    }

    if(config.TunerGain < 0 || config.TunerGain > 76)
    {
        printf("\nError: The tuner gain value should be in the range [0,76]\n\n");
        exit(EXIT_FAILURE);
    }

//    if(strcasecmp(config.USRPdevice, "serial=319474A") != 0)
//    {
//        printf("\nError: Enter the correct serial number of the USRP\n\n");
//        exit(EXIT_FAILURE);
//    }

    printf("Done\n");
}

void show_config_params(struct config_struct config)
{
    // Total reception time in HH:MM:SS format
    double hours, minutes, seconds;
    hours = floor(config.ReceptionTime/3600.0);
    minutes = floor(((unsigned long int) config.ReceptionTime % 3600)/60.0);
    seconds = ((unsigned long int) config.ReceptionTime % 3600) % 60;

    printf("\nProgram operation parameters:\n");
    printf("\n  Total reception time        : %.2d:%.2d:%.2d (%.0f seconds)", (int) hours, (int) minutes, (int) seconds, config.ReceptionTime);
    printf("\n  Sensing period (Ts)         : %g ms", config.Ts*1000);
    printf("\n  Sensing time (tau)          : %g ms", config.tau*1000);
    printf("\n  Decision threshold          : %g", config.DecisionThreshold);
    printf("\n  Test mode                   : %s", config.TestMode ? "ON" : "OFF");
    printf("\n  Display energy statistics   : %s", config.DisplayEnergyStats ? "ON" : "OFF");
    printf("\n  Save energy values          : %s", config.SaveEnergyValues ? "ON" : "OFF");
    printf("\n  Show detected periods       : %s", config.ShowDetectedPeriods ? "ON" : "OFF");

    printf("\n\nUSRP parameters:\n");
    printf("\n  Center frequency            : %g MHz", ((double) config.CenterFrequency)/1e6);
    printf("\n  Bandwidth                   : %g MHz", ((double) config.Bandwidth)/1e6);
    printf("\n  Sample rate                 : %g MHz", ((double) config.SampleRate)/1e6);
    printf("\n  Tuner gain                  : %g dB", (double) config.TunerGain);
    printf("\n  USRP serial number          : %s", config.USRPdevice);

    printf("\n\nOther parameters:\n");
    printf("\n  Samples per sensing event   : %d", config.N);
}

void open_plotter(FILE **gnuPlotter)
{
    printf("\n\nOpening plotter ... ");

    // Open pipe to communicate with gnuplot (gnuplot should be installed in the system)
    if((*gnuPlotter = popen("gnuplot", "w")) == NULL)
    {
        printf("\nError: Failed to open plotter (try sudo apt-get install gnuplot)\n\n");
        exit(EXIT_FAILURE);
    }

    // Sets the window's title
    fprintf(*gnuPlotter, "set term qt title \"USRP_RX (B200mini)\"\n");
    // Sets the background colour to black
    fprintf(*gnuPlotter, "set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb\"black\" behind\n");
    // Sets the properties of the lines in the axes
    fprintf(*gnuPlotter, "set border linewidth 2 linecolor rgb \"grey\"\n");
    // Sets the properties of the numbers in the axes
    fprintf(*gnuPlotter, "set tics font \"Sans,11\" textcolor rgb \"grey\"\n");
    // Style for plotting data: lines
    fprintf(*gnuPlotter, "set style data lines\n");
    // Sets the properties of linestyle #1
    fprintf(*gnuPlotter, "set linetype 1 linewidth 1 linecolor rgb \"yellow\"\n");
    // Sets linestyle #1 as the default type
    fprintf(*gnuPlotter, "set style line 1 default\n");
    // Removes the key/legend (will not be shown)
    fprintf(*gnuPlotter, "unset key\n");

    fprintf(*gnuPlotter, "set terminal qt size 959, 650\n");

    fprintf(*gnuPlotter, "set multiplot\n");

    // Flush commands
    fflush(*gnuPlotter);

    printf("Done");
}

void close_plotter(FILE *gnuPlotter)
{
    printf("\nClosing plotter ... ");

    // Open pipe to communicate with gnuplot (gnuplot should be installed in the system)
    if(pclose(gnuPlotter) == -1)
    {
        printf("\nError: Failed to close plotter\n\n");
        exit(EXIT_FAILURE);
    }

    printf("Done\n");
}



void sense_channel(struct config_struct config, struct channel_sensing_results_struct *channel_sensing_results, FILE *gnuPlotter, int *return_code)
{
    int i, nof_read_samples;
    struct timespec ReceptionStartTime, LastSensingTime, PeriodStartTime, CurrentTime;
    double current_energy, current_period_duration, *EnergyValues, *detected_IdlePeriods, *detected_BusyPeriods;
    unsigned long long int nof_EnergyValues, nof_detected_Periods, max_nof_elements;
    int CurrentPeriodType; // 0 = idle, 1 = busy
    double max_elapsed_time_between_sensing_events_ms = -1.0;
    double gnuPlotter_x_vector[NOF_PLOT_POINTS], gnuPlotter_y_vector[NOF_PLOT_POINTS];

    double *gnuPlotter_x_vector_mean_idle, *gnuPlotter_y_vector_mean_idle;
    double *gnuPlotter_x_vector_mean_busy, *gnuPlotter_y_vector_mean_busy, *gnuPlotter_y_vector_DC;
    double *gnuPlotter_y_vector_distribution;
    //double gnuPlotter_y_vector_distribution[1000]={0};

    /////////////////////////////////////////////////////////////USRP configuration////////////////////////////////////////////////////////////////////////////
    size_t channel = 0;
    char error_string[512];
//    size_t n_samples=config.N;//=2040

    // Create USRP
    uhd_usrp_handle usrp;
    fprintf(stderr, "\n\nConnecting to USRP \"%s\"...\n", config.USRPdevice);
    EXECUTE_OR_GOTO(free_option_strings,
                    uhd_usrp_make(&usrp, config.USRPdevice)
                    )

    // Create RX streamer
    uhd_rx_streamer_handle rx_streamer;
    EXECUTE_OR_GOTO(free_usrp,
                    uhd_rx_streamer_make(&rx_streamer)
                    )

    // Create RX metadata
    uhd_rx_metadata_handle md;
    EXECUTE_OR_GOTO(free_rx_streamer,
                    uhd_rx_metadata_make(&md)
                    )

    // Create other necessary structs
    uhd_tune_request_t tune_request = {
        .target_freq = config.CenterFrequency,
        .rf_freq_policy = UHD_TUNE_REQUEST_POLICY_AUTO,
        .dsp_freq_policy = UHD_TUNE_REQUEST_POLICY_AUTO,
    };
    uhd_tune_result_t tune_result;

    uhd_stream_args_t stream_args = {
        .cpu_format = "fc32",
        .otw_format = "sc16",
        .args = "",
        .channel_list = &channel,
        .n_channels = 1
    };

    uhd_stream_cmd_t stream_cmd = {
        .stream_mode = UHD_STREAM_MODE_NUM_SAMPS_AND_DONE,
        .num_samps = config.N,
        .stream_now = true
    };

    size_t samps_per_buff;
    float *buff = NULL;
    void **buffs_ptr = NULL;
//    FILE *fp = NULL;
    size_t num_acc_samps = 0;

    // Set rate
    fprintf(stderr, "Setting RX Rate: %f...\n", config.SampleRate);
    EXECUTE_OR_GOTO(free_rx_metadata,
                    uhd_usrp_set_rx_rate(usrp, config.SampleRate, channel)
                    )

    // See what rate actually is
    EXECUTE_OR_GOTO(free_rx_metadata,
                    uhd_usrp_get_rx_rate(usrp, channel, &config.SampleRate)
                    )
    fprintf(stderr, "Actual RX Rate: %f...\n", config.SampleRate);

    // Set gain
    fprintf(stderr, "Setting RX Gain: %f dB...\n", config.TunerGain);
    EXECUTE_OR_GOTO(free_rx_metadata,
                    uhd_usrp_set_rx_gain(usrp, config.TunerGain, channel, "")
                    )

    // See what gain actually is
    EXECUTE_OR_GOTO(free_rx_metadata,
                    uhd_usrp_get_rx_gain(usrp, channel, "", &config.TunerGain)
                    )
    fprintf(stderr, "Actual RX Gain: %f...\n", config.TunerGain);

    // Set frequency
    fprintf(stderr, "Setting RX frequency: %f MHz...\n", config.CenterFrequency/1e6);
    EXECUTE_OR_GOTO(free_rx_metadata,
                    uhd_usrp_set_rx_freq(usrp, &tune_request, channel, &tune_result)
                    )

    // See what frequency actually is
    EXECUTE_OR_GOTO(free_rx_metadata,
                    uhd_usrp_get_rx_freq(usrp, channel, &config.CenterFrequency)
                    )
    fprintf(stderr, "Actual RX frequency: %f MHz...\n", config.CenterFrequency / 1e6);

    // Set up streamer
    stream_args.channel_list = &channel;
    EXECUTE_OR_GOTO(free_rx_streamer,
                    uhd_usrp_get_rx_stream(usrp, &stream_args, rx_streamer)
                    )

    // Set up buffer
    EXECUTE_OR_GOTO(free_rx_streamer,
                    uhd_rx_streamer_max_num_samps(rx_streamer, &samps_per_buff)
                    )
    fprintf(stderr, "Buffer size in samples: %zu\n", samps_per_buff);
    buff = malloc(samps_per_buff * 2 * sizeof(float));
    buffs_ptr = (void**)&buff;

    uhd_usrp_set_rx_bandwidth(usrp, config.Bandwidth ,channel);

//    EXECUTE_OR_GOTO(free_rx_streamer,
//                    uhd_usrp_set_rx_agc(usrp, false, channel)
//                    )

    // Issue stream command
    //    fprintf(stderr, "Issuing stream command.\n");
    //    EXECUTE_OR_GOTO(free_buffer,
    //        uhd_rx_streamer_issue_stream_cmd(rx_streamer, &stream_cmd)
    //    )
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Allocate memory to the vectors that will store the energy values and the detected periods
    // One energy value per sensing event (plus 1000 extra)
    max_nof_elements = ((unsigned long long int) ceil(config.ReceptionTime/config.Ts)) + 1000;
    EnergyValues = (double*) malloc(max_nof_elements * sizeof(double));
    // One detected period per sensing event maximum, so this size should be enough
    detected_IdlePeriods = (double*) malloc(max_nof_elements * sizeof(double));
    detected_BusyPeriods = (double*) malloc(max_nof_elements * sizeof(double));
    nof_EnergyValues = 0;
    nof_detected_Periods = 0;

    gnuPlotter_x_vector_mean_idle=(double*) malloc(max_nof_elements * sizeof(double));
    gnuPlotter_y_vector_mean_idle=(double*) malloc(max_nof_elements * sizeof(double));
    gnuPlotter_x_vector_mean_busy=(double*) malloc(max_nof_elements * sizeof(double));
    gnuPlotter_y_vector_mean_busy=(double*) malloc(max_nof_elements * sizeof(double));
    gnuPlotter_y_vector_DC=(double*) malloc(max_nof_elements * sizeof(double));
    gnuPlotter_y_vector_distribution=(double*) malloc(max_nof_elements * sizeof(double));


    // Reset period type register (0 = idle / 1 = busy). Initially set to an unknown value.
    CurrentPeriodType = -1;

    // Reset plot vector
    if(config.TestMode)
    {
        for(i=0; i<NOF_PLOT_POINTS; i++)
        {
            gnuPlotter_x_vector[i] = i + 1.0;
            gnuPlotter_y_vector[i] = 0.0;
        }
    }


    printf("\n***** The transmitter can now be started *****\n");

    // Start sensing the channel
    printf("\nSensing center frequency %g MHz at %g MSamples/sec ...\n", config.CenterFrequency/1e6, config.SampleRate/1e6);
    fflush(stdout);

    clock_gettime(CLOCK_MONOTONIC, &ReceptionStartTime);
    LastSensingTime = ReceptionStartTime;
    PeriodStartTime = ReceptionStartTime;
    CurrentTime = ReceptionStartTime;

    // Start the reception
    while(elapsed_time_ms(ReceptionStartTime, CurrentTime) < config.ReceptionTime*1000.0)
    {
        // Check if a sensing event is required
        if(elapsed_time_ms(LastSensingTime, CurrentTime) >= config.Ts*1000.0)
        {
            // Debug (time error should be less than 1 microsecond in most computers)
            //			 printf("\nTime elapsed since last sensing event = %.6f ms\n", elapsed_time_ms(LastSensingTime, CurrentTime));
            if(elapsed_time_ms(LastSensingTime, CurrentTime) > max_elapsed_time_between_sensing_events_ms)
                max_elapsed_time_between_sensing_events_ms = elapsed_time_ms(LastSensingTime, CurrentTime);

            // Check if sensing event is delayed by more than 1 ms
            if(elapsed_time_ms(LastSensingTime, CurrentTime) > config.Ts*1000.0 + 1.0)
            {
                printf("\n  Warning: Real-time error -> Sensing event #%llu was delayed by %g ms", nof_EnergyValues, elapsed_time_ms(LastSensingTime, CurrentTime) - config.Ts*1000.0);
            }

            // Update time of last sensing event
            LastSensingTime = CurrentTime;

            // Debug sensing events' times
            // fprintf(stderr, "Issuing stream command...sensing event at time=%g ms.\n",elapsed_time_ms(ReceptionStartTime, CurrentTime));

            // Read data from the device
            EXECUTE_OR_GOTO(free_buffer,
                            uhd_rx_streamer_issue_stream_cmd(rx_streamer, &stream_cmd)
                            )
            current_energy = 0.0;
            num_acc_samps=0;
            while (num_acc_samps < config.N)
                {
                    // USRP receives the IQ samples and stores them in the memory buffer (of size 2*2040) as IQIQIQIQ...
                    // Each sample (I and Q) is a float value
                    size_t num_rx_samps = 0;
                    EXECUTE_OR_GOTO(close_file,
                                    uhd_rx_streamer_recv(rx_streamer, buffs_ptr, samps_per_buff, &md, 3.0, false, &num_rx_samps)
                                    )

                    uhd_rx_metadata_error_code_t error_code;
                    EXECUTE_OR_GOTO(close_file,
                                    uhd_rx_metadata_error_code(md, &error_code)
                                    )
                    if(error_code != UHD_RX_METADATA_ERROR_CODE_NONE){
                        fprintf(stderr, "Error code 0x%x was returned during streaming. Aborting.\n", *return_code);
                        goto close_file;
                    }

    //                if(num_rx_samps < samps_per_buff)
    //                {
    //                    printf("\nWarning: %zu samples lost (sensing phase)\n", samps_per_buff - num_rx_samps);
    //                }

                    // Energy calculation of the samples (energy of a single sensing event = I^2+Q^2 + I^2+Q^2 + ...)
                    for(i=0; i<config.N-num_acc_samps && i<num_rx_samps; i++)
                    {
                        current_energy += pow(((double) buff[2*i]), 2.0) + pow(((double) buff[2*i+1]), 2.0);
                    }

                    num_acc_samps += num_rx_samps;
                }

                // Save energy value and update counter (will also be used to check real-time operation)
                EnergyValues[nof_EnergyValues++] = current_energy;

//                gnuPlotter_x_vector[nof_EnergyValues-1] = nof_EnergyValues;
//                gnuPlotter_y_vector[nof_EnergyValues-1] = 0.002;


                // If in test mode, plot energy values
                if(config.TestMode)
                {
                    // Add new energy value to the vector to be plotted
                    if(nof_EnergyValues <= NOF_PLOT_POINTS)
                    {
                        gnuPlotter_x_vector[nof_EnergyValues-1] = nof_EnergyValues;
                        gnuPlotter_y_vector[nof_EnergyValues-1] = current_energy;
                    }
                    else
                    {
                        for(i=0; i<NOF_PLOT_POINTS-1; i++)
                        {
                            gnuPlotter_x_vector[i] = gnuPlotter_x_vector[i+1];
                            gnuPlotter_y_vector[i] = gnuPlotter_y_vector[i+1];
                        }
                        gnuPlotter_x_vector[NOF_PLOT_POINTS-1] = nof_EnergyValues;
                        gnuPlotter_y_vector[NOF_PLOT_POINTS-1] = current_energy;
                    }

                    // Plot latest history of energy values
                    plot_energy_history(gnuPlotter, config.DisplayEnergyStats, gnuPlotter_x_vector, gnuPlotter_y_vector, NOF_PLOT_POINTS);
                }

                // Debug
                // printf("\nEnergy = %g, Decision = %d, CurrentPeriodType = %d, Change = %d", current_energy, current_energy >= config.DecisionThreshold, CurrentPeriodType, (current_energy >= config.DecisionThreshold) != CurrentPeriodType);

                // Determine if the end of the current period has been reached
                if((current_energy >= config.DecisionThreshold) != CurrentPeriodType)
                {
                    // Update current time
                    clock_gettime(CLOCK_MONOTONIC, &CurrentTime);

                    // Calculate current period duration in seconds
                    current_period_duration = elapsed_time_ms(PeriodStartTime, CurrentTime)/1e3;

                    // Display duration of the detected period
                    if(config.ShowDetectedPeriods)
                    {
                        switch(CurrentPeriodType)
                        {
                            case -1:
                                printf("\n  Detected ???? period = %g seconds\n", current_period_duration);
                                break;
                            case 0:
                                printf("\n  Detected IDLE period = %g seconds", current_period_duration);
                                break;
                            case 1:
                                printf("\n  Detected BUSY period = %g seconds\n", current_period_duration);
                                break;
                        }
                        fflush(stdout);
                    }

                    // Add the period duration to the right vector
                    switch(CurrentPeriodType)
                    {
                        case 0:
                            detected_IdlePeriods[nof_detected_Periods] = current_period_duration;
                            if(config.DisplayEnergyStats)
                            {
                                //printf("\n this is x idle=%lld\n", nof_detected_Periods);
                                if(nof_detected_Periods == 0)//because the first detected idle period is at nof_detected_Periods=1
                                {
                                    gnuPlotter_x_vector_mean_idle[0] = 0;
                                    gnuPlotter_y_vector_mean_idle[0] = INFINITY;//=NAN;
                                    gnuPlotter_x_vector_mean_idle[nof_detected_Periods+1] = nof_detected_Periods+1;
                                    gnuPlotter_y_vector_mean_idle[nof_detected_Periods+1] = current_period_duration;
                                }
                                else
                                {
                                    if(gnuPlotter_y_vector_mean_idle[0] == INFINITY)
                                    {
                                        gnuPlotter_x_vector_mean_idle[nof_detected_Periods+1] = nof_detected_Periods+1;
                                        gnuPlotter_y_vector_mean_idle[nof_detected_Periods+1] = ((nof_detected_Periods)*gnuPlotter_y_vector_mean_idle[nof_detected_Periods]+current_period_duration)/(nof_detected_Periods+1);
                                    }
                                    else
                                    {
                                        gnuPlotter_x_vector_mean_idle[0] = 0;
                                        gnuPlotter_y_vector_mean_idle[0] = INFINITY;
                                        gnuPlotter_x_vector_mean_idle[1] = 1;
                                        gnuPlotter_y_vector_mean_idle[1] = INFINITY;
                                        gnuPlotter_x_vector_mean_idle[nof_detected_Periods+1] = nof_detected_Periods+1;
                                        gnuPlotter_y_vector_mean_idle[nof_detected_Periods+1] = current_period_duration;
                                    }
                                }
                                //plot_mean_history(gnuPlotter, gnuPlotter_x_vector_mean_idle, gnuPlotter_y_vector_mean_idle, nof_detected_Periods+1);
                            }
                            break;
                        case 1:
                            detected_BusyPeriods[nof_detected_Periods] = current_period_duration;
                            if(config.DisplayEnergyStats)
                            {
                                //printf("\n this is x busy=%lld\n", nof_detected_Periods);
                                if(nof_detected_Periods == 0)//because the first detected idle period is at nof_detected_Periods=1
                                {
                                    gnuPlotter_x_vector_mean_busy[0] = 0;
                                    gnuPlotter_y_vector_mean_busy[0] = INFINITY;
                                    gnuPlotter_x_vector_mean_busy[nof_detected_Periods+1] = nof_detected_Periods+1;
                                    gnuPlotter_y_vector_mean_busy[nof_detected_Periods+1] = current_period_duration;
                                }
                                else
                                {
                                    gnuPlotter_x_vector_mean_busy[nof_detected_Periods+1] = nof_detected_Periods+1;
                                    gnuPlotter_y_vector_mean_busy[nof_detected_Periods+1] = ((nof_detected_Periods)*gnuPlotter_y_vector_mean_busy[nof_detected_Periods]+current_period_duration)/(nof_detected_Periods+1);
                                }
                                plot_mean_history(gnuPlotter, gnuPlotter_x_vector_mean_busy, gnuPlotter_y_vector_mean_busy, nof_detected_Periods+2);

                                if(nof_detected_Periods > 0)
                                {
                                gnuPlotter_y_vector_DC[0] = INFINITY;
                                gnuPlotter_y_vector_DC[nof_detected_Periods]=gnuPlotter_y_vector_mean_busy[nof_detected_Periods+1]/(gnuPlotter_y_vector_mean_busy[nof_detected_Periods+1]+gnuPlotter_y_vector_mean_idle[nof_detected_Periods+1]);
                                plot_DC_history(gnuPlotter, gnuPlotter_x_vector_mean_busy, gnuPlotter_y_vector_DC, nof_detected_Periods+1);

                                //gnuPlotter_y_vector_distribution[(int)round(current_period_duration)]++;
                                gnuPlotter_y_vector_distribution[(int) round(current_period_duration/config.Ts)]++;
                                //printf("\nThis is yy=%d",(int) round(current_period_duration/config.Ts));
                                //printf("\nThis is y=%g\n",gnuPlotter_y_vector_distribution[(int) round(current_period_duration/config.Ts)]);
                                plot_distribution_history(gnuPlotter, gnuPlotter_y_vector_distribution, config.Ts, nof_detected_Periods);
                                }
                            }
                            // Update the counter for the next pair of idle-busy periods
                            nof_detected_Periods++;
                            break;
                    }

                    // Update the last observed period type
                    if(CurrentPeriodType == -1)
                    {
                        // This is the first period detected in this sequence
                        CurrentPeriodType = (current_energy >= config.DecisionThreshold);
                    }
                    else
                    {
                        // Switch period type
                        CurrentPeriodType = 1 - CurrentPeriodType;
                    }

                    // Start timing of new period
                    PeriodStartTime = CurrentTime;
                }

        }//end of if for Ts

        // Update current time
        clock_gettime(CLOCK_MONOTONIC, &CurrentTime);
    }//end of the main while

    printf("\n\nChannel sensing completed");

    // Check if operation was in real-time
    printf("\n\nReal-time execution check:\n");
    if(nof_EnergyValues < ((unsigned long long int) (config.ReceptionTime/config.Ts)) - 1)
    {
        printf("\n  Warning: The receiver did not operate in real-time");
        printf("\n  The duration of some estimated periods might be inaccurate");
    }
    else
    {
        printf("\n  The receiver operated in real-time");
    }
    printf("\n  Number of expected sensing events : %d", ((int) (config.ReceptionTime/config.Ts)) - 1);
    printf("\n  Number of actual   sensing events : %llu", nof_EnergyValues);
    printf("\n  Sensing events scheduled on time  : %.2f%%", 100.0 * ((double) nof_EnergyValues)/(((int) (config.ReceptionTime/config.Ts)) - 1));
    printf("\n  Selected sensing period (Ts)      : %.6f ms", config.Ts*1e3);
    printf("\n  Max time between sensing events   : %.6f ms", max_elapsed_time_between_sensing_events_ms);

    // Save the channel sensing results
    channel_sensing_results->EnergyValues = EnergyValues;
    channel_sensing_results->detected_IdlePeriods = detected_IdlePeriods;
    channel_sensing_results->detected_BusyPeriods = detected_BusyPeriods;
    channel_sensing_results->nof_EnergyValues = nof_EnergyValues;
    channel_sensing_results->nof_detected_Periods = nof_detected_Periods;

    free(gnuPlotter_x_vector_mean_idle);
    free(gnuPlotter_y_vector_mean_idle);
    free(gnuPlotter_x_vector_mean_busy);
    free(gnuPlotter_y_vector_mean_busy);
    free(gnuPlotter_y_vector_DC);
    free(gnuPlotter_y_vector_distribution);


close_file:
//    fclose(fp);

free_buffer:
    if(buff){
        free(buff);
    }
    buff = NULL;
    buffs_ptr = NULL;

free_rx_streamer:
    uhd_rx_streamer_free(&rx_streamer);

free_rx_metadata:
    uhd_rx_metadata_free(&md);

free_usrp:
    if(*return_code != EXIT_SUCCESS && usrp != NULL){
        uhd_usrp_last_error(usrp, error_string, 512);
        fprintf(stderr, "USRP reported the following error: %s\n", error_string);
    }
    uhd_usrp_free(&usrp);

free_option_strings:
    fprintf(stderr, (*return_code ? "**Failure!\n" : "\nSuccess!\n"));
    //        if(config.USRPdevice) {
    //            free(config.USRPdevice);
    //        }

}

void display_energy_statistics(struct channel_sensing_results_struct channel_sensing_results)
{
    unsigned long long int i;
    double min_energy, max_energy, average_energy;
    double cdf_x_vector[NOF_CDF_POINTS], cdf_y_vector[NOF_CDF_POINTS];

    min_energy = 1e50;
    max_energy = -1e50;
    average_energy = 0.0;

    for(i=0; i<channel_sensing_results.nof_EnergyValues; i++)
    {
        if(channel_sensing_results.EnergyValues[i] < min_energy)
            min_energy = channel_sensing_results.EnergyValues[i];

        if(channel_sensing_results.EnergyValues[i] > max_energy)
            max_energy = channel_sensing_results.EnergyValues[i];

        average_energy += channel_sensing_results.EnergyValues[i];
    }

    average_energy /= channel_sensing_results.nof_EnergyValues;

    calculate_cdf(channel_sensing_results.EnergyValues, channel_sensing_results.nof_EnergyValues, cdf_x_vector, cdf_y_vector, NOF_CDF_POINTS);

    printf("\nEnergy statistics:\n");
    printf("\n  No. energy samples     : %llu", channel_sensing_results.nof_EnergyValues);
    printf("\n  Min energy             : %.6f", min_energy);
    printf("\n  Avg energy             : %.6f", average_energy);
    printf("\n  Max energy             : %.6f", max_energy);
    printf("\n  Energy percentile 90%%  : %.6f", calculate_percentile(cdf_x_vector, cdf_y_vector, NOF_CDF_POINTS, 90.0));
    printf("\n  Energy percentile 95%%  : %.6f", calculate_percentile(cdf_x_vector, cdf_y_vector, NOF_CDF_POINTS, 95.0));
    printf("\n  Energy percentile 99%%  : %.6f", calculate_percentile(cdf_x_vector, cdf_y_vector, NOF_CDF_POINTS, 99.0));
    printf("\n");
}

void save_energy_values(struct channel_sensing_results_struct channel_sensing_results)
{
    unsigned long long int i;
    FILE *pFile;

    printf("\nSaving measured energy values to file energy_values.txt ... ");

    pFile = fopen("energy_values.txt", "w");

    if (pFile == NULL)
    {
        printf("\nError: Cannot create file to save energy values\n\n");
        exit(EXIT_FAILURE);
    }

    for(i=0; i<channel_sensing_results.nof_EnergyValues; i++)
        fprintf(pFile, "%lf\n", channel_sensing_results.EnergyValues[i]);

    fclose(pFile);

    free(channel_sensing_results.EnergyValues); // Release block of memory

    printf("Done\n");
}

void save_detected_periods(struct channel_sensing_results_struct channel_sensing_results)
{
    unsigned long long int i;
    FILE *pFile;

    printf("\nSaving detected periods to file RX_periods.txt ... ");

    pFile = fopen("RX_periods.txt", "w");

    if (pFile == NULL)
    {
        printf("\nError: Cannot create file to save detected periods\n\n");
        exit(EXIT_FAILURE);
    }

    // The first pair (i=0) will usually contain an incomplete period -> Discard
    for(i=1; i<channel_sensing_results.nof_detected_Periods; i++)
        fprintf(pFile, "%lf\t%lf\n", channel_sensing_results.detected_IdlePeriods[i], channel_sensing_results.detected_BusyPeriods[i]);

    fclose(pFile);

    free(channel_sensing_results.detected_IdlePeriods); // Release block of memory
    free(channel_sensing_results.detected_BusyPeriods); // Release block of memory

    printf("Done\n");
}

int continue_program(void)
{
    int key = '\0';
    static struct termios old_terminal_properties, new_terminal_properties;

    // Get the configuration of the current terminal
    tcgetattr(STDIN_FILENO, &old_terminal_properties);

    // Copy the configuration
    new_terminal_properties = old_terminal_properties;

    // Deactivate ICANON (waits for Enter key) and ECHO (displays the pressed key)
    new_terminal_properties.c_lflag &= ~(ICANON | ECHO);

    // Set the new terminal configuration immediately (TCSANOW)
    tcsetattr(STDIN_FILENO, TCSANOW, &new_terminal_properties);

    // Wait for user decision
    printf("\n\nDo you want to continue [y/n]? ");
    fflush(stdin);
    while(toupper(key = getchar()) != 'Y' && toupper(key) != 'N');

    // Restore the old terminal configuration immediately (TCSANOW)
    tcsetattr(STDIN_FILENO, TCSANOW, &old_terminal_properties);

    // Return user choice
    switch(toupper(key))
    {
        case 'Y':
            key = 1;
            break;
        case 'N':
            key = 0;
            break;
    }

    return key;
}

inline double elapsed_time_ms(struct timespec start_time, struct timespec end_time)
{
    return ((double) (end_time.tv_sec  - start_time.tv_sec) )*1e3 + ((double) (end_time.tv_nsec - start_time.tv_nsec))/1e6;
}

inline void plot_energy_history(FILE *gnuPlotter, int DisplayEnergyStats, double *x, double *y, int nof_points)
{
    if(DisplayEnergyStats==0)
    {
        fprintf(gnuPlotter, "set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb\"black\" behind\n");
        fprintf(gnuPlotter, "set size 1,1\n");
        fprintf(gnuPlotter, "set origin 0,0\n");
    }
    else
    {
        fprintf(gnuPlotter, "set object 1 rectangle from screen 0,0.5 to screen 0.5,1 fillcolor rgb\"black\" behind\n");
        if(y[1]==0)
        fprintf(gnuPlotter, "set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb\"black\" behind\n");

        fprintf(gnuPlotter, "set size 0.5,0.5\n");
        fprintf(gnuPlotter, "set origin 0,0.5\n");
        //fprintf(gnuPlotter, "clear\n");
    }

    int i;
    static double ymin = 1e50, ymax = -1e50;

    for(i=0; i<nof_points; i++)
    {
        if(y[i] > ymax)
        {
            ymax = y[i];
            ymax *= 1.10;
        }
    }

    // This adjustment provides better appearance for the graph of energy values
    ymin = -ymax*0.10;

    // Set the abscissa range
    fprintf(gnuPlotter, "set xrange [%g:%g]\n", x[0], x[nof_points-1]);

    // Set the ordinate range
    fprintf(gnuPlotter, "set yrange [%g:%g]\n", 0.0, ymax);

    // Sets the label for the X axis
    fprintf(gnuPlotter, "set xlabel \"Sensing event number\" font \"Sans,12\" textcolor rgb \"cyan\"\n");
    // Sets the label for the Y axis
    fprintf(gnuPlotter, "set ylabel \"Energy value\" font \"Sans,12\" textcolor rgb \"cyan\"\n");
    // Sets the title for the plot
    fprintf(gnuPlotter, "set title \"History of sensed energy values\" font \"Sans,14\" textcolor rgb \"web-blue\"\n");

    // Send to gnuplot the values to be plotted (one by one)
    fprintf(gnuPlotter, "plot '-'\n");
    for (i = 0; i < nof_points; i++)
    {
        fprintf(gnuPlotter, "%g %g\n", x[i], y[i]);
    }
    fprintf(gnuPlotter, "e\n");

    // Flush commands
    fflush(gnuPlotter);
}

inline void plot_mean_history(FILE *gnuPlotter_mean, double *x, double *y, int nof_points)
{
    fprintf(gnuPlotter_mean, "set object 1 rectangle from screen 0.5,0.5 to screen 1,1 fillcolor rgb\"black\" behind\n");
    fprintf(gnuPlotter_mean, "set size 0.5,0.5\n");
    fprintf(gnuPlotter_mean, "set origin 0.5,0.5\n");
    int i;
    static double ymax = 0;

    for(i=1; i<nof_points; i++)
    {
        if(y[i] > ymax)
        {
            ymax = y[i];
            ymax *= 1.10;
        }
    }

    // Set the abscissa range
    fprintf(gnuPlotter_mean, "set xrange [%d:%d]\n", 0, nof_points);

    // Set the ordinate range
    fprintf(gnuPlotter_mean, "set yrange [%d:%g]\n", 0, ymax);

    fprintf(gnuPlotter_mean, "set xlabel \"Number of detected busy periods\" font \"Sans,12\" textcolor rgb \"cyan\"\n");
    // Sets the label for the Y axis
    fprintf(gnuPlotter_mean, "set ylabel \"Mean of the busy periods\" font \"Sans,12\" textcolor rgb \"cyan\"\n");
    // Sets the title for the plot
    fprintf(gnuPlotter_mean, "set title \"Mean of the busy periods (real-time)\" font \"Sans,14\" textcolor rgb \"web-blue\"\n");


    // Send to gnuplot the values to be plotted (one by one)
    fprintf(gnuPlotter_mean, "plot '-'\n");
    for (i = 1; i < nof_points; i++)
    {
        fprintf(gnuPlotter_mean, "%g %g\n", x[i], y[i]);
    }
    fprintf(gnuPlotter_mean, "e\n");

    // Flush commands
    fflush(gnuPlotter_mean);
}

inline void plot_DC_history(FILE *gnuPlotter_DC, double *x, double *y, int nof_points)
{
    fprintf(gnuPlotter_DC, "set object 1 rectangle from screen 0,0 to screen 0.5,0.5 fillcolor rgb\"black\" behind\n");
    fprintf(gnuPlotter_DC, "set size 0.5,0.5\n");
    fprintf(gnuPlotter_DC, "set origin 0,0\n");

    // Set the abscissa range
    fprintf(gnuPlotter_DC, "set xrange [%d:%d]\n", 0, nof_points);

    // Set the ordinate range
    fprintf(gnuPlotter_DC, "set yrange [%d:%d]\n", 0, 1);

    fprintf(gnuPlotter_DC, "set xlabel \"Number of periods\" font \"Sans,12\" textcolor rgb \"cyan\"\n");
    // Sets the label for the Y axis
    fprintf(gnuPlotter_DC, "set ylabel \"Duty cycle of the periods\" font \"Sans,12\" textcolor rgb \"cyan\"\n");
    // Sets the title for the plot
    fprintf(gnuPlotter_DC, "set title \"Duty Cycle of the periods (real-time)\" font \"Sans,14\" textcolor rgb \"web-blue\"\n");


    // Send to gnuplot the values to be plotted (one by one)
    fprintf(gnuPlotter_DC, "plot '-'\n");
    for (int i = 0; i < nof_points; i++)
    {
        fprintf(gnuPlotter_DC, "%g %g\n", x[i], y[i]);
    }
    fprintf(gnuPlotter_DC, "e\n");

    // Flush commands
    fflush(gnuPlotter_DC);
}

inline void plot_distribution_history(FILE *gnuPlotter_distribution, double *y, double ts, int nof_points)
{
    fprintf(gnuPlotter_distribution, "set object 1 rectangle from screen 0.5,0 to screen 1,0.5 fillcolor rgb\"black\" behind\n");
    fprintf(gnuPlotter_distribution, "set size 0.5,0.5\n");
    fprintf(gnuPlotter_distribution, "set origin 0.5,0\n");
//    int xmin=0, xmax=1;
//    for (int i = 0; i < 1000; i++)
//    {
//        if(y[i]==0 && xmin==i)
//            xmin++;
//         if(y[i]>0)
//            xmax=i;
//    }

    int cdf=0,xmax=0,i=0;
    cdf+=y[i];
    while(cdf<nof_points)
    {
        //printf("%f\n",y[i]);
        xmax=i;
        i++;
        cdf+=y[i];
    }
    xmax=i;

    // Set the abscissa range
    //fprintf(gnuPlotter_distribution, "set xrange [%d:%d]\n", xmin-5, xmax+5);
    fprintf(gnuPlotter_distribution, "set xrange [%d:%f]\n", 0, (xmax+5)*ts);

    // Set the ordinate range
    fprintf(gnuPlotter_distribution, "set yrange [%d:%d]\n", 0, 1);

    fprintf(gnuPlotter_distribution, "set xlabel \"Period duration (seconds)\" font \"Sans,12\" textcolor rgb \"cyan\"\n");
    // Sets the label for the Y axis
    fprintf(gnuPlotter_distribution, "set ylabel \"PDF of the busy periods\" font \"Sans,12\" textcolor rgb \"cyan\"\n");
    // Sets the title for the plot
    fprintf(gnuPlotter_distribution, "set title \"Distribution of the busy periods (real-time)\" font \"Sans,14\" textcolor rgb \"web-blue\"\n");

    fprintf(gnuPlotter_distribution, "set boxwidth %f\n", ts);
    fprintf(gnuPlotter_distribution, "set style fill solid border -1\n");

    // Send to gnuplot the values to be plotted (one by one)
    fprintf(gnuPlotter_distribution, "plot '-' with boxes\n");
    for (int i = 0; i < xmax+5; i++)
    {
        fprintf(gnuPlotter_distribution, "%g %g\n", i*ts, y[i]/nof_points);
    }
    fprintf(gnuPlotter_distribution, "e\n");

    // Flush commands
    fflush(gnuPlotter_distribution);
}

void calculate_cdf(double *vector, unsigned long long int vector_size, double *cdf_x_vector, double *cdf_y_vector, int cdf_size)
{
    unsigned long long int i, counter;
    int j;
    double x_min, x_max;

    // Find the minimum and maximum values of the input vector
    x_min = 1e50;
    x_max = -1e50;
    for(i=0; i<vector_size; i++)
    {
        if(vector[i] < x_min)
            x_min = vector[i];
        if(vector[i] > x_max)
            x_max = vector[i];
    }

    // Calculate CDF
    for(j=0; j<cdf_size; j++)
    {
        // Create vector of values for the abscissa
        cdf_x_vector[j] = x_min + j*(x_max - x_min)/(cdf_size-1);

        // Count the frequency of each energy value (value for the ordinate)
        counter = 0;
        for(i=0; i<vector_size; i++)
        {
            if(vector[i] <= cdf_x_vector[j])
                counter++;
        }
        cdf_y_vector[j] = ((double) counter) / ((double) vector_size);
    }
}

double calculate_percentile(double *cdf_x_vector, double *cdf_y_vector, int cdf_size, double percentile)
{
    int i, min_index = 0;
    double abs_difference[cdf_size], min_abs_difference = 1e3;

    for(i=0; i<cdf_size; i++)
    {
        abs_difference[i] = fabs(cdf_y_vector[i] - percentile/100.0);
        if(abs_difference[i] < min_abs_difference)
        {
            min_abs_difference = abs_difference[i];
            min_index = i;
        }
    }

    return cdf_x_vector[min_index];
}


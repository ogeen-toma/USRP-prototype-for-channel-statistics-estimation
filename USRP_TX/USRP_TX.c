
/***********************************************************************
 *
 * USRP-Based Prototype for Real-Time Estimation of Channel Activity Statistics 
 *				in Smart Spectrum Sharing
 *
 ***********************************************************************
 *
 *                 		TRANSMITTER PROGRAM
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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <termios.h>
#include <time.h>
#include <unistd.h>
#include <ctype.h>

#include <uhd.h>
#include <signal.h>

#define EXECUTE_OR_GOTO(label, ...) \
if(__VA_ARGS__){ \
*return_code = EXIT_FAILURE; \
goto label; \
}

// Current version
#define VERSION 1.00

// List of distributions (see Table I of IEEE TVT 2013 paper)
#define TEST			0	// Used for testing (1 sec ON + 1 sec OFF)
#define EXPONENTIAL		1
#define GEN_EXPONENTIAL	2
#define PARETO			3
#define GEN_PARETO		4
#define LOG_NORMAL		5
#define GAMMA			6
#define WEIBULL			7

// Maximum number of characters per line in the configuration file
#define MAX_LINE_LENGTH_IN_CONFIG_FILE	100

// Structure with all configuration parameters
struct config_struct
{
    int Load_or_Generate;
    int SavePeriods;
    unsigned long long int nofPeriods;
    int Distribution_T0;
    double mu_0;
    double lambda_0;
    double alpha_0;
    int Distribution_T1;
    double mu_1;
    double lambda_1;
    double alpha_1;
    int ShowTransmittedPeriods;
    char USRPdevice[20];
    double CenterFrequency;
    double SampleRate;
    double TunerGain;
};

// Function prototypes
void display_header(void);
void set_default_config_params(struct config_struct *config);
void load_config_file(int argc, char *argv[], struct config_struct *config);
void check_config_params(struct config_struct config);
void show_config_params(struct config_struct config);
void allocate_memory_for_periods(double **idle_periods, double **busy_periods, unsigned long long int nofPeriods);
void free_memory_for_periods(double *idle_periods, double *busy_periods);
//void reset_transmitters(void);
void reset_vector(double *vector, unsigned long long int vector_size, double value);
void load_periods(double *idle_periods, double *busy_periods, unsigned long long int nofPeriods);
void generate_periods(double *idle_periods, double *busy_periods, struct config_struct config);
void save_periods(double *idle_periods, double *busy_periods, unsigned long long int nofPeriods);
void show_statistics(const double *idle_periods, const double *busy_periods, unsigned long long int nofPeriods);
void transmit_sequence(const double *idle_periods, const double *busy_periods, unsigned long long int nofPeriods, int ShowTransmittedPeriods, char USRPdevice[20], double CenterFrequency, double SampleRate, double TunerGain, int *return_code);
double elapsed_time_ms(struct timespec start_time, struct timespec end_time);
int continue_program(void);
// double erf_inv_maclaurin(double x);
double erf_inv(double y);
double lower_gammainc(double a, double x);
double lower_gammainc_inv(double a, double y);

bool stop_signal_called = false;
void sigint_handler(int code){
    (void)code;
    stop_signal_called = true;
}


int main(int argc, char *argv[])
{
    // Period durations in SECONDS
    double *idle_periods, *busy_periods;

    // Configuration structure
    struct config_struct config;

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


    // Allocate memory for the idle and busy periods
    allocate_memory_for_periods(&idle_periods, &busy_periods, config.nofPeriods);


    // Reset period durations
    reset_vector(idle_periods, config.nofPeriods, 0.0);
    reset_vector(busy_periods, config.nofPeriods, 0.0);

    // Load or generate periods
    switch(config.Load_or_Generate)
    {
        case 0:		// Load periods from "TX_periods.txt" file
            load_periods(idle_periods, busy_periods, config.nofPeriods);
            break;

        case 1:		// Generate random periods
            generate_periods(idle_periods, busy_periods, config);
            if(config.SavePeriods)
                save_periods(idle_periods, busy_periods, config.nofPeriods);
            break;

        default:
            printf("\nError: Wrong Load_or_Generate option\n\n");
            exit(EXIT_FAILURE);
    }

    /*
     // Used for debugging
     unsigned long long int i;
     for(i=0; i<config.nofPeriods; i++)
     printf("Idle[%llu/%llu] = %lf\tBusy[%llu/%llu] = %lf\n", i+1, config.nofPeriods, idle_periods[i], i+1, config.nofPeriods, busy_periods[i]);
     */

    // Show statistics of the sequence to be transmitted
    show_statistics(idle_periods, busy_periods, config.nofPeriods);

    // Final confirmation before starting experiment
    printf("\nReady to start transmission with USRP %s...", config.USRPdevice);

    int return_code = EXIT_SUCCESS;
    //if(continue_program())
      if(1)
    {
        // Transmit sequence
        transmit_sequence(idle_periods, busy_periods, config.nofPeriods, config.ShowTransmittedPeriods, config.USRPdevice, config.CenterFrequency, config.SampleRate, config.TunerGain,&return_code);

        // Free memory for the idle and busy periods
        free_memory_for_periods(idle_periods, busy_periods);

        if (return_code == EXIT_FAILURE){
            printf("\nProgram execution failed.\n\n");
            return 0;
        }

        printf("\nProgram execution finished successfully\n\n");
    }
    else
    {
        // Free memory for the idle and busy periods
        free_memory_for_periods(idle_periods, busy_periods);

        printf("\nProgram execution finished by the user\n\n");
    }

    return 0;
}

void display_header(void)
{
    printf("\n########################################################");
    printf("\n#                                                      #");
    printf("\n#                USRP_TX (B200mini) v%.2f              #", VERSION);
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

    config->Load_or_Generate = 1;			// Generate random periods
    config->SavePeriods = 0;				// Do not save the generated random periods
    config->nofPeriods = 10000;				// Number of periods to be generated
    config->Distribution_T0 = GEN_PARETO;	// Distribution for idle periods
    config->mu_0 = 3.5310;					// Location parameter for idle periods
    config->lambda_0 = 2.6272;				// Scale parameter for idle periods
    config->alpha_0 = 0.2119;				// Shape parameter for idle periods
    config->Distribution_T1 = GEN_PARETO;	// Distribution for busy periods
    config->mu_1 = 3.5470;					// Location parameter for busy periods
    config->lambda_1 = 10.7968;				// Scale parameter for busy periods
    config->alpha_1 = 0.1929;				// Shape parameter for busy periods
    config->ShowTransmittedPeriods = 1;		// Show the transmitted periods
    //    config->USRPdevice[10] = "b200";            // Specify the USRP device type
    strncpy(config->USRPdevice,"serial=319473B",16);  // Specify the USRP device serial number
    config->CenterFrequency = 433.92e6;     // Select the centre frequency in Hz
    config->SampleRate = 1e6;               // Select the sample rate in Hz
    config->TunerGain = 15;                 // Select the gain
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

            if(strcasecmp(parameter_name, "Load_or_Generate") == 0)
            {
                config->Load_or_Generate = atoi(parameter_value);
            }

            if(strcasecmp(parameter_name, "SavePeriods") == 0)
            {
                config->SavePeriods = atoi(parameter_value);
            }

            if(strcasecmp(parameter_name, "nofPeriods") == 0)
            {
                config->nofPeriods = strtoull(parameter_value, NULL, 10);
            }

            if(strcasecmp(parameter_name, "Distribution_T0") == 0)
            {
                config->Distribution_T0 = atoi(parameter_value);
            }

            if(strcasecmp(parameter_name, "mu_0") == 0)
            {
                config->mu_0 = atof(parameter_value);
            }

            if(strcasecmp(parameter_name, "lambda_0") == 0)
            {
                config->lambda_0 = atof(parameter_value);
            }

            if(strcasecmp(parameter_name, "alpha_0") == 0)
            {
                config->alpha_0 = atof(parameter_value);
            }

            if(strcasecmp(parameter_name, "Distribution_T1") == 0)
            {
                config->Distribution_T1 = atoi(parameter_value);
            }

            if(strcasecmp(parameter_name, "mu_1") == 0)
            {
                config->mu_1 = atof(parameter_value);
            }

            if(strcasecmp(parameter_name, "lambda_1") == 0)
            {
                config->lambda_1 = atof(parameter_value);
            }

            if(strcasecmp(parameter_name, "alpha_1") == 0)
            {
                config->alpha_1 = atof(parameter_value);
            }

            if(strcasecmp(parameter_name, "ShowTransmittedPeriods") == 0)
            {
                config->ShowTransmittedPeriods = atoi(parameter_value);
            }

            if(strcasecmp(parameter_name, "USRPdevice") == 0)
            {
                strncpy(config->USRPdevice,"serial=",16);
                strcat(config->USRPdevice,parameter_value);
            }

            if(strcasecmp(parameter_name, "CenterFrequency") == 0)
            {
                config->CenterFrequency = atof(parameter_value);
            }

            if(strcasecmp(parameter_name, "SampleRate") == 0)
            {
                config->SampleRate = atof(parameter_value);
            }

            if(strcasecmp(parameter_name, "TunerGain") == 0)
            {
                config->TunerGain = atof(parameter_value);
            }
        }
    }

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

    if(config.Load_or_Generate != 0 && config.Load_or_Generate != 1)
    {
        printf("\nError: Load_or_Generate must be either 0 or 1\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.SavePeriods != 0 && config.SavePeriods != 1)
    {
        printf("\nError: SavePeriods must be either 0 or 1\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.nofPeriods <= 0)
    {
        printf("\nError: The number of periods (nofPeriods) must be greater than zero\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.Distribution_T0 < 0 || config.Distribution_T0 > 7)
    {
        printf("\nError: The distribution for idle periods must be an integer number from 0 to 7\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.mu_0 <= 0.0)
    {
        printf("\nError: mu_0 must be greater than zero\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.lambda_0 <= 0.0)
    {
        printf("\nError: lambda_0 must be greater than zero\n\n");
        exit(EXIT_FAILURE);
    }

    if((config.Distribution_T0 == GEN_EXPONENTIAL && config.alpha_0 <= 0.0) ||
       (config.Distribution_T0 == PARETO          && config.alpha_0 <= 2.0) ||
       (config.Distribution_T0 == GEN_PARETO      && config.alpha_0 >= 0.5) ||
       (config.Distribution_T0 == GAMMA           && config.alpha_0 <= 0.0) ||
       (config.Distribution_T0 == WEIBULL         && config.alpha_0 <= 0.0))
    {
        printf("\nError: alpha_0 is incorrect for the selected distribution\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.Distribution_T1 < 0 || config.Distribution_T1 > 7)
    {
        printf("\nError: The distribution for busy periods must be an integer number from 0 to 7\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.mu_1 <= 0.0)
    {
        printf("\nError: mu_1 must be greater than zero\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.lambda_1 <= 0.0)
    {
        printf("\nError: lambda_1 must be greater than zero\n\n");
        exit(EXIT_FAILURE);
    }

    if((config.Distribution_T1 == GEN_EXPONENTIAL && config.alpha_1 <= 0.0) ||
       (config.Distribution_T1 == PARETO          && config.alpha_1 <= 2.0) ||
       (config.Distribution_T1 == GEN_PARETO      && config.alpha_1 >= 0.5) ||
       (config.Distribution_T1 == GAMMA           && config.alpha_1 <= 0.0) ||
       (config.Distribution_T1 == WEIBULL         && config.alpha_1 <= 0.0))
    {
        printf("\nError: alpha_1 is incorrect for the selected distribution\n\n");
        exit(EXIT_FAILURE);
    }

    if((config.Distribution_T0 == TEST && config.Distribution_T1 != TEST) ||
       (config.Distribution_T1 == TEST && config.Distribution_T0 != TEST))
    {
        printf("\nError: Only one distribution is in TEST mode. For TEST mode, both distributions should be TEST.\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.ShowTransmittedPeriods != 0 && config.ShowTransmittedPeriods != 1)
    {
        printf("\nError: ShowTransmittedPeriods must be either 0 or 1\n\n");
        exit(EXIT_FAILURE);
    }

//    if(strcasecmp(config.USRPdevice, "serial=319473B") != 0)
//    {
//        printf("\nError: Enter the correct serial number of the USRP\n\n");
//        exit(EXIT_FAILURE);
//    }

    if(config.CenterFrequency < 70e6 || config.CenterFrequency > 6e9)
    {
        printf("\nError: USRP B200mini center frequency must be between [70 MHz – 6 GHz]\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.SampleRate <= 200e3 || config.SampleRate > 56e6)
    {
        printf("\nError: USRP B200mini sample rate must be between [200 kHz – 56 MHz]\n\n");
        exit(EXIT_FAILURE);
    }

    if(config.TunerGain < 0 || config.TunerGain > 89.8)
    {
        printf("\nError: USRP B200mini transmit frontend gain must be between [0 dB – 89.8 dB]\n\n");
        exit(EXIT_FAILURE);
    }
    printf("Done\n");
}

void show_config_params(struct config_struct config)
{
    printf("\nProgram operation parameters:\n");
    printf("\n  Load or generate periods    : %s", config.Load_or_Generate ? "Generate" : "Load");
    if(config.Load_or_Generate)
    {
        printf("\n  Save periods to file        : %s", config.SavePeriods ? "YES" : "NO");
        printf("\n  No. periods to be generated : %llu", config.nofPeriods);
        switch(config.Distribution_T0)
        {
            case TEST:
                printf("\n  Distribution IDLE periods   : Test mode");
                break;
            case EXPONENTIAL:
                printf("\n  Distribution IDLE periods   : Exponential");
                break;
            case GEN_EXPONENTIAL:
                printf("\n  Distribution IDLE periods   : Generalised exponential");
                break;
            case PARETO:
                printf("\n  Distribution IDLE periods   : Pareto");
                break;
            case GEN_PARETO:
                printf("\n  Distribution IDLE periods   : Generalised Pareto");
                break;
            case LOG_NORMAL:
                printf("\n  Distribution IDLE periods   : Log-normal");
                break;
            case GAMMA:
                printf("\n  Distribution IDLE periods   : Gamma");
                break;
            case WEIBULL:
                printf("\n  Distribution IDLE periods   : Weibull");
                break;
        }
        if(config.Distribution_T0 != TEST)
        {
            printf("\n  ├──────────────── mu_0      : %g", config.mu_0);
            printf("\n  ├──────────────── lambda_0  : %g", config.lambda_0);
            printf("\n  └──────────────── alpha_0   : %g", config.alpha_0);
        }
        switch(config.Distribution_T1)
        {
            case TEST:
                printf("\n  Distribution BUSY periods   : Test mode");
                break;
            case EXPONENTIAL:
                printf("\n  Distribution BUSY periods   : Exponential");
                break;
            case GEN_EXPONENTIAL:
                printf("\n  Distribution BUSY periods   : Generalised exponential");
                break;
            case PARETO:
                printf("\n  Distribution BUSY periods   : Pareto");
                break;
            case GEN_PARETO:
                printf("\n  Distribution BUSY periods   : Generalised Pareto");
                break;
            case LOG_NORMAL:
                printf("\n  Distribution BUSY periods   : Log-normal");
                break;
            case GAMMA:
                printf("\n  Distribution BUSY periods   : Gamma");
                break;
            case WEIBULL:
                printf("\n  Distribution BUSY periods   : Weibull");
                break;
        }
        if(config.Distribution_T1 != TEST)
        {
            printf("\n  ├──────────────── mu_1      : %g", config.mu_1);
            printf("\n  ├──────────────── lambda_1  : %g", config.lambda_1);
            printf("\n  └──────────────── alpha_1   : %g", config.alpha_1);
        }
    }
    else
    {
        printf("\n  No. periods to be loaded    : %llu", config.nofPeriods);
    }
    printf("\n  Show transmitted periods    : %s", config.ShowTransmittedPeriods ? "YES" : "NO");
    printf("\n  USRP serial number          : %s", config.USRPdevice);
    printf("\n  Selected centre frequency   : %g", config.CenterFrequency);
    printf("\n  Selected sample rate        : %g", config.SampleRate);
    printf("\n  Selected gain               : %g", config.TunerGain);
}

void allocate_memory_for_periods(double **idle_periods, double **busy_periods, unsigned long long int nofPeriods)
{
    printf("\nAllocating memory for the periods ... ");

    *idle_periods = (double *) malloc(nofPeriods * sizeof(double));
    if(*idle_periods == NULL)
    {
        printf("\nError: Failed to allocate memory for idle periods\n\n");
        exit(EXIT_FAILURE);
    }

    *busy_periods = (double *) malloc(nofPeriods * sizeof(double));
    if(*busy_periods == NULL)
    {
        printf("\nError: Failed to allocate memory for busy periods\n\n");
        exit(EXIT_FAILURE);
    }

    printf("Done\n");
}

void free_memory_for_periods(double *idle_periods, double *busy_periods)
{
    printf("\n\nReleasing allocated memory ... ");

    free(idle_periods);
    free(busy_periods);

    printf("Done\n");
}

void reset_vector(double *vector, unsigned long long int vector_size, double value)
{
    unsigned long long int i;

    for(i=0; i<vector_size; i++)
        vector[i] = value;
}

void load_periods(double *idle_periods, double *busy_periods, unsigned long long int nofPeriods)
{
    unsigned long long int i;
    FILE *pFile;

    printf("\nLoading periods from file TX_periods.txt ... ");

    // Attempt to open the file
    pFile = fopen("TX_periods.txt", "r");

    if(pFile == NULL)
    {
        printf("\nError: Failed to open TX_periods.txt\n\n");
        exit(EXIT_FAILURE);
    }

    for(i=0; i<nofPeriods; i++)
    {
        if(fscanf(pFile, "%lf\t%lf\n", &idle_periods[i], &busy_periods[i]) == EOF)
        {
            printf("\nWarning: End of file TX_periods.txt was reached before loading the expected number of periods.");
            printf("\nLoaded %llu out of %llu periods (remaining periods will be set to zero)\n", i, nofPeriods);
        }
    }

    // Close the file
    fclose(pFile);

    printf("Done\n");
}

void generate_periods(double *idle_periods, double *busy_periods, struct config_struct config)
{
    // This function generates random periods according to the
    // desired distribution based on the inversion method

    unsigned long long int i;
    double unif_random_number;
    time_t time_info;

    // Generate periods
    printf("\nGenerating periods ... ");
    fflush(stdout);

    // Use a different seed every time a new sequence is generated
    srand((unsigned int) time(&time_info));

    for(i=0; i<config.nofPeriods; i++)
    {
        // Generate IDLE period
        unif_random_number = ((double) rand()) / ((double) RAND_MAX);

        switch(config.Distribution_T0)
        {
            case TEST:
                idle_periods[i] = 1.0;
                break;

            case EXPONENTIAL:
                idle_periods[i] = config.mu_0 - (1.0/config.lambda_0)*log(1.0 - unif_random_number);
                break;

            case GEN_EXPONENTIAL:
                idle_periods[i] = config.mu_0 - (1.0/config.lambda_0)*log(1.0 - pow(unif_random_number, 1.0/config.alpha_0));
                break;

            case PARETO:
                idle_periods[i] = config.lambda_0 / pow(1.0 - unif_random_number, 1.0/config.alpha_0);
                break;

            case GEN_PARETO:
                idle_periods[i] = config.mu_0 + (config.lambda_0/config.alpha_0)*(pow(1.0 - unif_random_number,-config.alpha_0) - 1.0);
                break;

            case LOG_NORMAL:
                idle_periods[i] = exp(config.mu_0 + sqrt(2.0)*config.lambda_0*erf_inv(2.0*unif_random_number - 1.0));
                break;

            case GAMMA:
                idle_periods[i] = config.mu_0 + config.lambda_0*lower_gammainc_inv(config.alpha_0, unif_random_number*tgamma(config.alpha_0));
                break;

            case WEIBULL:
                idle_periods[i] = config.mu_0 + config.lambda_0*pow(-log(1.0 - unif_random_number), 1.0/config.alpha_0);
                break;

            default:
                printf("\nError: Wrong distribution for IDLE periods\n\n");
                exit(EXIT_FAILURE);
        }

        // Generate BUSY period
        unif_random_number = ((double) rand()) / ((double) RAND_MAX);

        switch(config.Distribution_T1)
        {
            case TEST:
                busy_periods[i] = 1.0;
                break;

            case EXPONENTIAL:
                busy_periods[i] = config.mu_1 - (1.0/config.lambda_1)*log(1.0 - unif_random_number);
                break;

            case GEN_EXPONENTIAL:
                busy_periods[i] = config.mu_1 - (1.0/config.lambda_1)*log(1.0 - pow(unif_random_number, 1.0/config.alpha_1));
                break;

            case PARETO:
                busy_periods[i] = config.lambda_1 / pow(1.0 - unif_random_number, 1.0/config.alpha_1);
                break;

            case GEN_PARETO:
                busy_periods[i] = config.mu_1 + (config.lambda_1/config.alpha_1)*(pow(1.0 - unif_random_number,-config.alpha_1) - 1.0);
                break;

            case LOG_NORMAL:
                busy_periods[i] = exp(config.mu_1 + sqrt(2.0)*config.lambda_1*erf_inv(2.0*unif_random_number - 1.0));
                break;

            case GAMMA:
                busy_periods[i] = config.mu_1 + config.lambda_1*lower_gammainc_inv(config.alpha_1, unif_random_number*tgamma(config.alpha_1));
                break;

            case WEIBULL:
                busy_periods[i] = config.mu_1 + config.lambda_1*pow(-log(1.0 - unif_random_number), 1.0/config.alpha_1);
                break;

            default:
                printf("ERROR: Wrong distribution for IDLE periods");
                exit(1);
        }
    }

    printf("Done\n");
}

void save_periods(double *idle_periods, double *busy_periods, unsigned long long int nofPeriods)
{
    unsigned long long int i;
    FILE *pFile;

    printf("\nSaving generated periods to file TX_periods.txt ... ");

    // Attempt to open the file
    pFile = fopen("TX_periods.txt", "w");

    if(pFile == NULL)
    {
        printf("\nError: Failed to open TX_periods.txt\n\n");
        exit(EXIT_FAILURE);
    }

    for(i=0; i<nofPeriods; i++)
    {
        fprintf(pFile, "%lf\t%lf\n", idle_periods[i], busy_periods[i]);
    }

    // Close the file
    fclose(pFile);

    printf("Done\n");
}

void show_statistics(const double *idle_periods, const double *busy_periods, unsigned long long int nofPeriods)
{
    // The selected values for the distributions' parameters lead to
    // certain theoretical values for the minimum/maximum period
    // durations as well as their average and standard deviation, etc.
    // However, the sequence has a limited number of periods and may
    // have slightly different statistical properties (the higher the
    // number of periods, the closer to the theoretical values will be).
    // This function shows statistics for the real transmitted sequence.
    // Notice that the best estimation the receiver can perform is these
    // values, not the ones selected in the configuration parameters.
    // The estimation of the distribution parameters is more complex and
    // can be done following different methods, which is not done here.

    unsigned long long int i;
    double min_T0 = 1.0e100, max_T0 = 0.0, avg_T0 = 0.0, std_T0 = 0.0;
    double min_T1 = 1.0e100, max_T1 = 0.0, avg_T1 = 0.0, std_T1 = 0.0;
    double hours, minutes, seconds, total_seconds;

    // Minimum, maximum and average
    for(i=0; i<nofPeriods; i++)
    {
        // Minimum
        if(idle_periods[i] < min_T0)
            min_T0 = idle_periods[i];
        if(busy_periods[i] < min_T1)
            min_T1 = busy_periods[i];

        // Maximum
        if(idle_periods[i] > max_T0)
            max_T0 = idle_periods[i];
        if(busy_periods[i] > max_T1)
            max_T1 = busy_periods[i];

        // Average
        avg_T0 = avg_T0 + idle_periods[i];
        avg_T1 = avg_T1 + busy_periods[i];
    }

    // Average
    avg_T0 = avg_T0 / ((double) nofPeriods);
    avg_T1 = avg_T1 / ((double) nofPeriods);

    // Standard deviation
    for(i=0; i<nofPeriods; i++)
    {
        std_T0 = std_T0 + pow(idle_periods[i] - avg_T0, 2.0);
        std_T1 = std_T1 + pow(busy_periods[i] - avg_T1, 2.0);
    }
    std_T0 = sqrt(std_T0 / (((double) nofPeriods) - 1.0));
    std_T1 = sqrt(std_T1 / (((double) nofPeriods) - 1.0));

    // Calculate the total duration of the transmission
    total_seconds = (avg_T0 + avg_T1) * nofPeriods;
    hours = floor(total_seconds/3600.0);
    minutes = floor(((unsigned long int) total_seconds % 3600)/60.0);
    seconds = ((unsigned long int) total_seconds % 3600) % 60;

    // Show statistics
    printf("\nStatistics of the sequence:\n\n");
    printf("  Total number of periods = %llu\n", nofPeriods);
    printf("  IDLE periods:\n");
    printf("       Minimum T0 = %g seconds\n", min_T0);
    printf("       Maximum T0 = %g seconds\n", max_T0);
    printf("       Average T0 = %g seconds\n", avg_T0);
    printf("       Std dev T0 = %g seconds\n", std_T0);
    printf("       Var     T0 = %g seconds^2\n", std_T0*std_T0);
    printf("  BUSY periods:\n");
    printf("       Minimum T1 = %g seconds\n", min_T1);
    printf("       Maximum T1 = %g seconds\n", max_T1);
    printf("       Average T1 = %g seconds\n", avg_T1);
    printf("       Std dev T1 = %g seconds\n", std_T1);
    printf("       Var     T1 = %g seconds^2\n", std_T1*std_T1);
    printf("  Average duty cycle = %.2f\n", avg_T1 / (avg_T0 + avg_T1));
    printf("  Total transmission duration [HH:MM:SS] = %.2d:%.2d:%.2d (%.0f seconds)\n", (int) hours, (int) minutes, (int) seconds, total_seconds);
}

void transmit_sequence(const double *idle_periods, const double *busy_periods, unsigned long long int nofPeriods, int ShowTransmittedPeriods, char USRPdevice[20], double CenterFrequency, double SampleRate, double TunerGain, int *return_code)
{
    unsigned long long int i;
    struct timespec TransmissionStartTime, TransmissionEndTime;
    double hours, minutes, seconds, expected_total_seconds, actual_total_seconds;
    struct timespec *tm_idle_periods;

    /////////////////////////////////////////////////// USRP configuration ////////////////////////////////////////////////
    size_t channel = 0;
    char error_string[512];

    if(uhd_set_thread_priority(uhd_default_thread_priority, true)){
        fprintf(stderr, "Unable to set thread priority. Continuing anyway.\n");
    }

    //    if (!USRPdevice)
    //        USRPdevice = strdup("");

    // Create USRP
    uhd_usrp_handle usrp;
    fprintf(stderr, "\n\nConnecting to USRP \"%s\"...\n", USRPdevice);
    EXECUTE_OR_GOTO(free_option_strings,
                    uhd_usrp_make(&usrp, USRPdevice)
                    )

    // Create TX streamer
    uhd_tx_streamer_handle tx_streamer;
    EXECUTE_OR_GOTO(free_usrp,
                    uhd_tx_streamer_make(&tx_streamer)
                    )

    // Create TX metadata
    uhd_tx_metadata_handle md;
    EXECUTE_OR_GOTO(free_tx_streamer,
                    uhd_tx_metadata_make(&md, false, 0, 0.1, true, false)
                    )

    // Create other necessary structs
    uhd_tune_request_t tune_request = {
        .target_freq = CenterFrequency,
        .rf_freq_policy = UHD_TUNE_REQUEST_POLICY_AUTO,
        .dsp_freq_policy = UHD_TUNE_REQUEST_POLICY_AUTO
    };
    uhd_tune_result_t tune_result;

    uhd_stream_args_t stream_args = {
        .cpu_format = "fc32",
        .otw_format = "sc16",
        .args = "",
        .channel_list = &channel,
        .n_channels = 1
    };

    size_t samps_per_buff;
    float* buff_busy = NULL;
    const void** buffs_ptr_busy = NULL;

    // Set rate
    fprintf(stderr, "Setting TX Rate: %f...\n", SampleRate);
    EXECUTE_OR_GOTO(free_tx_metadata,
                    uhd_usrp_set_tx_rate(usrp, SampleRate, channel)
                    )

    // See what rate actually is
    EXECUTE_OR_GOTO(free_tx_metadata,
                    uhd_usrp_get_tx_rate(usrp, channel, &SampleRate)
                    )
    fprintf(stderr, "Actual TX Rate: %f...\n\n", SampleRate);

    // Set gain
    fprintf(stderr, "Setting TX Gain: %f db...\n", TunerGain);
    EXECUTE_OR_GOTO(free_tx_metadata,
                    uhd_usrp_set_tx_gain(usrp, TunerGain, 0, "")
                    )

    // See what gain actually is
    EXECUTE_OR_GOTO(free_tx_metadata,
                    uhd_usrp_get_tx_gain(usrp, channel, "", &TunerGain)
                    )
    fprintf(stderr, "Actual TX Gain: %f...\n", TunerGain);

    // Set frequency
    fprintf(stderr, "Setting TX frequency: %f MHz...\n", CenterFrequency / 1e6);
    EXECUTE_OR_GOTO(free_tx_metadata,
                    uhd_usrp_set_tx_freq(usrp, &tune_request, channel, &tune_result)
                    )

    // See what frequency actually is
    EXECUTE_OR_GOTO(free_tx_metadata,
                    uhd_usrp_get_tx_freq(usrp, channel, &CenterFrequency)
                    )
    fprintf(stderr, "Actual TX frequency: %f MHz...\n", CenterFrequency / 1e6);

    // Set up streamer
    stream_args.channel_list = &channel;
    EXECUTE_OR_GOTO(free_tx_streamer,
                    uhd_usrp_get_tx_stream(usrp, &stream_args, tx_streamer)
                    )

    // Set up buffer for idle and busy state
    EXECUTE_OR_GOTO(free_tx_streamer,
                    uhd_tx_streamer_max_num_samps(tx_streamer, &samps_per_buff)
                    )
    fprintf(stderr, "Buffer size in samples: %zu...\n", samps_per_buff);
    buff_busy = malloc(samps_per_buff * 2 * sizeof(float));
    buffs_ptr_busy = (const void**)&buff_busy;
    size_t j = 0;
    for(j = 0; j < (samps_per_buff*2); j+=2){
        buff_busy[j]   = 1.0f;//0.1f
        buff_busy[j+1] = 0;
    }

    // Transmit the sequence of periods
    printf("\n\nTransmitting the sequence of periods ... ");
    fflush(stdout);

    // Ctrl+C will exit loop
    signal(SIGINT, &sigint_handler);
    fprintf(stderr, "\n(Press Ctrl+C to stop streaming)\n");
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // It is better to duplicate the code with and without showing the period
    // durations, rather than checking in every period whether the duration
    // should be shown in the screen. Printing the durations of the transmitted
    // periods in the screen will lead to some additional delay that will add up
    // to the duration of the generated periods. At short time scales this extra
    // delay will introduce noticeable inaccuracies. If the user has opted to show
    // the periods in the screen they know that this will probably lead to inaccurate
    // period durations (this has been warned in the config.txt file). If the user has
    // opted not to show the periods (for the final experiment) then this part of the
    // code becomes a real-time critical part and it is better to only change the state
    // of the pins and generate the delay (without additional checks in every period).

    tm_idle_periods = (struct timespec *) malloc(nofPeriods * sizeof(struct timespec));
    if(tm_idle_periods == NULL)
    {
        printf("\nError: Failed to allocate memory for idle period structs\n\n");
        exit(EXIT_FAILURE);
    }
    for(i=0; i<nofPeriods; i++)
    {
        tm_idle_periods[i].tv_sec = (time_t) idle_periods[i];
        tm_idle_periods[i].tv_nsec = (idle_periods[i] - tm_idle_periods[i].tv_sec)*1e9;
    }

    size_t num_samps_sent = 0;
    double buffer_time=samps_per_buff/SampleRate; // 2040/56e6=36.42e-6 s
    double buffer_time_acc=0;

    if(ShowTransmittedPeriods)
    {
        clock_gettime(CLOCK_MONOTONIC, &TransmissionStartTime);

        for(i=0; i<nofPeriods; i++)
        {
            // Transmit IDLE period
            printf("\nTransmitting IDLE period (TX USRP %s)... %llu/%llu (%.1f%%) = %lf seconds\n", USRPdevice, i + 1, nofPeriods, 100.0*(((double) (i + 1))/((double) nofPeriods)), idle_periods[i]);
            fflush(stdout);
            nanosleep(&tm_idle_periods[i], NULL);

            // Transmit BUSY period
            printf("Transmitting BUSY period (TX USRP %s)... %llu/%llu (%.1f%%) = %lf seconds\n", USRPdevice, i + 1, nofPeriods, 100.0*(((double) (i + 1))/((double) nofPeriods)), busy_periods[i]);
            fflush(stdout);
            buffer_time_acc=0;
            while(buffer_time_acc < busy_periods[i]) {
                if (stop_signal_called) break;
                EXECUTE_OR_GOTO(free_tx_streamer,
                                uhd_tx_streamer_send(tx_streamer, buffs_ptr_busy, samps_per_buff, &md, 0.1, &num_samps_sent)
                                )
                buffer_time_acc += buffer_time;
            }
        }

        clock_gettime(CLOCK_MONOTONIC, &TransmissionEndTime);
    }
    else
    {
        clock_gettime(CLOCK_MONOTONIC, &TransmissionStartTime);

        for(i=0; i<nofPeriods; i++)
        {
            // Transmit IDLE period
            nanosleep(&tm_idle_periods[i], NULL);

            // Transmit BUSY period
            buffer_time_acc=0;
            while(buffer_time_acc < busy_periods[i]) {
                if (stop_signal_called) break;
                EXECUTE_OR_GOTO(free_tx_streamer,
                                uhd_tx_streamer_send(tx_streamer, buffs_ptr_busy, samps_per_buff, &md, 0.1, &num_samps_sent)
                                )
                buffer_time_acc += buffer_time;
            }
        }

        clock_gettime(CLOCK_MONOTONIC, &TransmissionEndTime);
    }

free_tx_streamer:
    uhd_tx_streamer_free(&tx_streamer);

free_tx_metadata:
    uhd_tx_metadata_free(&md);

free_usrp:
    if(*return_code != EXIT_SUCCESS && usrp != NULL){
        uhd_usrp_last_error(usrp, error_string, 512);
        fprintf(stderr, "USRP reported the following error: %s\n", error_string);
    }
    uhd_usrp_free(&usrp);

free_option_strings:
    fprintf(stderr, (*return_code ? "**Failure!\n" : "\nSuccess!\n"));
    if (*return_code == EXIT_FAILURE)
        return;


    // Check real-time operation
    printf("\n\nReal-time execution check:\n");
    expected_total_seconds = 0.0;
    for(i=0; i<nofPeriods; i++)
    {
        // Update total transmission time counter
        expected_total_seconds += idle_periods[i] + busy_periods[i];
    }
    hours = floor(expected_total_seconds/3600.0);
    minutes = floor(((unsigned long int) expected_total_seconds % 3600)/60.0);
    seconds = ((unsigned long int) expected_total_seconds % 3600) % 60;
    printf("\n  Expected transmission duration [HH:MM:SS] = %.2d:%.2d:%.2d (%.6f seconds)", (int) hours, (int) minutes, (int) seconds, expected_total_seconds);
    actual_total_seconds = elapsed_time_ms(TransmissionStartTime, TransmissionEndTime)/1e3;
    hours = floor(actual_total_seconds/3600.0);
    minutes = floor(((unsigned long int) actual_total_seconds % 3600)/60.0);
    seconds = ((unsigned long int) actual_total_seconds % 3600) % 60;
    printf("\n  Actual   transmission duration [HH:MM:SS] = %.2d:%.2d:%.2d (%.6f seconds)", (int) hours, (int) minutes, (int) seconds, actual_total_seconds);
    printf("\n  Actual - Expected transmission time       = %.6f ms", 1000.0 * (actual_total_seconds - expected_total_seconds));
}

double elapsed_time_ms(struct timespec start_time, struct timespec end_time)
{
    return ((double) (end_time.tv_sec  - start_time.tv_sec) )*1e3 + ((double) (end_time.tv_nsec - start_time.tv_nsec))/1e6;
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
    fprintf(stderr,"\nDo you want to continue [y/n]? ");
    //printf("\nDo you want to continue [y/n]? ");
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

/*
 double erf_inv_maclaurin(double x)
 {
 // The approximation used here to compute the inverse error function
 // is based on its Maclaurin series. The accuracy provided by this
 // approximation is good in the interval ~[0.1,0.9] and degrades for
 // smaller/larger arguments out of this interval. In general, the
 // provided accuracy could be sufficient in practice.

 // [OEIS A092676] Numerators of coefficients in the series for inverf(2x/sqrt(Pi))
 double numerators[15] = {1.0, 1.0, 7.0, 127.0, 4369.0, 34807.0, 20036983.0, 2280356863.0, 49020204823.0, 65967241200001.0, 15773461423793767.0, 655889589032992201.0, 94020690191035873697.0, 655782249799531714375489.0, 44737200694996264619809969.0};

 // [OEIS A132467] Denominators associated with Taylor series expansion of inverse error function
 double denominators[15] = {1.0, 1.0, 6.0, 90.0, 2520.0, 16200.0, 7484400.0, 681080400.0, 11675664000.0, 12504636144000.0, 2375880867360000.0, 78404068622880000.0, 8910391798788480000.0, 49229914688306352000000.0, 2658415393168543008000000.0};

 int k;
 double c_k;
 double result = 0.0;

 // Compute the Maclaurin series approximation
 for(k=0; k<15; k++)
 {
 c_k = numerators[k]/denominators[k];
 result = result + ((c_k/(2.0*k+1))*pow((sqrt(M_PI)/2.0)*x, 2.0*k+1));
 }

 // Return the result
 return result;
 }
 */

double erf_inv(double y)
{
    // The funcion will be arbitrarily accurate within the interval
    // from erf(x_min) to erf(x_max)
    double x_min = -5.0, x_max = +5.0;
    double x_mid = (x_min + x_max)/2.0;

    // Tolerance (maximum error)
    double tolerance = 1e-10;

    // Check the input argument
    if(y >= +1.0)
    {
        printf("\nError: erf_inv() input argument y >= 1\n\n");
        exit(EXIT_FAILURE);
    }
    if(y <= -1.0)
    {
        printf("\nError: erf_inv() input argument y <= -1\n\n");
        exit(EXIT_FAILURE);
    }

    // When x is out of the considered interval
    if(y >= erf(x_max))
    {
        printf("\nError: erf_inv() exceeded the maximum value (increase x_max)\n\n");
        exit(EXIT_FAILURE);
        // return x_max;
    }
    if(y <= erf(x_min))
    {
        printf("\nError: erf_inv() exceeded the minimum value (decrease x_min)\n\n");
        exit(EXIT_FAILURE);
        // return x_min;
    }

    // When x is within the considered interval
    while(fabs(y - erf(x_mid)) > tolerance)
    {
        // Adjust search interval as appropriate
        if(y > erf(x_mid))
        {
            x_min = x_mid;
        }
        else
        {
            x_max = x_mid;
        }
        // Recompute middle point
        x_mid = (x_min + x_max)/2.0;
    }

    // Return the result
    return x_mid;
}

double lower_gammainc(double a, double x)
{
    // The approximation used here to compute the lower incomplete gamma
    // function is based on expressions 6.5.4 and 6.5.29 from the
    // "Handbook of Mathematical Functions" by Abramowitz and Stegun.

    int n, nof_terms = 30; // Number of terms to compute the series approximation
    double result = 0.0;

    // Check the input arguments
    if(a <= 0.0)
    {
        printf("\nError: lower_gammainc() input argument a <= 0\n\n");
        exit(EXIT_FAILURE);
    }
    if(x < 0.0)
    {
        printf("\nError: lower_gammainc() input argument x < 0\n\n");
        exit(EXIT_FAILURE);
    }

    // Compute the series approximation
    for(n=0; n<nof_terms; n++)
    {
        result = result + pow(x, (double) n) / tgamma(a + (double) n + 1.0);
    }

    // Return the result
    return pow(x,a)*tgamma(a)*exp(-x)*result;
}

double lower_gammainc_inv(double a, double y)
{
    // The funcion will be arbitrarily accurate within the interval
    // from lower_gammainc(a, x_min) to lower_gammainc(a, x_max)
    double x_min = 0.0, x_max = +3.0;
    double x_mid = (x_min + x_max)/2.0;

    // Tolerance (maximum error)
    double tolerance = 1e-10;

    // Check the input arguments
    if(a <= 0.0)
    {
        printf("\nError: lower_gammainc_inv() input argument a <= 0\n\n");
        exit(EXIT_FAILURE);
    }
    if(y < 0.0)
    {
        printf("\nError: lower_gammainc_inv() input argument y < 0\n\n");
        exit(EXIT_FAILURE);
    }

    // When x is out of the considered interval
    if(y >= lower_gammainc(a, x_max))
    {
        // An automatic readjustement of x_max is preferred since
        // there is no single x_max value adequate for all cases.
        // However, this might hang up if the lower incomplete gamma
        // function is not computed with sufficient accuracy (nof_terms)
        printf("\nWarning: lower_gammainc_inv() exceeded the maximum value. Readjusting x_max.");
        printf("\nIf the program hangs up, increase the value of nof_terms in lower_gammainc().");
        printf("\nNote: The generation of gamma random periods is slower than other distributions (please be patient before deciding to stop the program).\n");
        while(y >= lower_gammainc(a, x_max))
        {
            x_max += 1.0;
        }
    }
    if(y <= lower_gammainc(a, x_min))
    {
        printf("\nError: lower_gammainc_inv() exceeded the minimum value\n\n");
        exit(EXIT_FAILURE);
        // return x_min;
    }

    // When x is within the considered interval
    while(fabs(y - lower_gammainc(a, x_mid)) > tolerance)
    {
        // Adjust search interval as appropriate
        if(y > lower_gammainc(a, x_mid))
        {
            x_min = x_mid;
        }
        else
        {
            x_max = x_mid;
        }
        // Recompute middle point
        x_mid = (x_min + x_max)/2.0;
    }

    // Return the result
    return x_mid;
}

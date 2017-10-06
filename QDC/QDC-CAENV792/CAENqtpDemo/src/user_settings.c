/////////////////////////////////////////////////////////////////////////////
// Name:        user_settings.c
// Purpose:     user setting loading and storing
// Author:      NDA
// Modified by:
// Created:     14/05/07
// RCS-ID:
// Copyright:   CAEN S.p.A. All rights reserved
// Licence:
/////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////
// File includes
////////////////////////////////////////////
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#ifdef WIN32
#include <time.h>
#include <sys/timeb.h>
#include <conio.h>
#include <process.h>
#define popen  _popen   // redefine 'deprecated' popen as _popen
#define pclose _pclose  // redefine 'deprecated' pclose as _pclose
#else
#include <unistd.h>
#include <sys/time.h>
#define Sleep(t)  usleep((t)*1000)
#include <memory.h>
#endif
#include "common_defs.h"
#include "user_settings.h"

////////////////////////////////////////////
// File local defines
////////////////////////////////////////////
#define DEF_RAW_OUTPUT_FILENAME					"raw_values.dat"					/*!< Raw board values default filename */
#define DEF_PARSED_OUTPUT_FILENAME				"parsed_values.txt"					/*!< Parsed board values default filename */
#define DEF_ANALYSIS_OUTPUT_FILENAME				"analysis_values.dat"					/*!< Parsed board values default filename */

#define DEF_NUM_EVENTS						-1							/*!< Default number of events */
#define DEF_BASE_ADDRESS					0x0011							/*!< Default boards' base address */
#define DEF_BOARD_TYPE_STR					"V792A"							/*!< Default V792's type string */
#define DEF_V792_TYPE						CVT_V792_TYPE_A						/*!< Default V792's type */
#define DEF_ZERO_SUPPRESSION_PARAMS				"0"							/*!< Default zero suppression thresholds */
#define DEF_CONFIG_FILENAME					"CAENQDC.conf"					/*!< Configuration file default filename */
#define DEF_GNU_PLOT_PATH					""							/*!< GnuPlot folder default path */
#define DEF_GNU_PLOT_REFRESH					500							/*!< Gnu plot refresh rate (msec) */
#define DEF_GNU_PLOT_X_SCALE					0.3f							/*!< GnuPlot X conversion factor */

////////////////////////////////////////////
// File local variables declaration
////////////////////////////////////////////

////////////////////////////////////////////
// File local methods declaration
////////////////////////////////////////////
//static BOOL parse_acquisition_mode_param( const char* line, acquisition_mode_param *target);
static BOOL parse_zero_suppression_thresholds( const char* line, UINT16 *thresholds_buff);
static BOOL parse_config_file( user_setting_data* p_data);
static char* trim( char charset[], int n_chars, const char* string);

////////////////////////////////////////////
// Global visible variables declaration
////////////////////////////////////////////

/**************************************************
**************************************************/
BOOL user_settings_open( user_setting_data* p_data) {
    // module variable init
    p_data->m_raw_output_filename= (char*)malloc( MAX_FILENAME_LENGHT);
    p_data->m_parsed_output_filename= (char*)malloc( MAX_FILENAME_LENGHT);
    p_data->m_analysis_output_filename= (char*)malloc( MAX_FILENAME_LENGHT);

    p_data->m_gnu_plot_path= (char*)malloc( MAX_FILENAME_LENGHT);
    p_data->m_config_filename= (char*)malloc( MAX_FILENAME_LENGHT);

    p_data->m_gnu_plot_refresh= DEF_GNU_PLOT_REFRESH;
    p_data->m_gnu_plot_x_scale= DEF_GNU_PLOT_X_SCALE;

    // setting default values ...
	if ( p_data->m_raw_output_filename== NULL) {
        TRACE("user setting module: insufficient memory !\n");
	} else {
        strcpy( p_data->m_raw_output_filename, DEF_RAW_OUTPUT_FILENAME);
	}

	if ( p_data->m_parsed_output_filename== NULL) {
        TRACE("user setting module: insufficient memory !\n");
	} else {
        strcpy( p_data->m_parsed_output_filename, DEF_PARSED_OUTPUT_FILENAME);
	}

	if ( p_data->m_analysis_output_filename== NULL) {
        TRACE("user setting module: insufficient memory !\n");
	} else {
        strcpy( p_data->m_analysis_output_filename, DEF_ANALYSIS_OUTPUT_FILENAME);
	}

	if ( p_data->m_gnu_plot_path== NULL) {
        TRACE("user setting module: insufficient memory !\n");
	} else {
        strcpy( p_data->m_gnu_plot_path, DEF_GNU_PLOT_PATH);
	}

	if ( p_data->m_config_filename== NULL) {
        TRACE("user setting module: insufficient memory !\n");
	} else {
        strcpy( p_data->m_config_filename, DEF_CONFIG_FILENAME);
	}

    p_data->m_num_events= DEF_NUM_EVENTS;
    p_data->m_base_address= DEF_BASE_ADDRESS;
    p_data->m_qtp_type= DEF_V792_TYPE;

	parse_zero_suppression_thresholds( DEF_ZERO_SUPPRESSION_PARAMS, p_data->m_zero_suppression_param.m_thresholds_buff);
	parse_zero_suppression_thresholds( DEF_ZERO_SUPPRESSION_PARAMS, p_data->m_zero_suppression_param.m_high_thresholds_buff);

    p_data->m_acquisition_mode_param.m_count_all_events= FALSE;
    p_data->m_acquisition_mode_param.m_empty_enable= FALSE;
    p_data->m_acquisition_mode_param.m_overflow_suppression_enable= FALSE;
    p_data->m_acquisition_mode_param.m_sliding_scale_enable= FALSE;
    p_data->m_acquisition_mode_param.m_zero_suppression_enable= FALSE;
    p_data->m_acquisition_mode_param.m_valid_suppression_enable= FALSE;
    p_data->m_acquisition_mode_param.m_common_stop_enable= FALSE;

    return TRUE;
}
/**************************************************
**************************************************/
BOOL user_settings_close( user_setting_data* p_data) {
    // module variable free
    if ( p_data->m_raw_output_filename)
        free( p_data->m_raw_output_filename);
    if ( p_data->m_parsed_output_filename)
        free( p_data->m_parsed_output_filename);
    if ( p_data->m_analysis_output_filename)
        free( p_data->m_analysis_output_filename);

    // setting invalid values ...
    p_data->m_raw_output_filename= NULL;
    p_data->m_parsed_output_filename= NULL;
    p_data->m_analysis_output_filename= NULL;

    return TRUE;
}

/**************************************************
**************************************************/
BOOL user_settings_parse_input_param( user_setting_data* p_data, int argc, char* argv[]) {

    BOOL help_to_show = FALSE;

    int i;
    for ( i= 1; i< argc; i++) {
        if ( !strcmp( argv[i], "-h")||
                !strcmp( argv[i], "-H")) {
            help_to_show= TRUE;
        } else if ( !strncmp( argv[i], "-f", 2)||
                    !strncmp( argv[i], "-F", 2)) {
            // config filename
            strcpy( p_data->m_config_filename, argv[i]+ 2);
        } else {
            // Unknown parameter
            TRACE1("\nuser settings: ignored command line parameter '%s'", argv[i]);
            TRACE("\n  restart with -h flag for help on usage\n");
            help_to_show= TRUE;

        }
    }
    if ( help_to_show) {
        // show help screen
        TRACE ( "\n");
        TRACE ( "CAEN qtp Demo usage:    CAENqtpDemo {params_list}\n");
        TRACE ( "params_list:\n");
        TRACE ( "-ffilename | -Ffilename\n");
        TRACE1( "      filename= configuration filename. If not specified '%s' will be assumed.\n", DEF_CONFIG_FILENAME);
        TRACE ( "           example: CAENqtpDemo -f\"CAENqtpDemo.conf\"\n");
        TRACE ( "                    CAENqtpDemo -FCAENqtpDemo.conf\n");
        TRACE ( "-h | -H\n");
        TRACE ( "      Show this help screen\n");
        TRACE ( "           example: CAENqtpDemo -h\n");
        TRACE ( "\n");
        TRACE ( "Hit a key to continue ... \n");
        while ( !kbhit());
        return FALSE;
    }
    if ( !parse_config_file( p_data))
        return FALSE;
    return TRUE;
}
/**************************************************
**************************************************/
BOOL parse_config_file( user_setting_data* p_data) {
    // char line[ 200];
    char org_line[ 200];
    char *line;
    char org_str[200];
    char *str;
    FILE *p_conf_file= NULL;
    BOOL ret_val= TRUE;
    char trim_space_charset[]= {' ', '\t', '\n'};
    char trim_text_charset[]= {' ', '\t', '"', '\n'};


    if ( ( p_conf_file= fopen( p_data->m_config_filename, "r"))== NULL) {
        TRACE1(" Cannot open configuration file %s\n", p_data->m_config_filename);
        return FALSE;
    }
    TRACE1(" Reading configuration file '%s' ...", p_data->m_config_filename);
    while (!feof( p_conf_file)) {
        fgets( org_line, sizeof( org_line), p_conf_file);  // read the whole line
	line= trim( trim_space_charset, sizeof( trim_space_charset), org_line);;
        if ( !strlen( line)) {
            continue;
        }
	sscanf( line, "%s", org_str);  // read one string from the file
        str= trim( trim_space_charset, sizeof( trim_space_charset), org_str);
	line+= strlen( str);
        if ( !strlen( str)) {
            continue;
        }
        if ( str[0] == '#') {
            continue;
        }

        // The board base address (16 bytes MSW hex format)
        if ( !stricmp(str, "BASE_ADDRESS")) {
            int value= 0;
            sscanf( line, "%x", &value);
            p_data->m_base_address= value;
            continue;
        }

        // The link number (dec format)
        if ( !stricmp(str, "LINK_NUMBER")) {
            int value= 0;
            sscanf( line, "%d", &value);
			p_data->m_link_number= value;
            continue;
        }

        // The board number (dec format)
        if ( !stricmp(str, "BOARD_NUMBER")) {
            int value= 0;
            sscanf( line, "%d", &value);
			p_data->m_board_number= value;
            continue;
        }

        // The QTP's type: Valid values are: V792A, V792N
        if ( !stricmp(str, "QTP_TYPE")) {
            char value[ 100];
            sscanf( line, "%s", value);
            if ( !stricmp( value, "v792a")) {
                p_data->m_qtp_type= CVT_V792_TYPE_A;
            } else if ( !stricmp( value, "v792n")) {
                p_data->m_qtp_type= CVT_V792_TYPE_N;
            } else if ( !stricmp( value, "v862a")) {
                p_data->m_qtp_type= CVT_V862_TYPE_A;
            } else if ( !stricmp( value, "v965")) {
                p_data->m_qtp_type= CVT_V965;
            } else if ( !stricmp( value, "v965a")) {
                p_data->m_qtp_type= CVT_V965_TYPE_A;
            } else if ( !stricmp( value, "v785a")) {
                p_data->m_qtp_type= CVT_V785_TYPE_A;
            } else if ( !stricmp( value, "v785n")) {
                p_data->m_qtp_type= CVT_V785_TYPE_N;
            } else if ( !stricmp( value, "v1785")) {
                p_data->m_qtp_type= CVT_V1785;
            } else if ( !stricmp( value, "v775a")) {
                p_data->m_qtp_type= CVT_V775_TYPE_A;
            } else if ( !stricmp( value, "v775n")) {
                p_data->m_qtp_type= CVT_V775_TYPE_N;
            } else {
                // Bad parameter
                TRACE1("\nuser settings: bad QTP type '%s'", value);
            }
            continue;
        }
        // The number of events to acquire:
        // A non positive number means don't care (i.e. acquire until key pressed)
        if ( !stricmp(str, "NUM_EVENTS")) {
            long value= 0;
            sscanf( line, "%ld", &value);
            p_data->m_num_events= value;
            continue;
        }
        // The zero suppression threshold resolution
        if ( !stricmp(str, "STEP_THRESHOLD")) {
            int value= 0;
            sscanf( line, "%d", &value);
            p_data->m_zero_suppression_param.m_step_threshold= value;
            continue;
        }
        // The thresholds buffer (value and kill bit) (comma separated values)
		// for dual range boards (V965, V1785) these are the low range thresholds
        if ( !stricmp(str, "THRESHOLDS")) {
            char value[ 400];
            sscanf( line, "%s", value);
			parse_zero_suppression_thresholds( value, p_data->m_zero_suppression_param.m_thresholds_buff);
            continue;
        }
        // The high range thresholds buffer (value and kill bit) (comma separated values) just for dual range boards (V965, V1785)
        if ( !stricmp(str, "THRESHOLDS_HIGH")) {
            char value[ 100];
            sscanf( line, "%s", value);
			parse_zero_suppression_thresholds( value, p_data->m_zero_suppression_param.m_high_thresholds_buff);
            continue;
        }
        // Enable/Disable the sliding scale feature ( set a non zero value to enable)
        if ( !stricmp(str, "SLIDING_SCALE_ENABLE")) {
            int value= 0;
            sscanf( line, "%d", &value);
            p_data->m_acquisition_mode_param.m_sliding_scale_enable= value!= 0;
            continue;
        }
        // Enable/Disable the zero suppression ( set a non zero value to enable)
        if ( !stricmp(str, "ZERO_SUPPRESSION_ENABLE")) {
            int value= 0;
            sscanf( line, "%d", &value);
            p_data->m_acquisition_mode_param.m_zero_suppression_enable= value!= 0;
            continue;
        }
        // Enable/Disable the overflow suppression feature ( set a non zero value to enable)
        if ( !stricmp(str, "OVERFLOW_SUPPRESSION_ENABLE")) {
            int value= 0;
            sscanf( line, "%d", &value);
            p_data->m_acquisition_mode_param.m_overflow_suppression_enable= value!= 0;
            continue;
        }
        // Enable/Disable the valid suppression feature ( set a non zero value to enable)
	// V775 ONLY!!!!
        if ( !stricmp(str, "VALID_SUPPRESSION_ENABLE")) {
            int value= 0;
            sscanf( line, "%d", &value);
            p_data->m_acquisition_mode_param.m_valid_suppression_enable= value!= 0;
            continue;
        }
        // Enable/Disable the common stop acquisition mode ( set a non zero value to enable)
	// V775 ONLY!!!!
        if ( !stricmp(str, "COMMON_STOP_ENABLE")) {
            int value= 0;
            sscanf( line, "%d", &value);
            p_data->m_acquisition_mode_param.m_common_stop_enable= value!= 0;
            continue;
        }
        // Enable/Disable empty event storing feature ( set a non zero value to enable)
        if ( !stricmp(str, "EMPTY_ENABLE")) {
            int value= 0;
            sscanf( line, "%d", &value);
            p_data->m_acquisition_mode_param.m_empty_enable= value!= 0;
            continue;
        }
        // Enable/Disable counting all triggers or accepted one only ( set a non zero value to enable)
        if ( !stricmp(str, "COUNT_ALL_EVENTS")) {
            int value= 0;
            sscanf( line, "%d", &value);
            p_data->m_acquisition_mode_param.m_count_all_events= value!= 0;
            continue;
        }
        // enable debugging
        if (!stricmp(str, "DEBUG")) {
            int value= 0;
            sscanf(line, "%d", &value);
            p_data->m_debug= value!= 0;
            continue;
        }
        // The path of the for where the gnuplot executable is installed
        if (!stricmp(str, "GNUPLOT_PATH")) {
            char value[ MAX_FILENAME_LENGHT];
			value[0]= '\0';
            sscanf( line, "%s", value);
            strcpy( p_data->m_gnu_plot_path, trim( trim_text_charset, sizeof( trim_text_charset), value));
            continue;
        }
        // Gnu plot refresh rate (msec)
        if (!stricmp(str, "GNU_PLOT_REFRESH")) {
            int value= 0;
            sscanf( line, "%d", &value);
            p_data->m_gnu_plot_refresh= value;
            continue;
        }
        // GnuPlot X conversion factor
        if (!stricmp(str, "GNU_PLOT_X_SCALE")) {
            float value= 0;
            sscanf(line, "%f", &value);
            p_data->m_gnu_plot_x_scale= value;
            continue;
        }
        // Raw readout data log filename: set to empty string to disable log
        if (!stricmp(str, "RAW_LOG_FILENAME")) {
            char value[ MAX_FILENAME_LENGHT];
			value[0]= '\0';
            sscanf( line, "%s", value);
            strcpy( p_data->m_raw_output_filename, trim( trim_text_charset, sizeof( trim_text_charset), value));
            continue;
        }
        // Parsed readout data log filename: set to empty string to disable log
        if (!stricmp(str, "PARSED_LOG_FILENAME")) {
            char value[ MAX_FILENAME_LENGHT];
			value[0]= '\0';
            sscanf( line, "%s", value);
            strcpy( p_data->m_parsed_output_filename, trim( trim_text_charset, sizeof( trim_text_charset), value));
            continue;
        }
        if (!stricmp(str, "ANALYSIS_OUTPUT_FILENAME")) {
            char value[ MAX_FILENAME_LENGHT];
			value[0]= '\0';
            sscanf( line, "%s", value);
            strcpy( p_data->m_analysis_output_filename, trim( trim_text_charset, sizeof( trim_text_charset), value));
            continue;
        }
    }
    TRACE(  " Done \n");
    fclose( p_conf_file);
    return ret_val;
}

char* trim( char charset[], int n_chars, const char* string) {
    int i;
    char* str= (char* )string;
    BOOL found_char;
    // Trim leading chars
    found_char= TRUE;
    while ( found_char&& strlen( str) ) {
        found_char= FALSE;
        for ( i= 0; i< n_chars; i++) {
            if ( str[0]== charset[ i]) {
                found_char= TRUE;
                ++str;
                break;
            }
        }
    }
    // Trim trailing chars
    found_char= TRUE;
    while ( found_char&& strlen( str) ) {
        found_char= FALSE;
        for ( i= 0; i< n_chars; i++) {
            if ( str[ strlen( str)- 1]== charset[ i]) {
                found_char= TRUE;
                str[ strlen( str)- 1]= '\0';
                break;
            }
        }
    }
    return str;
}
BOOL parse_zero_suppression_thresholds( const char* line, UINT16 *thresholds_buff) {
    int i;
    const char* tok;
    char buff[ 400];
	char trim_text_charset[]= {'"'};

    strncpy( buff, trim( trim_text_charset, 400, line), 400);
    memset( thresholds_buff, 0, MAX_THRESHOLD_NUM* sizeof( thresholds_buff[ 0]));
    tok= strtok( buff, ", ");
    for ( i= 0; ( i< MAX_THRESHOLD_NUM)&& tok; i++) {
        int tmp= 0;
        if ( !strnicmp( tok, "0x", 2)) {
            // hex format
            sscanf( tok+ 2, "%x", &tmp);
        } else {
            // decimal format
            tmp= atoi( tok);
        }
        thresholds_buff[ i]= (UINT16)tmp;
        tok= strtok( NULL, ", ");
    }
    return TRUE;
}
BOOL parse_acquisition_mode_param( const char* line, acquisition_mode_param *target) {
    const char* tok;
    char buff[ 300];

    strncpy( buff, line, 300);

    tok= strtok( buff, ", ");
    if ( !tok)
        return FALSE;
    target->m_sliding_scale_enable= atoi( tok)? TRUE: FALSE;

    tok= strtok( NULL, ", ");
    if ( !tok)
        return FALSE;
    target->m_zero_suppression_enable= atoi( tok)? TRUE: FALSE;

    tok= strtok( NULL, ", ");
    if ( !tok)
        return FALSE;
    target->m_overflow_suppression_enable= atoi( tok)? TRUE: FALSE;

    tok= strtok( NULL, ", ");
    if ( !tok)
        return FALSE;
    target->m_empty_enable= atoi( tok)? TRUE: FALSE;

    return TRUE;
}


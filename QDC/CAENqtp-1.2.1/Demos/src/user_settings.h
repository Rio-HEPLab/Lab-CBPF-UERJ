/////////////////////////////////////////////////////////////////////////////
// Name:        user_settings.h
// Purpose:     user settng definitions
// Author:      NDA
// Modified by:
// Created:     14/05/07
// RCS-ID:
// Copyright:   CAEN S.p.A. All rights reserved
// Licence:
/////////////////////////////////////////////////////////////////////////////
#ifndef __USER_SETTINGS_DEF_H
#define __USER_SETTINGS_DEF_H

////////////////////////////////////////////
// File includes
////////////////////////////////////////////
#include "common_defs.h"
#include "../include/CAENVMETool/cvt_board_commons.h"
#include "../include/CAENVMETool/cvt_common_defs.h"
#include "../include/CAENVMETool/cvt_V792.h"

////////////////////////////////////////////
// Global defines
////////////////////////////////////////////
#define MAX_FILENAME_LENGHT					255					/*!< Maximum filenam length */
#define MAX_THRESHOLD_NUM					32					/*!< \brief Maximun Number of threshold registers */

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \struct  acquisition_mode_param
*   \brief   Acquisition mode parameters
*
*            This structure stores the relevant parameter for the acquisition mode settings
*/
////////////////////////////////////////////////////////////////////////////////////////////////
typedef struct {
    BOOL m_sliding_scale_enable;						/*!< Enable/disable sliding scale. */
    BOOL m_zero_suppression_enable;						/*!< Enable/disable zero suppression. */
    BOOL m_overflow_suppression_enable;					        /*!< Enable/disable overflow suppression. */
    BOOL m_empty_enable;							/*!< Enable/disable empty events. */
    BOOL m_count_all_events;							/*!< Enable/disable counting all triggers or accepted one only. */
    BOOL m_valid_suppression_enable;                                            /*!< Enable/disable valid data check (V775 only). */
    BOOL m_common_stop_enable;                                                  /*!< Enable/disable common stop mode (V775 only). */
} acquisition_mode_param;

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \struct  zero_suppression_params
*   \brief   Zero suppression parameters
*
*            This structure stores the relevant parameter for the zero suppression settings
*/
////////////////////////////////////////////////////////////////////////////////////////////////
typedef struct {
    UINT16 m_thresholds_buff[ MAX_THRESHOLD_NUM];				/*!< The zero suppression threshold buffer : for dual range boards (V965, V1785) these are the low range thresholds */
    UINT16 m_high_thresholds_buff[ MAX_THRESHOLD_NUM];			/*!< The high range zero suppression threshold buffer (just for dual range boards (V965, V1785)) */
    BOOL m_step_threshold;										/*!< The zero suppression threshold resolution */
} zero_suppression_param;


////////////////////////////////////////////////////////////////////////////////////////////////
/*! \struct  user_setting_data
*   \brief   User settings data structure
*
*            This structure stores user defined application settings
*/
////////////////////////////////////////////////////////////////////////////////////////////////
typedef struct {
    //
    // parameters passed through configuration file or command line
    char* m_config_filename;									/*!< The configuration filename */
    char* m_parsed_output_filename;								/*!< The name for parsed data output data file */
    char* m_raw_output_filename;								/*!< The name for raw output data file */
    UINT16 m_base_address;										/*!< The board base address */

    long m_num_events;											/*!< Number of events to acquire */
    CVT_QTP_TYPES m_qtp_type;									/*!< The qtp type */
    acquisition_mode_param m_acquisition_mode_param;			/*!< The acquisition mode parameters */
    zero_suppression_param m_zero_suppression_param;			/*!< The zero suppression parameters */
    float m_gnu_plot_x_scale;									/*!< The GnuPlot X conversion factor */
    int m_gnu_plot_refresh;										/*!< Gnu plot refresh rate (msec) */
    char *m_gnu_plot_path;										/*!< The path of the for where the gnuplot executable is installed */
    BOOL m_debug;												/*!< Enable debugging features */

    //
    // parameters NOT parsed from configuration file: just for application wide visibility
    CVBoardTypes m_board_type;							/*!< board Type */
    int m_board_number;									/*!< board number */
    int m_link_number;									/*!< link number */
} user_setting_data;

////////////////////////////////////////////
// Global variables declaration
////////////////////////////////////////////

////////////////////////////////////////////
// Global methods declaration
////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL user_settings_open( user_setting_data* p_data)
*   \brief   init user setting module
*
*   \param   *p_data data buffer pointer
*   \return  TRUE procedure completed correctly
*   \note    call before any other methods
*/
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL user_settings_open( user_setting_data* p_data);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL user_settings_close( user_setting_data* p_data)
*   \brief   closes user setting module (frees allocated resources)
*
*   \param   *p_data data buffer pointer
*   \return  TRUE procedure completed correctly
*   \note    call after any other methods
*/
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL user_settings_close( user_setting_data* p_data);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL user_settings_parse_input_param( user_setting_data* p_data, int argc, char* argv[])
*   \brief   Scan<s user input parameters
*
*   \param   *p_data data buffer pointer
*   \param   argc number of input parameters
*   \param   *argv[] string parameters list
*   \return  TRUE procedure completed correctly
*/
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL user_settings_parse_input_param( user_setting_data* p_data, int argc, char* argv[]);

#endif

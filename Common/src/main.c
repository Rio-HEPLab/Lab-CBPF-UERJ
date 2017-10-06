/////////////////////////////////////////////////////////////////////////////
// Name:        main.c
// Purpose:     
// Author:      NDA
// Modified by: 
// Created:     10/04/07
// RCS-ID:      
// Copyright:   CAEN S.p.A. All rights reserved
// Licence:     
/////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////
// File includes
////////////////////////////////////////////
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <time.h>
#if defined (LINUX)
	#include <memory.h>
	#include <ctype.h>
#endif
#include "common_defs.h"
#include "cvt_board_commons.h"
#include "cvt_common_defs.h"
#include "cvt_V1190.h"
#include "cvt_V792.h"
#include "user_settings_QDC.h"
#include "user_settings_TDC.h"
#include "console.h"

////////////////////////////////////////////
// File local defines - TDC
////////////////////////////////////////////
#define GLB_HDR_STR		"GLB_HDR   - EVT COUNT   : %08x GEO      : %08x \n"
#define GLB_TRL_STR		"GLB_TRL   - STATUS      : %08x WCOUNT   : %08x GEO     : %08x \n"
#define TDC_HDR_STR		" TDC_HDR  - TDC         : %08x EVT ID   : %08x BUNCH ID: %08x \n"
#define TDC_MSR_STR		"  TDC_MSR - TRAILING    : %08x CH       : %08x MEASURE : %08x \n"
#define TDC_ERR_STR		"  TDC_ERR - TDC         : %08x ERR FLAGS: %08x \n"
#define TDC_TRL_STR		" TDC_TRL  - TDC         : %08x EVT ID   : %08x WCOUNT  : %08x \n"
#define TDC_TRG_STR		"  TDC_TRG - TRG TIME TAG: %08x \n"
#define FILLER_STR		"  FILLER  - \n"
#define UNKNOWN_STR		"\n??? UNKNOWN TAG ??? -          READ WORD: %08x \n\n"
#define GLB_HDR_STR_TXT		"%d\n"
#define TDC_CHANNEL_COUNT_STR	"%d %d\n"
////////////////////////////////////////////
// File local defines - QDC
////////////////////////////////////////////
#define QDC_HDR_STR			"\n HDR              - GEO  : %08x CRATE         : %08x CH COUNT : %d \n"
#define QDC_DATUM_STR			"\n  DATUM           - GEO  : %08x CHANNEL       : %d ADC      : %d UN : %08x OV : %08x \n"
#define QDC_NOT_VALID_DATUM_STR      	"\n  NOT VALID DATUM - \n"
#define QDC_EOB_STR			"\n EOB              - GEO  : %08x EVENT COUNTER : %d \n"
#define QDC_UNKNOWN_STR			"\n??? UNKNOWN TAG ??? -          READ WORD     : %08x \n\n"

#define EVENT_COUNTER_STR	"Ev %d\n"
#define QDC_CHANNEL_COUNT_STR	"%3d %5d\n"

////////////////////////////////////////////
// File local variables declaration
////////////////////////////////////////////

////////////////////////////////////////////
// Global visible variables declaration
////////////////////////////////////////////

////////////////////////////////////////////
// File local methods declaration
////////////////////////////////////////////


/**************************************************
**************************************************/

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      int main(int argc, char **argv) 
*   \brief   CAENVMETool demo usage for the V1190 board
*            
*            Setups a V1190 board, reads some events and stores into output file
*   \param   argc number of command line arguments
*   \param   *argv[] command line arguments' list 
*   \return  = 0 procedure completed correctly: < 0 some error occurred
*/
////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) 

{
   // TDC
   int ret_val= 0;					// procedure exit value
   cvt_V1190_data board_data_TDC;			// board data
   user_setting_data_TDC user_setting_TDC;		// user settings
   UINT8 *data_buff_TDC = NULL;				// read data buffer
   UINT32 data_size = 0;
   int32_t vme_handle = -1;				// The CAENVMELib handle
   long read_events = 0;
   const int DATA_BUFF_SIZE_TDC = 1024;			// The data buffer size
   // QDC
   cvt_V792_data board_data_QDC;			// board data
   user_setting_data_QDC user_setting_QDC;		// user settings
   UINT8 *data_buff_QDC = NULL;				// read data buffer
   UINT32 DATA_BUFF_SIZE_QDC = 1024;			// The data buffer size

   FILE* raw_TDC_file = NULL;				// raw output file
   FILE* raw_QDC_file = NULL;				// raw output file
   FILE* parsed_out_file = NULL;				// parsed output file
   FILE* analysis_out_file = NULL;			// parsed output file

   /////////////////////////////////////////
   // Application specific
   /////////////////////////////////////////

   memset( &user_setting_TDC, 0, sizeof( user_setting_TDC ));
   memset( &user_setting_QDC, 0, sizeof( user_setting_QDC ));
   memset( &board_data_TDC, 0, sizeof( board_data_TDC ));
   memset( &board_data_QDC, 0, sizeof( board_data_QDC ));

   //
   // init the console module
   con_init( );

   //
   // print header
   con_printf( "\n");
   con_printf( "-------------------------------------------------------\n");
   con_printf( "-           C A E N    V 1 1 9 0 / V 7 9 2            -\n");
   con_printf( "-                                                     -\n");
   con_printf( "-   www.caen.it                                       -\n");
   con_printf( "-------------------------------------------------------\n");
   con_printf( "\n");

   //
   // Init user setting module
   if( !user_settings_TDC_open( &user_setting_TDC ))
   {
      ret_val= -1;
      goto exit_point;
   }
   if( !user_settings_QDC_open( &user_setting_QDC ))
   {
      ret_val= -1;
      goto exit_point;
   }

   //
   // Input parameter check
   if( !user_settings_TDC_parse_input_param( &user_setting_TDC, argc, (char**)argv))
   {
      ret_val= -2;
      goto exit_point;
   }
   //
   // Input parameter check
   if( !user_settings_QDC_parse_input_param( &user_setting_QDC, argc, (char**)argv))
   {
      ret_val= -2;
      goto exit_point;
   }
   //
   // Vme handle initialization
   while( TRUE)
   {
      if( CAENVME_Init( cvV1718, 0, 0, &vme_handle)== cvSuccess)
      {
	 user_setting_TDC.m_board_type= cvV1718;
	 user_setting_TDC.m_board_number= 0;
	 user_setting_TDC.m_link_number= 0;
	 user_setting_QDC.m_board_type= cvV1718;
	 user_setting_QDC.m_board_number= 0;
	 user_setting_QDC.m_link_number= 0;
	 break;
      }
      if( CAENVME_Init( cvV2718, 0, 0, &vme_handle)== cvSuccess)
      {
	 user_setting_TDC.m_board_type= cvV2718;
	 user_setting_TDC.m_board_number= 0;
	 user_setting_TDC.m_link_number= 0;
	 user_setting_QDC.m_board_type= cvV2718;
	 user_setting_QDC.m_board_number= 0;
	 user_setting_QDC.m_link_number= 0;
	 break;
      }

      TRACE("VME INIT ERROR :  press 'Q' to quit or any other key to retry\n");
      if( toupper( con_getch()) == 'Q')
      {
	 ret_val= -3;
	 goto exit_point;
      }
   }

   /////////////////////////////////////////
   // Library specific
   /////////////////////////////////////////

   //
   // Init V1190 board data
   TRACE(  " Initializing V1190 board data ... ");
   if( !cvt_V1190_open( &board_data_TDC, user_setting_TDC.m_base_address, vme_handle, user_setting_TDC.m_V1190_type))
   {	
      TRACE( "\nError executing cvt_V1190_open \n");
      ret_val= -4;
      goto exit_point;
   }
   TRACE(  " Ok \n");
   // Init V792 board data
   TRACE(  " Initializing QTP board data ... ");
   if ( !cvt_V792_open( &board_data_QDC, user_setting_QDC.m_base_address, vme_handle, user_setting_QDC.m_qtp_type)) {
      TRACE( "\nError executing cvt_V792_open \n");
      ret_val= -4;
      goto exit_point;
   }
   printf(" Base address: %08x",user_setting_TDC.m_base_address);
   printf(" Base address: %08x",user_setting_QDC.m_base_address);
   TRACE(  " Ok \n");

   //
   // Get system information - TDC
   {
      UINT32 tdc_id_buff[ MAX_V1190_TDC_COUNT];
      UINT16 firmware_rev;
      UINT16 micro_firmware_rev;
      UINT16 serial_number;
      int i; 

      TRACE(  " Getting system informations ... ");
      if( !cvt_V1190_get_system_info( &board_data_TDC, &firmware_rev, tdc_id_buff, &micro_firmware_rev, &serial_number))
      {
	 TRACE( "\nError executing cvt_V1190_get_system_info \n");
	 ret_val= -5;
	 goto exit_point;
      }
      TRACE(  " Ok \n\n");

      // Show system infos
      TRACE1( "   Firmware Rev.       : %04x\n", firmware_rev);
      TRACE1( "   Micro Firmware Rev. : %04x\n", micro_firmware_rev);
      TRACE1( "   Serial Number       : %04x\n", serial_number);
      for( i= 0; i < user_setting_TDC.m_tdc_count; i++)
      {
	 TRACE2( "   TDC %d               : %08x\n", i, tdc_id_buff[ i]);
      }
   }
   //
   // Get system information - QDC
   {
      UINT16 firmware_rev;
      UINT8  piggy_back_type;
      UINT16 serial_number;

      TRACE(  " Getting system informations ... ");
      if ( !cvt_V792_get_system_info( &board_data_QDC, &firmware_rev, &piggy_back_type, &serial_number)) {
	 TRACE( "\nError executing cvt_V792_get_system_info \n");
	 ret_val= -5;
	 goto exit_point;
      }
      TRACE(  " Ok \n\n");

      // Show system infos
      TRACE1( "   Firmware Rev.       : %04x\n", firmware_rev);
      TRACE1( "   Piggy back type     : %s\n", cvt_V792_get_piggy_back_str( &board_data_QDC, piggy_back_type));
      TRACE1( "   Serial Number       : %04x\n", serial_number);
      TRACE(  "\n");
   }
   //
   // Data clear - TDC
   TRACE(  " Sending data clear ... ");
   if( !cvt_V1190_data_clear( &board_data_TDC))	
   {	
      TRACE( "\nError executing cvt_V1190_data_clear \n");
      ret_val= -5;
      goto exit_point;
   }
   TRACE(  " Ok \n");
   //
   // Data clear - QDC
   TRACE(  " Sending data clear ... ");
   if ( !cvt_V792_data_clear( &board_data_QDC)) {
      TRACE( "\nError executing cvt_V792_data_clear \n");
      ret_val= -5;
      goto exit_point;
   }
   TRACE(  " Ok \n");

   //
   // Sliding constant (QDC)
   TRACE(  " Setting sliding constant... ");
   if ( !cvt_V792_set_sliding_scale( &board_data_QDC,
	    user_setting_QDC.m_acquisition_mode_param.m_sliding_scale_enable)) {
      TRACE( "Error executing cvt_V792_set_sliding_scale\n");
      ret_val= -5;
      goto exit_point;
   }
   TRACE(  " Ok\n");

   //
   // Zero suppression (QDC)
   TRACE(  " Setting zero suppression... ");
   switch( board_data_QDC.m_type) {
      case CVT_V965:
      case CVT_V965_TYPE_A:
      case CVT_V1785:
	 // Dual range board
	 if ( !cvt_V792_set_dual_range_zero_suppression( &board_data_QDC,
		  user_setting_QDC.m_acquisition_mode_param.m_zero_suppression_enable,
		  user_setting_QDC.m_zero_suppression_param.m_step_threshold,
		  &user_setting_QDC.m_zero_suppression_param.m_high_thresholds_buff[0],
		  &user_setting_QDC.m_zero_suppression_param.m_thresholds_buff[0])) {
	    TRACE( "Error executing cvt_V792_set_dual_range_zero_suppression\n");
	    ret_val= -5;
	    goto exit_point;
	 }
	 break;
      default:
	 if ( !cvt_V792_set_zero_suppression( &board_data_QDC,
		  user_setting_QDC.m_acquisition_mode_param.m_zero_suppression_enable,
		  user_setting_QDC.m_zero_suppression_param.m_step_threshold,
		  &user_setting_QDC.m_zero_suppression_param.m_thresholds_buff[0])) {
	    TRACE( "Error executing cvt_V792_set_zero_suppression\n");
	    ret_val= -5;
	    goto exit_point;
	 }
	 break;
   }
   TRACE(  " Ok\n");

   //
   // Acquisition mode - TDC
   TRACE(  " Setting acquisition mode ... ");
   switch( user_setting_TDC.m_acquisition_mode.m_mode)
   {
      case AM_CONTINUOUS:
	 if( !cvt_V1190_set_continuous_acquisition_mode( &board_data_TDC, 
		  user_setting_TDC.m_acquisition_mode.m_params.m_continuos.m_edge_detection,
		  user_setting_TDC.m_acquisition_mode.m_params.m_continuos.m_res_width,
		  &user_setting_TDC.m_acquisition_mode.m_params.m_continuos.m_enable_msk[0]))
	 {	
	    TRACE( "Error executing cvt_V1190_set_continuous_acquisition_mode \n");
	    ret_val= -5;
	    goto exit_point;
	 }
	 break;
      case AM_TRIGGER_MATCH:
	 TRACE("\n Setting Trigger Matching Mode ...");
	 if( !cvt_V1190_set_trigger_matching_acquisition_mode(	&board_data_TDC, 
		  user_setting_TDC.m_acquisition_mode.m_params.m_trigger_matching.m_window_width,
		  user_setting_TDC.m_acquisition_mode.m_params.m_trigger_matching.m_window_offset,
		  user_setting_TDC.m_acquisition_mode.m_params.m_trigger_matching.m_extra_search_margin,
		  user_setting_TDC.m_acquisition_mode.m_params.m_trigger_matching.m_reject_margin,
		  user_setting_TDC.m_acquisition_mode.m_params.m_trigger_matching.m_edge_detection,
		  user_setting_TDC.m_acquisition_mode.m_params.m_trigger_matching.m_res_width,
		  &user_setting_TDC.m_acquisition_mode.m_params.m_trigger_matching.m_enable_msk[0],
		  user_setting_TDC.m_acquisition_mode.m_params.m_trigger_matching.m_header_trailer_enable,
		  user_setting_TDC.m_acquisition_mode.m_params.m_trigger_matching.m_empty_event_enable,
		  user_setting_TDC.m_acquisition_mode.m_params.m_trigger_matching.m_trigger_time_tag_enable))
	 {	
	    TRACE( "Error executing cvt_V1190_set_trigger_matching_acquisition_mode \n");
	    ret_val= -5;
	    goto exit_point;
	 }
	 break;
      case AM_NONE:
	 break;
   }
   TRACE(  " Ok\n");
   //
   // Acquisition mode - QDC
   TRACE(  " Setting acquisition mode ... ");
   if ( !cvt_V792_set_acquisition_mode( &board_data_QDC,
	    user_setting_QDC.m_acquisition_mode_param.m_sliding_scale_enable,
	    user_setting_QDC.m_acquisition_mode_param.m_zero_suppression_enable,
	    user_setting_QDC.m_acquisition_mode_param.m_overflow_suppression_enable,
	    user_setting_QDC.m_acquisition_mode_param.m_valid_suppression_enable,
	    user_setting_QDC.m_acquisition_mode_param.m_common_stop_enable,
	    user_setting_QDC.m_acquisition_mode_param.m_empty_enable,
	    user_setting_QDC.m_acquisition_mode_param.m_count_all_events
	    )) {
      TRACE( "Error executing cvt_V792_set_continuous_acquisition_mode \n");
      ret_val= -5;
      goto exit_point;
   }
   TRACE(  " Ok\n");

   INT16 tdc_lsb = 3;
   // TDC resolution - implemented HN
   TRACE(" Setting TDC resolution ...");
   if (!cvt_V1190_set_TDC_LSB(&board_data_TDC, tdc_lsb))
   {
      TRACE("Error executing cvt_V1190_set_TDC_LSB \n");
      ret_val = -5;
      goto exit_point;
   }
   TRACE("Ok\n");

   //
   // Readout mode - TDC
   TRACE(  " Setting readout mode (V1190)...");
   if( !cvt_V1190_set_readout_mode( &board_data_TDC, TRUE, TRUE, 1))
   {	
      TRACE( "Error executing cvt_V1190_set_readout_mode \n");
      ret_val= -5;
      goto exit_point;
   }
   TRACE(  " Ok\n");
   //
   // Readout mode - QDC
   TRACE(  " Setting readout mode (V792)...");
   if ( !cvt_V792_set_readout_mode( &board_data_QDC, TRUE, TRUE, TRUE)) {
      TRACE( "Error executing cvt_V792_set_readout_mode \n");
      ret_val= -5;
      goto exit_point;
   }
   TRACE(  " Ok\n");

   // Allocate buffer storage
   data_buff_TDC = malloc( DATA_BUFF_SIZE_TDC );
   if( data_buff_TDC == NULL)
   {
      // Insufficient memory
      TRACE1( "Error allocating events' buffer (%li bytes)", DATA_BUFF_SIZE_TDC );
      ret_val= -5;
      goto exit_point;
   }
   // Allocate buffer storage
   data_buff_QDC = (UINT8*)malloc( DATA_BUFF_SIZE_QDC );
   if ( data_buff_QDC == NULL ) {
      // Insufficient memory
      TRACE1( "Error allocating events' buffer (%i bytes)", DATA_BUFF_SIZE_QDC );
      ret_val= -5;
      goto exit_point;
   }

   // Create output files
   if( ( raw_TDC_file = fopen( user_setting_TDC.m_raw_output_filename, "w+b"))== NULL)
   {
      TRACE1( "Error creating raw output file '%s' \n", user_setting_TDC.m_raw_output_filename);
      ret_val= -5;
      goto exit_point;
   }
   if( ( raw_QDC_file = fopen( user_setting_QDC.m_raw_output_filename, "w+b"))== NULL)
   {
      TRACE1( "Error creating raw output file '%s' \n", user_setting_QDC.m_raw_output_filename);
      ret_val= -5;
      goto exit_point;
   }
   if ( strlen( user_setting_QDC.m_parsed_output_filename)) {
      if ( ( parsed_out_file= fopen( user_setting_QDC.m_parsed_output_filename, "wt"))== NULL) {
	 TRACE1( "Error creating parsed output file '%s'", user_setting_QDC.m_parsed_output_filename);
	 ret_val= -5;
	 goto exit_point;
      }
   }
   if ( strlen( user_setting_QDC.m_analysis_output_filename)) {
      if ( ( analysis_out_file= fopen( user_setting_QDC.m_analysis_output_filename, "wt"))== NULL) {
	 TRACE1( "Error creating analysis output file '%s'", user_setting_QDC.m_analysis_output_filename);
	 ret_val= -5;
	 goto exit_point;
      }
   }

   //
   // Start acquisition 
   time_t time_begin;
   time(&time_begin);
   struct tm * time_info_begin;
   time_info_begin = localtime(&time_begin);

   TRACE (  " Getting events: hit any key to abort ...\n");
   TRACE1(  " Local time %s",asctime(time_info_begin));
   TRACE (  "\n ");

   while( ( (read_events < user_setting_TDC.m_num_events) || 
            (user_setting_TDC.m_num_events <= 0) ) && !con_kbhit() )
   {
      static long last_read_bytes= 0;
      static long read_bytes= 0;
      data_size = DATA_BUFF_SIZE_TDC;
      //
      // Read from MEB .....
      if( !cvt_V1190_read_MEB( &board_data_TDC, data_buff_TDC, &data_size ) )
      {
	 TRACE( " \nError executing cvt_V1190_read_MEB \n");
	 ret_val = -5;
	 goto exit_point;
      }
      if( !data_size )
	 continue;
      //
      // .... and store raw data into file
      if( fwrite( data_buff_TDC, 1, data_size, raw_TDC_file) != data_size)
      {
	 TRACE( " \nError writing raw data file (TDC) \n");
	 ret_val = -5;
	 goto exit_point;
      }
      read_bytes+= data_size;
      if( ( read_bytes >> 10 ) != ( last_read_bytes >> 10 ) )
      {
	 // Give user a feedback every 1KB data
	 TRACE( "." );
      }
      last_read_bytes= read_bytes;
      ++read_events;
      if( read_events % 100 == 0) TRACE1( "\n Analyzed %dth event",read_events); 
   }
   TRACE(  "\n Done \n");

   time_t time_end;
   time(&time_end);
   struct tm * time_info_end;
   time_info_end = localtime(&time_end);

   double time_run = difftime(time_end,time_begin);

   TRACE1(  " Local time %s",asctime(time_info_end));
   TRACE1(  " Total time ellapsed was %.0f seconds.\n",time_run);
   TRACE (  "\n");
   TRACE (  " Parsing events and writing to output file ... ");
/*
	//
	// Post process data : parse raw data and store events in clear text
	if( fflush( raw_out_file))
	{
		TRACE1( "\nError flushing raw output file '%s' \n", user_setting.m_raw_output_filename);
		ret_val= -5;
		goto exit_point;
	}
	if( fseek( raw_out_file, 0, SEEK_SET))
	{
		TRACE1( "\nError flushing raw output file '%s' \n", user_setting.m_raw_output_filename);
		ret_val= -5;
		goto exit_point;
	}

	while( ( data_size= ( UINT32)fread( data_buff, 4, DATA_BUFF_SIZE>> 2, raw_out_file)))
	{
		UINT32 *tmp_buff= (UINT32*)data_buff;
		//char line[ 400];
		char line[400], line_txt[400];
		size_t str_len;
		UINT16 data_line = 0;

		while( data_size-- > 0)
		{
			data_line = 0;
			UINT32 data= *(tmp_buff++);
			*line     = '\0';
			*line_txt = '\0';
			switch( data& CVT_V1190_DATA_TYPE_MSK)
			{
			case CVT_V1190_GLOBAL_HEADER:
				{
					// Global header
					UINT32 event_count= CVT_V1190_GET_GLB_HDR_EVENT_COUNT( data);
					UINT32 geo= CVT_V1190_GET_GLB_HDR_GEO( data);
					data_line = 1;

					sprintf( line, GLB_HDR_STR, event_count, geo);
					sprintf( line_txt, GLB_HDR_STR_TXT, event_count);
				}
				break;
			case CVT_V1190_GLOBAL_TRAILER:
				{
					// Global trailer
					UINT32 status= CVT_V1190_GET_GLB_TRL_STATUS( data);
					UINT32 wcount= CVT_V1190_GET_GLB_TRL_WCOUNT( data);
					UINT32 geo= CVT_V1190_GET_GLB_TRL_GEO( data);

					sprintf( line, GLB_TRL_STR, status, wcount, geo);
				}
				break;
			case CVT_V1190_TDC_HEADER:
				{
					// TDC header
					UINT32 tdc= CVT_V1190_GET_TDC_HDR_TDC( data);
					UINT32 event_id= CVT_V1190_GET_TDC_HDR_EVENT_ID( data);
					UINT32 bunch_id= CVT_V1190_GET_TDC_HDR_BUNCH_ID( data);

					sprintf( line, TDC_HDR_STR, tdc, event_id, bunch_id);
				}
				break;
			case CVT_V1190_TDC_MEASURE:
				{
					// TDC measure
					UINT32 trailing= CVT_V1190_GET_TDC_MSR_TRAILING( data);
					UINT32 channel= CVT_V1190_GET_TDC_MSR_CHANNEL( data);
					UINT32 measure= CVT_V1190_GET_TDC_HDR_MEASURE( data);
					//fix July 2013
					if (user_setting.m_V1190_type ==  CVT_V1190_TYPE_A || user_setting.m_V1190_type == CVT_V1190_TYPE_B) {
						channel= CVT_V1190_GET_TDC_MSR_CHANNEL( data);
						measure= CVT_V1190_GET_TDC_HDR_MEASURE( data);
					}
					else {
						channel= CVT_V1290_GET_TDC_MSR_CHANNEL( data);
						measure= CVT_V1290_GET_TDC_HDR_MEASURE( data);
					}

					//sprintf( line, TDC_MSR_STR, trailing, channel, measure);
					data_line = 1;
					sprintf(line, TDC_MSR_STR, trailing, channel, measure);
					sprintf(line_txt, TDC_MSR_STR_TXT, channel, measure);
				}
				break;
			case CVT_V1190_TDC_ERROR:
				{
					UINT32 tdc= CVT_V1190_GET_TDC_ERR_TDC( data);
					UINT32 err_flags= CVT_V1190_GET_TDC_ERR_ERR_FLAGS( data);

					sprintf( line, TDC_ERR_STR, tdc, err_flags);
				}
				break;
			case CVT_V1190_TDC_TRAILER:
				{
					UINT32 tdc= CVT_V1190_GET_TDC_TRL_TDC( data);
					UINT32 event_id= CVT_V1190_GET_TDC_TRL_EVENT_ID( data);
					UINT32 wcount= CVT_V1190_GET_TDC_TRL_WCOUNT( data);

					sprintf( line, TDC_TRL_STR, tdc, event_id, wcount);
				}
				break;
			case CVT_V1190_GLOBAL_TRIGGER_TIME:
				{
					UINT32 trg_time_tag= CVT_V1190_GET_GLB_TRG_TIME_TAG( data);

					sprintf( line, TDC_TRG_STR, trg_time_tag);
				}
				break;
			case CVT_V1190_FILLER:
				{
					sprintf( line, FILLER_STR);
				}
				break;
			default:
				{
					sprintf( line, UNKNOWN_STR, data );
					sprintf( line_txt, UNKNOWN_STR, data );
				}
				break;
			}
			if( (str_len= strlen( line))> 0)
			{
				if( fwrite( line, 1, str_len, parsed_out_file)!= str_len)
				{
					// error writing file
					TRACE1( "\nError writing parsed output file '%s' \n", user_setting.m_parsed_output_filename);
					ret_val= -5;
					goto exit_point;
				}
			}
			if ((str_len = strlen(line_txt)) > 0)
			{
				if (data_line && (fwrite(line_txt, 1, str_len, parsed_txt_file) != str_len))
				{
					// error writing file
					TRACE1("\nError writing parsed output file '%s' \n", user_setting.m_parsed_txt_filename);
					ret_val = -5;
 					goto exit_point;
 				}
 			}
		}
	}	
	TRACE(  " Ok\n");
*/

exit_point:	
	/////////////////////////////////////////
	// Library specific
	/////////////////////////////////////////
	
	//
	// Release board resources
	if( !cvt_V1190_close( &board_data_TDC ) )
	{	
		TRACE( "\nError executing cvt_V1190_close\n");
	}

	/////////////////////////////////////////
	// Demo application specific
	/////////////////////////////////////////

	if( raw_TDC_file!= NULL)
	{
	   fclose( raw_TDC_file);
	}
	if( raw_QDC_file!= NULL)
	{
	   fclose( raw_QDC_file);
	}
	if( parsed_out_file != NULL )
	{
	   fclose( parsed_out_file );
	}
	if( analysis_out_file != NULL )
	{
	   fclose( analysis_out_file );
	}

	if( data_buff_TDC != NULL)
	{
	   free( data_buff_TDC );
	}
	if( data_buff_QDC != NULL)
	{
	   free( data_buff_QDC );
	}

	// close modules
	user_settings_TDC_close( &user_setting_TDC );
	user_settings_QDC_close( &user_setting_QDC );

	TRACE( " Hit a key to exit ...");
	con_getch();
	//
	// close the console module
	con_end( );

	return ret_val;
}


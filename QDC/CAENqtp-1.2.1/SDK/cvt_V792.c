////////////////////////////////////////////////////////////////////////////////////////////////
/*! \file    cvt_V792.c
*   \brief   V792 VME board implementation
*   \author  NDA
*   \version 1.0
*   \date    03/2006
*            
*            Provides methods, properties and defines to handle V792 VME boards
*/
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
// File includes
////////////////////////////////////////////////////////////////////////////////////////////////
#include <memory.h>
#include <string.h>
#include "cvt_V792.h"

////////////////////////////////////////////////////////////////////////////////////////////////
// File local defines
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
// Global variables declaration
////////////////////////////////////////////////////////////////////////////////////////////////

const int CVT_QTP_NUM_CHANNELS[]= {
	32,							/*!< \brief The board is V792A */
	16,							/*!< \brief The board is V792N */
	32,							/*!< \brief The board is V862A */
	16,							/*!< \brief The board is V965  */
	8,							/*!< \brief The board is V965A */
	32,							/*!< \brief The board is V785A */
	16,							/*!< \brief The board is V785N */
	8,							/*!< \brief The board is V1785 */
	32,							/*!< \brief The board is V775A */
	16,							/*!< \brief The board is V775N */
};

////////////////////////////////////////////////////////////////////////////////////////////////
// Static variables declaration
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \var     static const cvt_reg_table CVT_V792_REG_TABLE[]
*   \brief   The board registers table
*            
*            Provides an entry for each declared register: keep synched with \ref CVT_V792_REG_INDEX
*/
////////////////////////////////////////////////////////////////////////////////////////////////
static const cvt_reg_table CVT_V792_REG_TABLE[]=
{
	{CVT_V792_OUT_BUFFER_ADD,			CVT_V792_OUT_BUFFER_AM,				CVT_V792_OUT_BUFFER_DATA_SIZE},					/*!< \brief Output buffer */
	{CVT_V792_FW_REV_ADD,				CVT_V792_FW_REV_AM,					CVT_V792_FW_REV_DATA_SIZE},						/*!< \brief Firmware Revision register */
	{CVT_V792_GEO_ADDRESS_ADD,			CVT_V792_GEO_ADDRESS_AM,			CVT_V792_GEO_ADDRESS_DATA_SIZE},				/*!< \brief Geo Address register */
	{CVT_V792_MCST_CBLT_ADDRESS_ADD,	CVT_V792_MCST_CBLT_ADDRESS_AM,		CVT_V792_MCST_CBLT_ADDRESS_DATA_SIZE},			/*!< \brief MCST/CBLT Address register */
	{CVT_V792_BIT_SET_1_ADD,			CVT_V792_BIT_SET_1_AM,				CVT_V792_BIT_SET_1_DATA_SIZE},					/*!< \brief Bit Set 1 register */
	{CVT_V792_BIT_CLEAR_1_ADD,			CVT_V792_BIT_CLEAR_1_AM,			CVT_V792_BIT_CLEAR_1_DATA_SIZE},				/*!< \brief Bit Clear 1 register */
	{CVT_V792_INT_LEVEL_ADD,			CVT_V792_INT_LEVEL_AM,				CVT_V792_INT_LEVEL_DATA_SIZE},					/*!< \brief Interrupt level register */
	{CVT_V792_INT_VECTOR_ADD,			CVT_V792_INT_VECTOR_AM,				CVT_V792_INT_VECTOR_DATA_SIZE},					/*!< \brief Interrupt vector register */
	{CVT_V792_STATUS_1_ADD,				CVT_V792_STATUS_1_AM,				CVT_V792_STATUS_1_DATA_SIZE},					/*!< \brief Status 1 register */
	{CVT_V792_CONTROL_1_ADD,			CVT_V792_CONTROL_1_AM,				CVT_V792_CONTROL_1_DATA_SIZE},					/*!< \brief Control 1 register */
	{CVT_V792_ADER_HIGH_ADD,			CVT_V792_ADER_HIGH_AM,				CVT_V792_ADER_HIGH_DATA_SIZE},					/*!< \brief Address decoder HIGH register */
	{CVT_V792_ADER_LOW_ADD,				CVT_V792_ADER_LOW_AM,				CVT_V792_ADER_LOW_DATA_SIZE},					/*!< \brief Address decoder LOW register */
	{CVT_V792_SINGLE_SHOT_RESET_ADD,	CVT_V792_SINGLE_SHOT_RESET_AM,		CVT_V792_SINGLE_SHOT_RESET_DATA_SIZE},			/*!< \brief Single Shot Reset register */
	{CVT_V792_MCST_CBLT_CTRL_ADD,		CVT_V792_MCST_CBLT_CTRL_AM,			CVT_V792_MCST_CBLT_CTRL_DATA_SIZE},				/*!< \brief MCST/CBLT control register */
	{CVT_V792_EVENT_TRG_ADD,			CVT_V792_EVENT_TRG_AM,				CVT_V792_EVENT_TRG_DATA_SIZE},					/*!< \brief Event trigger register */
	{CVT_V792_STATUS_2_ADD,				CVT_V792_STATUS_2_AM,				CVT_V792_STATUS_2_DATA_SIZE},					/*!< \brief Status 2 register */
	{CVT_V792_EVENT_COUNTER_LOW_ADD,	CVT_V792_EVENT_COUNTER_LOW_AM,		CVT_V792_EVENT_COUNTER_LOW_DATA_SIZE},			/*!< \brief Event counter low register */
	{CVT_V792_EVENT_COUNTER_HIGH_ADD,	CVT_V792_EVENT_COUNTER_HIGH_AM,		CVT_V792_EVENT_COUNTER_HIGH_DATA_SIZE},			/*!< \brief Event counter high register */
	{CVT_V792_INC_EVENT_ADD,			CVT_V792_INC_EVENT_AM,				CVT_V792_INC_EVENT_DATA_SIZE},					/*!< \brief Increment Event register */
	{CVT_V792_INC_OFFSET_ADD,			CVT_V792_INC_OFFSET_AM,				CVT_V792_INC_OFFSET_DATA_SIZE},					/*!< \brief Increment Offset register */
	{CVT_V792_LOAD_TEST_ADD,			CVT_V792_LOAD_TEST_AM,				CVT_V792_LOAD_TEST_DATA_SIZE},					/*!< \brief Load Test register */
	{CVT_V792_FCLR_WND_ADD,				CVT_V792_FCLR_WND_AM,				CVT_V792_FCLR_WND_DATA_SIZE},					/*!< \brief FCLR window register */
	{CVT_V792_BIT_SET_2_ADD,			CVT_V792_BIT_SET_2_AM,				CVT_V792_BIT_SET_2_DATA_SIZE},					/*!< \brief Bit Set 2 register */
	{CVT_V792_BIT_CLEAR_2_ADD,			CVT_V792_BIT_CLEAR_2_AM,			CVT_V792_BIT_CLEAR_2_DATA_SIZE},				/*!< \brief Bit Clear 2 register */
	{CVT_V792_W_MEMORY_TEST_ADDRESS_ADD,CVT_V792_W_MEMORY_TEST_ADDRESS_AM,	CVT_V792_W_MEMORY_TEST_ADDRESS_DATA_SIZE},		/*!< \brief W Memory Test Address register */
	{CVT_V792_MEMORY_TEST_WORD_HIGH_ADD,CVT_V792_MEMORY_TEST_WORD_HIGH_AM,	CVT_V792_MEMORY_TEST_WORD_HIGH_DATA_SIZE},		/*!< \brief Memory Test Word High register */
	{CVT_V792_MEMORY_TEST_WORD_LOW_ADD, CVT_V792_MEMORY_TEST_WORD_LOW_AM,	CVT_V792_MEMORY_TEST_WORD_LOW_DATA_SIZE},		/*!< \brief Memory Test Word Low register */
	{CVT_V792_CRATE_SELECT_ADD,			CVT_V792_CRATE_SELECT_AM,			CVT_V792_CRATE_SELECT_DATA_SIZE},				/*!< \brief Crate select register */
	{CVT_V792_TEST_EVENT_WRITE_ADD,		CVT_V792_TEST_EVENT_WRITE_AM,		CVT_V792_TEST_EVENT_WRITE_DATA_SIZE},			/*!< \brief Test Event register Write */
	{CVT_V792_EVENT_COUNT_RESET_ADD,	CVT_V792_EVENT_COUNT_RESET_AM,		CVT_V792_EVENT_COUNT_RESET_DATA_SIZE},			/*!< \brief Event Counter Reset register */
	{CVT_V792_IPED_ADD,					CVT_V792_IPED_AM,					CVT_V792_IPED_DATA_SIZE},						/*!< \brief Iped register */
	{CVT_V792_R_TEST_ADDRESS_ADD,		CVT_V792_R_TEST_ADDRESS_AM,			CVT_V792_R_TEST_ADDRESS_DATA_SIZE},				/*!< \brief R Test Address register */
	{CVT_V792_SW_COMM_ADD,				CVT_V792_SW_COMM_AM,				CVT_V792_SW_COMM_DATA_SIZE},					/*!< \brief SW Comm register */
	{CVT_V792_SLIDE_CONSTANT_ADD,		CVT_V792_SLIDE_CONSTANT_AM,			CVT_V792_SLIDE_CONSTANT_DATA_SIZE},				/*!< \brief Sliding constant register */
	{CVT_V792_AAD_ADD,					CVT_V792_AAD_AM,					CVT_V792_AAD_DATA_SIZE},						/*!< \brief AAD register */
	{CVT_V792_BAD_ADD,					CVT_V792_BAD_AM,					CVT_V792_BAD_DATA_SIZE},						/*!< \brief BAD register */
	{CVT_V792_THRESHOLD_0_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 0 register */
	{CVT_V792_THRESHOLD_1_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 1 register */
	{CVT_V792_THRESHOLD_2_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 2 register */
	{CVT_V792_THRESHOLD_3_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 3 register */
	{CVT_V792_THRESHOLD_4_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 4 register */
	{CVT_V792_THRESHOLD_5_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 5 register */
	{CVT_V792_THRESHOLD_6_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 6 register */
	{CVT_V792_THRESHOLD_7_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 7 register */
	{CVT_V792_THRESHOLD_8_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 8 register */
	{CVT_V792_THRESHOLD_9_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 9 register */
	{CVT_V792_THRESHOLD_10_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 10 register */
	{CVT_V792_THRESHOLD_11_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 11 register */
	{CVT_V792_THRESHOLD_12_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 12 register */
	{CVT_V792_THRESHOLD_13_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 13 register */
	{CVT_V792_THRESHOLD_14_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 14 register */
	{CVT_V792_THRESHOLD_15_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 15 register */
	{CVT_V792_THRESHOLD_16_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 16 register */
	{CVT_V792_THRESHOLD_17_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 17 register */
	{CVT_V792_THRESHOLD_18_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 18 register */
	{CVT_V792_THRESHOLD_19_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 19 register */
	{CVT_V792_THRESHOLD_20_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 20 register */
	{CVT_V792_THRESHOLD_21_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 21 register */
	{CVT_V792_THRESHOLD_22_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 22 register */
	{CVT_V792_THRESHOLD_23_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 23 register */
	{CVT_V792_THRESHOLD_24_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 24 register */
	{CVT_V792_THRESHOLD_25_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 25 register */
	{CVT_V792_THRESHOLD_26_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 26 register */
	{CVT_V792_THRESHOLD_27_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 27 register */
	{CVT_V792_THRESHOLD_28_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 28 register */
	{CVT_V792_THRESHOLD_29_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 29 register */
	{CVT_V792_THRESHOLD_30_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 30 register */
	{CVT_V792_THRESHOLD_31_ADD,			CVT_V792_THRESHOLD_AM,				CVT_V792_THRESHOLD_DATA_SIZE},				/*!< \brief Threshold 31 register */
	{CVT_V792_ROM_OUI_2_ADD,			CVT_V792_ROM_OUI_2_AM,				CVT_V792_ROM_OUI_2_DATA_SIZE},				/*!< \brief Manufacturer identifier (IEEE OUI) (MSB) register*/
	{CVT_V792_ROM_OUI_1_ADD,			CVT_V792_ROM_OUI_1_AM,				CVT_V792_ROM_OUI_1_DATA_SIZE},				/*!< \brief Manufacturer identifier (IEEE OUI) register*/
	{CVT_V792_ROM_OUI_0_ADD,			CVT_V792_ROM_OUI_0_AM,				CVT_V792_ROM_OUI_0_DATA_SIZE},				/*!< \brief Manufacturer identifier (IEEE OUI) (LSB) register*/
	{CVT_V792_ROM_VERSION_ADD,			CVT_V792_ROM_VERSION_AM,			CVT_V792_ROM_VERSION_DATA_SIZE},			/*!< \brief Purchased version of the Mod.V792 register*/
	{CVT_V792_ROM_BOARD_ID_2_ADD,		CVT_V792_ROM_BOARD_ID_2_AM,			CVT_V792_ROM_BOARD_ID_2_DATA_SIZE},			/*!< \brief Board identifier (MSB) register*/
	{CVT_V792_ROM_BOARD_ID_1_ADD,		CVT_V792_ROM_BOARD_ID_1_AM,			CVT_V792_ROM_BOARD_ID_1_DATA_SIZE},			/*!< \brief Board identifier register*/
	{CVT_V792_ROM_BOARD_ID_0_ADD,		CVT_V792_ROM_BOARD_ID_0_AM,		CVT_V792_ROM_BOARD_ID_0_DATA_SIZE},		/*!< \brief Board identifier (LSB) register*/
	{CVT_V792_ROM_REVISION_ADD,			CVT_V792_ROM_REVISION_AM,			CVT_V792_ROM_REVISION_DATA_SIZE},			/*!< \brief Hardware revision identifier register*/
	{CVT_V792_ROM_SERIAL_1_ADD,			CVT_V792_ROM_SERIAL_1_AM,			CVT_V792_ROM_SERIAL_1_DATA_SIZE},			/*!< \brief Serial number (MSB) register*/
	{CVT_V792_ROM_SERIAL_0_ADD,			CVT_V792_ROM_SERIAL_0_AM,			CVT_V792_ROM_SERIAL_0_DATA_SIZE},			/*!< \brief Serial number (LSB) register*/
};

////////////////////////////////////////////////////////////////////////////////////////////////
// Static methods declaration
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
// Global visible variables declaration
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
//
//     B O A R D S   H A N D L I N G
//
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_open( cvt_V792_data* p_data, UINT16 base_address, long vme_handle, CVT_QTP_TYPES type)
{
	memset( p_data, 0, sizeof( cvt_V792_data));
	// Check board type
	switch( type)
	{
	case CVT_V792_TYPE_A:							/*!< \brief The board is V792A */
	case CVT_V792_TYPE_N:							/*!< \brief The board is V792N */
	case CVT_V862_TYPE_A:							/*!< \brief The board is V862A */
	case CVT_V965:									/*!< \brief The board is V965  */
	case CVT_V965_TYPE_A:							/*!< \brief The board is V965A */
	case CVT_V785_TYPE_A:							/*!< \brief The board is V785A */
	case CVT_V785_TYPE_N:							/*!< \brief The board is V785N */
	case CVT_V1785:									/*!< \brief The board is V1785 */
	case CVT_V775_TYPE_A:							/*!< \brief The board is V775A */
	case CVT_V775_TYPE_N:							/*!< \brief The board is V775N */

		break;
	default:
		TRACE1( "cvt_V792_open: bad board type: %d\n", type);
		return FALSE;
	}

	// basic data initialization
	if( !cvt_board_open( &p_data->m_common_data, base_address, vme_handle, CVT_V792_REG_TABLE))
		return FALSE;
	// board specific data initialization
	p_data->m_type= type;
	p_data->m_common_data.set_MCST_CBLT= ( BOOL (*)( void*, UINT8, MCST_CBLT_board_pos))cvt_V792_set_MCST_CBLT;

#ifdef CVT_V792_USE_DATA_QUEUE
	p_data->m_queue= malloc( CVT_V792_QUEUE_SIZE_DWORD* sizeof( UINT32));	
	if( p_data->m_queue== NULL)
	{
		TRACE( "cvt_V792_open: cannot allocate queue buffer\n");
		return FALSE;
	}
	memset( p_data->m_queue, 0, CVT_V792_QUEUE_SIZE_DWORD* sizeof( UINT32));

#endif

	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_close( cvt_V792_data* p_data)
{
#ifdef CVT_V792_USE_DATA_QUEUE
	if( p_data->m_queue!= NULL)
	{
		free( p_data->m_queue);
	}
#endif
	if( !cvt_board_close( &p_data->m_common_data))
		return FALSE;
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//
//     L E V E L   0   A P I s
//
////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////
//
//     L E V E L   1   A P I s
//
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_geo_address( cvt_V792_data* p_data, UINT8 value)
{
	UINT16 reg_value= value;
	// write register
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_GEO_ADDRESS_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_GEO_ADDRESS write failed !\n");
		return FALSE;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_get_geo_address( cvt_V792_data* p_data, UINT8 *p_value)
{
	UINT16 reg_value= 0;
	// read register
	if( !cvt_read_reg( &p_data->m_common_data, CVT_V792_GEO_ADDRESS_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_GEO_ADDRESS read failed !\n");
		return FALSE;
	}
	*p_value= (UINT8)((UINT16)( reg_value& 0x1f));
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_MCST_CBLT_address( cvt_V792_data* p_data, UINT8 value)
{
	UINT16 reg_value= value;
	// write register
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_MCST_CBLT_ADDRESS_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_MCST_CBLT_ADDRESS write failed !\n");
		return FALSE;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_get_MCST_CBLT_address( cvt_V792_data* p_data, UINT8 *p_value)
{
	UINT16 reg_value= 0;
	// read register
	if( !cvt_read_reg( &p_data->m_common_data, CVT_V792_MCST_CBLT_ADDRESS_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_GEO_ADDRESS read failed !\n");
		return FALSE;
	}
	*p_value= (UINT8)((UINT16)( reg_value& 0xff));
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_bit_set_1( cvt_V792_data* p_data, UINT16 value)
{
	UINT16 reg_value= value;
	// write register
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_SET_1_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_BIT_SET_1 write failed !\n");
		return FALSE;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_bit_clear_1( cvt_V792_data* p_data, UINT16 value)
{
	UINT16 reg_value= value;
	// write register
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_CLEAR_1_INDEX, &reg_value))
	{
		TRACE( "V792 BIT_CLEAR_1 write failed !\n");
		return FALSE;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_interrupt_level( cvt_V792_data* p_data, UINT8 value)
{
	UINT16 reg_value= value;
	// write register
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_INT_LEVEL_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_INT_LEVEL write failed !\n");
		return FALSE;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_get_interrupt_level( cvt_V792_data* p_data, UINT8 *p_value)
{
	UINT16 reg_value= 0;
	// read register
	if( !cvt_read_reg( &p_data->m_common_data, CVT_V792_INT_LEVEL_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_INT_LEVEL read failed !\n");
		return FALSE;
	}
	*p_value= (UINT8)((UINT16)( reg_value& 0x07));
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_interrupt_vector( cvt_V792_data* p_data, UINT8 value)
{
	UINT16 reg_value= value;
	// write register
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_INT_VECTOR_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_INT_VECTOR write failed !\n");
		return FALSE;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_get_interrupt_vector( cvt_V792_data* p_data, UINT8 *p_value)
{
	UINT16 reg_value= 0;
	// read register
	if( !cvt_read_reg( &p_data->m_common_data, CVT_V792_INT_VECTOR_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_INT_VECTOR read failed !\n");
		return FALSE;
	}
	*p_value= (UINT8)((UINT16)( reg_value& 0xff));
	return TRUE;
}
	
////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_get_status_1( cvt_V792_data* p_data, UINT16* p_value)
{
	UINT16 reg_value= 0;
	// read register
	if( !cvt_read_reg( &p_data->m_common_data, CVT_V792_STATUS_1_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_STATUS_1 read failed !\n");
		return FALSE;
	}
	*p_value= reg_value;
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_control_1( cvt_V792_data* p_data, UINT16 value)
{
	UINT16 reg_value= value;
	// write register
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_CONTROL_1_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_CONTROL_1 write failed !\n");
		return FALSE;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_get_control_1( cvt_V792_data* p_data, UINT16* p_value)
{
	UINT16 reg_value= 0;
	// read register
	if( !cvt_read_reg( &p_data->m_common_data, CVT_V792_CONTROL_1_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_CONTROL_1 read failed !\n");
		return FALSE;
	}
	*p_value= reg_value;
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_bitmask_control_1( cvt_V792_data* p_data, CVT_V792_CONTROL_1_MSK value)
{
	UINT16 reg_value= value;
	// set bitmask to register
	if( !cvt_set_bitmask( &p_data->m_common_data, CVT_V792_CONTROL_1_ADD, &reg_value, CVT_V792_CONTROL_1_AM, CVT_V792_CONTROL_1_DATA_SIZE))
	{
		TRACE( "V792 CVT_V792_CONTROL_1_SET failed !\n");
		return FALSE;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_clear_bitmask_control_1( cvt_V792_data* p_data, CVT_V792_CONTROL_1_MSK value)
{
	UINT16 reg_value= value;
	// set bitmask to register
	if( !cvt_clear_bitmask( &p_data->m_common_data, CVT_V792_CONTROL_1_ADD, &reg_value, CVT_V792_CONTROL_1_AM, CVT_V792_CONTROL_1_DATA_SIZE))
	{
		TRACE( "V792 CVT_V792_CONTROL_1_CLEAR failed !\n");
		return FALSE;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_ader_high( cvt_V792_data* p_data, UINT8 value)
{
	UINT16 reg_value= value;
	// write register
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_ADER_HIGH_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_ADER_HIGH write failed !\n");
		return FALSE;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_get_ader_high( cvt_V792_data* p_data, UINT8 *p_value)
{
	UINT16 reg_value= 0;
	// read register
	if( !cvt_read_reg( &p_data->m_common_data, CVT_V792_ADER_HIGH_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_ADER_HIGH read failed !\n");
		return FALSE;
	}
	*p_value= (UINT8)((UINT16)( reg_value& 0xff));
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_ader_low( cvt_V792_data* p_data, UINT8 value)
{
	UINT16 reg_value= value;
	// write register
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_ADER_LOW_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_ADER_LOW write failed !\n");
		return FALSE;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_get_ader_low( cvt_V792_data* p_data, UINT8 *p_value)
{
	UINT16 reg_value= 0;
	// read register
	if( !cvt_read_reg( &p_data->m_common_data, CVT_V792_ADER_LOW_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_ADER_LOW read failed !\n");
		return FALSE;
	}
	*p_value= (UINT8)((UINT16)( reg_value& 0xff));
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_single_shot_reset( cvt_V792_data* p_data)
{
	UINT16 reg_value= 0;
	// write register
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_SINGLE_SHOT_RESET_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_SINGLE_SHOT_RESET write failed !\n");
		return FALSE;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_MCST_CBLT_control( cvt_V792_data* p_data, UINT8 value)
{
	UINT16 reg_value= value;
	// write register
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_MCST_CBLT_CTRL_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_MCST_CBLT_CTRL write failed !\n");
		return FALSE;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_get_MCST_CBLT_control( cvt_V792_data* p_data, UINT8 *p_value)
{
	UINT16 reg_value= 0;
	// read register
	if( !cvt_read_reg( &p_data->m_common_data, CVT_V792_MCST_CBLT_CTRL_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_MCST_CBLT_CTRL read failed !\n");
		return FALSE;
	}
	*p_value= (UINT8)((UINT16)( reg_value& 0x03));
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_event_trigger( cvt_V792_data* p_data, UINT8 value)
{
	UINT16 reg_value= value;
	// write register
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_EVENT_TRG_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_EVENT_TRG write failed !\n");
		return FALSE;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_get_event_trigger( cvt_V792_data* p_data, UINT8 *p_value)
{
	UINT16 reg_value= 0;
	// read register
	if( !cvt_read_reg( &p_data->m_common_data, CVT_V792_EVENT_TRG_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_EVENT_TRG read failed !\n");
		return FALSE;
	}
	*p_value= (UINT8)((UINT16)( reg_value& 0x1f));
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_get_status_2( cvt_V792_data* p_data, UINT16* p_value)
{
	UINT16 reg_value= 0;
	// read register
	if( !cvt_read_reg( &p_data->m_common_data, CVT_V792_STATUS_2_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_STATUS_2 read failed !\n");
		return FALSE;
	}
	*p_value= reg_value;
	return TRUE;
}


////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_threshold( cvt_V792_data* p_data, cvt_V792_threshold_id threshold_id, UINT16 value)
{
	UINT16 reg_value= value;
	int max_thresholds= 0;
	int address_shift= 0;
	switch( p_data->m_type)
	{
	case CVT_V792_TYPE_A:							/*!< \brief The board is V792A */
		address_shift= 0;
		break;
	case CVT_V792_TYPE_N:							/*!< \brief The board is V792N */
		address_shift= 1;
		break;
	case CVT_V862_TYPE_A:							/*!< \brief The board is V862A */
		address_shift= 0;
		break;
	case CVT_V785_TYPE_A:							/*!< \brief The board is V785A */
		address_shift= 0;
		break;
	case CVT_V785_TYPE_N:							/*!< \brief The board is V785N */
		address_shift= 1;
		break;
	case CVT_V775_TYPE_A:							/*!< \brief The board is V775A */
		address_shift= 0;
		break;
	case CVT_V775_TYPE_N:							/*!< \brief The board is V775N */
		address_shift= 1;
		break;
	default:
		TRACE1( "cvt_V792_set_threshold: bad board type: %d\n", p_data->m_type);
		return FALSE;
	}
	max_thresholds= CVT_QTP_NUM_CHANNELS[ p_data->m_type];

	// Parameter checking
	switch( threshold_id)
	{
	case CVT_V792_THRESHOLD_ALL:						// Catch all the thresholds
		{
			// write all the thresholds
			int i;
			for( i= CVT_V792_THRESHOLD_0; i< max_thresholds; i++)
			{
				if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_THRESHOLD_0_INDEX+ ((i- CVT_V792_THRESHOLD_0)<< address_shift), &reg_value))
				{
					TRACE1( "V792 THRESHOLD %i write failed !\n", i- CVT_V792_THRESHOLD_0);
					return FALSE;
				}
			}
		}
		break;
	case CVT_V792_THRESHOLD_0:							// Threshold 0 register
	case CVT_V792_THRESHOLD_1:							// Threshold 1 register
	case CVT_V792_THRESHOLD_2:							// Threshold 2 register
	case CVT_V792_THRESHOLD_3:							// Threshold 3 register
	case CVT_V792_THRESHOLD_4:							// Threshold 4 register
	case CVT_V792_THRESHOLD_5:							// Threshold 5 register
	case CVT_V792_THRESHOLD_6:							// Threshold 6 register
	case CVT_V792_THRESHOLD_7:							// Threshold 7 register
	case CVT_V792_THRESHOLD_8:							// Threshold 8 register
	case CVT_V792_THRESHOLD_9:							// Threshold 9 register
	case CVT_V792_THRESHOLD_10:							// Threshold 10 register
	case CVT_V792_THRESHOLD_11:							// Threshold 11 register
	case CVT_V792_THRESHOLD_12:							// Threshold 12 register
	case CVT_V792_THRESHOLD_13:							// Threshold 13 register
	case CVT_V792_THRESHOLD_14:							// Threshold 14 register
	case CVT_V792_THRESHOLD_15:							// Threshold 15 register
	case CVT_V792_THRESHOLD_16:							// Threshold 16 register (V792 only)
	case CVT_V792_THRESHOLD_17:							// Threshold 17 register (V792 only)
	case CVT_V792_THRESHOLD_18:							// Threshold 18 register (V792 only)
	case CVT_V792_THRESHOLD_19:							// Threshold 19 register (V792 only)
	case CVT_V792_THRESHOLD_20:							// Threshold 20 register (V792 only)
	case CVT_V792_THRESHOLD_21:							// Threshold 21 register (V792 only)
	case CVT_V792_THRESHOLD_22:							// Threshold 22 register (V792 only)
	case CVT_V792_THRESHOLD_23:							// Threshold 23 register (V792 only)
	case CVT_V792_THRESHOLD_24:							// Threshold 24 register (V792 only)
	case CVT_V792_THRESHOLD_25:							// Threshold 25 register (V792 only)
	case CVT_V792_THRESHOLD_26:							// Threshold 26 register (V792 only)
	case CVT_V792_THRESHOLD_27:							// Threshold 27 register (V792 only)
	case CVT_V792_THRESHOLD_28:							// Threshold 28 register (V792 only)
	case CVT_V792_THRESHOLD_29:							// Threshold 29 register (V792 only)
	case CVT_V792_THRESHOLD_30:							// Threshold 30 register (V792 only)
	case CVT_V792_THRESHOLD_31:							// Threshold 31 register (V792 only)
		if( ( threshold_id- CVT_V792_THRESHOLD_0)>= max_thresholds)
			return FALSE;
		// write specific threshold
		if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_THRESHOLD_0_INDEX+ (( threshold_id- CVT_V792_THRESHOLD_0)<< address_shift), &reg_value))
		{
			TRACE1( "V792 THRESHOLD %i write failed !\n", threshold_id- CVT_V792_THRESHOLD_0);
			return FALSE;
		}
		break;
	default:
		TRACE1( "V792 Unhandled threshold id %i !\n", threshold_id);
		return FALSE;
	}
	return TRUE;
}


////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_dual_range_threshold( cvt_V792_data* p_data, cvt_V792_threshold_id threshold_id, UINT16 value_high, UINT16 value_low)
{
	int max_thresholds= 0;
	int address_shift= 0;

	switch( p_data->m_type)
	{
	case CVT_V965:									/*!< \brief The board is V965  */
		address_shift= 1;
		break;
	case CVT_V965_TYPE_A:							/*!< \brief The board is V965A */
		address_shift= 2;
		break;
	case CVT_V1785:									/*!< \brief The board is V1785 */
		address_shift= 2;
		break;
	default:
		TRACE1( "cvt_V792_set_dual_range_threshold: bad board type: %d\n", p_data->m_type);
		return FALSE;
	}
	max_thresholds= CVT_QTP_NUM_CHANNELS[ p_data->m_type];

	// Parameter checking
	switch( threshold_id)
	{
	case CVT_V792_THRESHOLD_ALL:						// Catch all the thresholds
		{
			// write all the thresholds
			int i;
			for( i= CVT_V792_THRESHOLD_0; i< max_thresholds; i++)
			{
				if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_THRESHOLD_0_INDEX+ ((i- CVT_V792_THRESHOLD_0)<< address_shift), &value_high))
				{
					TRACE1( "V792 THRESHOLD %i write failed !\n", i- CVT_V792_THRESHOLD_0);
					return FALSE;
				}
				if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_THRESHOLD_0_INDEX+ ((i- CVT_V792_THRESHOLD_0)<< address_shift)+ ( 1<< ( address_shift- 1)), &value_low))
				{
					TRACE1( "V792 THRESHOLD %i write failed !\n", i- CVT_V792_THRESHOLD_0);
					return FALSE;
				}
			}
		}
		break;
	case CVT_V792_THRESHOLD_0:							// Threshold 0 register
	case CVT_V792_THRESHOLD_1:							// Threshold 1 register
	case CVT_V792_THRESHOLD_2:							// Threshold 2 register
	case CVT_V792_THRESHOLD_3:							// Threshold 3 register
	case CVT_V792_THRESHOLD_4:							// Threshold 4 register
	case CVT_V792_THRESHOLD_5:							// Threshold 5 register
	case CVT_V792_THRESHOLD_6:							// Threshold 6 register
	case CVT_V792_THRESHOLD_7:							// Threshold 7 register
	case CVT_V792_THRESHOLD_8:							// Threshold 8 register
	case CVT_V792_THRESHOLD_9:							// Threshold 9 register
	case CVT_V792_THRESHOLD_10:							// Threshold 10 register
	case CVT_V792_THRESHOLD_11:							// Threshold 11 register
	case CVT_V792_THRESHOLD_12:							// Threshold 12 register
	case CVT_V792_THRESHOLD_13:							// Threshold 13 register
	case CVT_V792_THRESHOLD_14:							// Threshold 14 register
	case CVT_V792_THRESHOLD_15:							// Threshold 15 register
	case CVT_V792_THRESHOLD_16:							// Threshold 16 register (V792 only)
	case CVT_V792_THRESHOLD_17:							// Threshold 17 register (V792 only)
	case CVT_V792_THRESHOLD_18:							// Threshold 18 register (V792 only)
	case CVT_V792_THRESHOLD_19:							// Threshold 19 register (V792 only)
	case CVT_V792_THRESHOLD_20:							// Threshold 20 register (V792 only)
	case CVT_V792_THRESHOLD_21:							// Threshold 21 register (V792 only)
	case CVT_V792_THRESHOLD_22:							// Threshold 22 register (V792 only)
	case CVT_V792_THRESHOLD_23:							// Threshold 23 register (V792 only)
	case CVT_V792_THRESHOLD_24:							// Threshold 24 register (V792 only)
	case CVT_V792_THRESHOLD_25:							// Threshold 25 register (V792 only)
	case CVT_V792_THRESHOLD_26:							// Threshold 26 register (V792 only)
	case CVT_V792_THRESHOLD_27:							// Threshold 27 register (V792 only)
	case CVT_V792_THRESHOLD_28:							// Threshold 28 register (V792 only)
	case CVT_V792_THRESHOLD_29:							// Threshold 29 register (V792 only)
	case CVT_V792_THRESHOLD_30:							// Threshold 30 register (V792 only)
	case CVT_V792_THRESHOLD_31:							// Threshold 31 register (V792 only)
		if( ( threshold_id- CVT_V792_THRESHOLD_0)>= max_thresholds)
			return FALSE;
		// write specific threshold
		if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_THRESHOLD_0_INDEX+ (( threshold_id- CVT_V792_THRESHOLD_0)<< address_shift), &value_high))
		{
			TRACE1( "V792 THRESHOLD %i write failed !\n", threshold_id- CVT_V792_THRESHOLD_0);
			return FALSE;
		}
		if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_THRESHOLD_0_INDEX+ (( threshold_id- CVT_V792_THRESHOLD_0)<< address_shift)+ ( 1<< ( address_shift- 1)), &value_low))
		{
			TRACE1( "V792 THRESHOLD %i write failed !\n", threshold_id- CVT_V792_THRESHOLD_0);
			return FALSE;
		}
		break;
	default:
		TRACE1( "V792 Unhandled threshold id %i !\n", threshold_id);
		return FALSE;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_bit_set_2( cvt_V792_data* p_data, UINT16 value)
{
	UINT16 reg_value= value;
	// write register
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_SET_2_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_BIT_SET_2 write failed !\n");
		return FALSE;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_bit_clear_2( cvt_V792_data* p_data, UINT16 value)
{
	UINT16 reg_value= value;
	// write register
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_CLEAR_2_INDEX, &reg_value))
	{
		TRACE( "V792 BIT_CLEAR_2 write failed !\n");
		return FALSE;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//
//     L E V E L   2   A P I s
//
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_sliding_scale( cvt_V792_data* p_data, BOOL enable)
{
	UINT16 reg_value;
	reg_value= CVT_V792_BSC2_SLIDE_EN_MSK;
	// set the sliding scale enable bit
	if( !cvt_write_reg( &p_data->m_common_data, enable? CVT_V792_BIT_SET_2_INDEX: CVT_V792_BIT_CLEAR_2_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_BSC2_SLIDE_EN_MSK bit set failed !\n");
		return FALSE;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_zero_suppression( cvt_V792_data* p_data, BOOL enable, BOOL step_threshold, const UINT16* thresholds_buff)
{
	UINT16 reg_value;
	int max_thresholds= 0;
	int address_shift= 0;
	switch( p_data->m_type)
	{
	case CVT_V792_TYPE_A:							/*!< \brief The board is V792A */
		address_shift= 0;
		break;
	case CVT_V792_TYPE_N:							/*!< \brief The board is V792N */
		address_shift= 1;
		break;
	case CVT_V862_TYPE_A:							/*!< \brief The board is V862A */
		address_shift= 0;
		break;
	case CVT_V785_TYPE_A:							/*!< \brief The board is V785A */
		address_shift= 0;
		break;
	case CVT_V785_TYPE_N:							/*!< \brief The board is V785N */
		address_shift= 1;
		break;
	case CVT_V775_TYPE_A:							/*!< \brief The board is V775A */
		address_shift= 0;
		break;
	case CVT_V775_TYPE_N:							/*!< \brief The board is V775N */
		address_shift= 1;
		break;
	default:
		TRACE1( "cvt_V792_set_zero_suppression: bad board type: %d\n", p_data->m_type);
		return FALSE;
	}
	max_thresholds= CVT_QTP_NUM_CHANNELS[ p_data->m_type];

	if( !enable)
	{
		// Set the low threshold disable bit
		reg_value= CVT_V792_BSC2_LOW_THR_DIS_MSK;
		if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_SET_2_INDEX, &reg_value))
		{
			TRACE( "V792 CVT_V792_BSC2_LOW_THR_DIS_MSK bit set failed !\n");
			return FALSE;
		}
	}
	else
	{
		// setup thresholds
		reg_value= CVT_V792_BSC2_LOW_THR_DIS_MSK;
		if( thresholds_buff!= NULL)
		{
			int i;
			// write all the thresholds
			for( i= CVT_V792_THRESHOLD_0; i< max_thresholds; i++)
			{
				if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_THRESHOLD_0_INDEX+ ((i- CVT_V792_THRESHOLD_0)<< address_shift), thresholds_buff++))
				{
					TRACE1( "V792 THRESHOLD %i write failed !\n", i- CVT_V792_THRESHOLD_0);
					return FALSE;
				}
			}
		}

		// Clear the step threshold bit
		reg_value= CVT_V792_BSC2_STEP_TH_MSK;
		if( step_threshold)
		{
			if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_SET_2_INDEX, &reg_value))
			{
				TRACE( "V792 CVT_V792_BSC2_STEP_TH_MSK bit set failed !\n");
				return FALSE;
			}
		}
		else
		{
			if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_CLEAR_2_INDEX, &reg_value))
			{
				TRACE( "V792 CVT_V792_BSC2_STEP_TH_MSK bit clear failed !\n");
				return FALSE;
			}
		}

		// clear the low threshold disable bit
		reg_value= CVT_V792_BSC2_LOW_THR_DIS_MSK;
		if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_CLEAR_2_INDEX, &reg_value))
		{
			TRACE( "V792 CVT_V792_BSC2_LOW_THR_DIS_MSK bit clear failed !\n");
			return FALSE;
		}
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_dual_range_zero_suppression( cvt_V792_data* p_data, BOOL enable, BOOL step_threshold, const UINT16* thresholds_buff_high, const UINT16* thresholds_buff_low)
{
	UINT16 reg_value;
	int max_thresholds= 0;
	int address_shift= 0;

	switch( p_data->m_type)
	{
	case CVT_V965:									/*!< \brief The board is V965  */
		address_shift= 1;
		break;
	case CVT_V965_TYPE_A:							/*!< \brief The board is V965A */
		address_shift= 2;
		break;
	case CVT_V1785:									/*!< \brief The board is V1785 */
		address_shift= 2;
		break;
	default:
		TRACE1( "cvt_V792_set_dual_range_threshold: bad board type: %d\n", p_data->m_type);
		return FALSE;
	}
	max_thresholds= CVT_QTP_NUM_CHANNELS[ p_data->m_type];
	if( !enable)
	{
		// Set the low threshold disable bit
		reg_value= CVT_V792_BSC2_LOW_THR_DIS_MSK;
		if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_SET_2_INDEX, &reg_value))
		{
			TRACE( "V792 CVT_V792_BSC2_LOW_THR_DIS_MSK bit set failed !\n");
			return FALSE;
		}
	}
	else
	{
		// setup thresholds
		reg_value= CVT_V792_BSC2_LOW_THR_DIS_MSK;
		if( thresholds_buff_high!= NULL)
		{
			int i;
			// write all the thresholds
			for( i= CVT_V792_THRESHOLD_0; i< max_thresholds; i++)
			{
				if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_THRESHOLD_0_INDEX+ (( i- CVT_V792_THRESHOLD_0)<< address_shift), thresholds_buff_high++))
				{
					TRACE1( "V792 THRESHOLD %i write failed !\n", threshold_id- CVT_V792_THRESHOLD_0);
					return FALSE;
				}
			}
		}
		if( thresholds_buff_low!= NULL)
		{
			int i;
			// write all the thresholds
			for( i= CVT_V792_THRESHOLD_0; i< max_thresholds; i++)
			{
				if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_THRESHOLD_0_INDEX+ (( i- CVT_V792_THRESHOLD_0)<< address_shift)+ ( 1<< ( address_shift- 1)), thresholds_buff_low++))
				{
					TRACE1( "V792 THRESHOLD %i write failed !\n", threshold_id- CVT_V792_THRESHOLD_0);
					return FALSE;
				}
			}
		}

		// Clear the step threshold bit
		reg_value= CVT_V792_BSC2_STEP_TH_MSK;
		if( step_threshold)
		{
			if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_SET_2_INDEX, &reg_value))
			{
				TRACE( "V792 CVT_V792_BSC2_STEP_TH_MSK bit set failed !\n");
				return FALSE;
			}
		}
		else
		{
			if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_CLEAR_2_INDEX, &reg_value))
			{
				TRACE( "V792 CVT_V792_BSC2_STEP_TH_MSK bit clear failed !\n");
				return FALSE;
			}
		}

		// clear the low threshold disable bit
		reg_value= CVT_V792_BSC2_LOW_THR_DIS_MSK;
		if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_CLEAR_2_INDEX, &reg_value))
		{
			TRACE( "V792 CVT_V792_BSC2_LOW_THR_DIS_MSK bit clear failed !\n");
			return FALSE;
		}
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_overflow_suppression( cvt_V792_data* p_data, BOOL enable)
{
	UINT16 reg_value;

	// Overflow suppression disable bit
	reg_value= CVT_V792_BSC2_OVER_RANGE_DIS_MSK;
	if( !enable)
	{
		// Set the Overflow suppression disable bit
		if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_SET_2_INDEX, &reg_value))
		{
			TRACE( "V792 CVT_V792_BSC2_OVER_RANGE_DIS_MSK bit set failed !\n");
			return FALSE;
		}
	}
	else
	{
		// Clear the Overflow suppression disable bit
		if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_CLEAR_2_INDEX, &reg_value))
		{
			TRACE( "V792 CVT_V792_BSC2_OVER_RANGE_DIS_MSK bit clear failed !\n");
			return FALSE;
		}
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_MEB( cvt_V792_data* p_data, BOOL empty_enable, BOOL auto_incr_enable)
{
	UINT16 reg_value;

	// Empty enable bit
	reg_value= CVT_V792_BSC2_EMPTY_EN_MSK;
	if( !empty_enable)
	{
		// Clear the empty enable bit
		if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_CLEAR_2_INDEX, &reg_value))
		{
			TRACE( "V792 CVT_V792_BSC2_EMPTY_EN_MSK bit clear failed !\n");
			return FALSE;
		}
	}
	else
	{
		// Set the empty enable bit
		if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_SET_2_INDEX, &reg_value))
		{
			TRACE( "V792 CVT_V792_BSC2_EMPTY_EN_MSK bit set failed !\n");
			return FALSE;
		}
	}
	// Auto increment enable bit
	reg_value= CVT_V792_BSC2_AUTO_INC_MSK;
	if( !auto_incr_enable)
	{
		// Clear the auto increment enable bit
		if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_CLEAR_2_INDEX, &reg_value))
		{
			TRACE( "V792 CVT_V792_BSC2_AUTO_INC_MSK bit clear failed !\n");
			return FALSE;
		}
	}
	else
	{
		// Set the empty enable bit
		if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_SET_2_INDEX, &reg_value))
		{
			TRACE( "V792 CVT_V792_BSC2_AUTO_INC_MSK bit set failed !\n");
			return FALSE;
		}
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_event_counter( cvt_V792_data* p_data, BOOL count_all_events, BOOL reset_count)
{
	UINT16 reg_value;

	// Count All events enable bit
	reg_value= CVT_V792_BSC2_ALL_TRG_MSK;
	if( !count_all_events)
	{
		// Clear the Count All events enable bit
		if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_CLEAR_2_INDEX, &reg_value))
		{
			TRACE( "V792 CVT_V792_BSC2_ALL_TRG_MSK bit clear failed !\n");
			return FALSE;
		}
	}
	else
	{
		// Set the Count All events enable bit
		if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_SET_2_INDEX, &reg_value))
		{
			TRACE( "V792 CVT_V792_BSC2_ALL_TRG_MSK bit set failed !\n");
			return FALSE;
		}
	}
	
	// Counter reset
	reg_value= 0;
	if( reset_count)
	{
		// Resets the counter
		if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_EVENT_COUNT_RESET_INDEX, &reg_value))
		{
			TRACE( "V792 CVT_V792_EVENT_COUNT_RESET write failed !\n");
			return FALSE;
		}
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_read_MEB( cvt_V792_data* p_data, void* p_buff, UINT32* p_buff_size)
{
	BOOL is_berr= FALSE;
	if( !cvt_FIFO_BLT_read( &p_data->m_common_data, CVT_V792_OUT_BUFFER_ADD, p_buff, *p_buff_size, p_buff_size, CVT_V792_OUT_BUFFER_AM, CVT_V792_OUT_BUFFER_DATA_SIZE, &is_berr))
	{
		return FALSE;
	}
	// V792 return an invalid datum if output buffer empty, so discard it
	if( *p_buff_size== sizeof( UINT32))
	{
		if( CVT_QTP_IS_NOT_VALID_DATUM( *((UINT32*)p_buff)))
			*p_buff_size= 0;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_fast_clear_window( cvt_V792_data* p_data, UINT16 value)
{
	UINT16 reg_value= value;
	// write register
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_FCLR_WND_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_FCLR_WND write failed !\n");
		return FALSE;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_acquisition_mode( cvt_V792_data* p_data, BOOL sliding_scale_enable, BOOL zero_suppression_enable, BOOL overflow_suppression_enable, BOOL valid_suppression_enable, BOOL common_stop_enable, BOOL empty_enable, BOOL count_all_events)
{
	UINT16 set_msk= 0;
	UINT16 clear_msk= 0;

	//
	// Sliding scale
	set_msk|= ( sliding_scale_enable)? CVT_V792_BSC2_SLIDE_EN_MSK: 0;
	clear_msk|= ( sliding_scale_enable)? 0: CVT_V792_BSC2_SLIDE_EN_MSK;
	//
	// zero suppression enable bit
	set_msk|= ( zero_suppression_enable)? 0: CVT_V792_BSC2_LOW_THR_DIS_MSK;
	clear_msk|= ( zero_suppression_enable)? CVT_V792_BSC2_LOW_THR_DIS_MSK: 0;
	//
	// Overflow suppression disable bit
	set_msk|= ( overflow_suppression_enable)? 0: CVT_V792_BSC2_OVER_RANGE_DIS_MSK;
	clear_msk|= ( overflow_suppression_enable)? CVT_V792_BSC2_OVER_RANGE_DIS_MSK: 0;


	//
	// valid suppression disable bit
	set_msk|= ( valid_suppression_enable)? 0: CVT_V792_BSC2_VALID_DIS_MSK;
	clear_msk|= ( valid_suppression_enable)? CVT_V792_BSC2_VALID_DIS_MSK: 0;

	//
	// common stop disable bit
	clear_msk|= ( common_stop_enable)? 0: CVT_V792_BSC2_COMMON_STOP_MSK;
	set_msk|= ( common_stop_enable)? CVT_V792_BSC2_COMMON_STOP_MSK: 0;
	//
	// Empty enable bit
	set_msk|= ( empty_enable)? CVT_V792_BSC2_EMPTY_EN_MSK: 0;
	clear_msk|= ( empty_enable)? 0: CVT_V792_BSC2_EMPTY_EN_MSK;
	//
	// All trigger bit (count all events)
	set_msk|= ( count_all_events)? CVT_V792_BSC2_ALL_TRG_MSK: 0;
	clear_msk|= ( count_all_events)? 0: CVT_V792_BSC2_ALL_TRG_MSK;
	
	//
	// bitmask set
	if( set_msk)
	{
		if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_SET_2_INDEX, &set_msk))
		{
			TRACE( "V792 CVT_V792_BIT_SET_2 write failed !\n");
			return FALSE;
		}
	}
	//
	// bitmask clear
	if( clear_msk)
	{
		// clear the enable bit
		if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_CLEAR_2_INDEX, &clear_msk))
		{
			TRACE( "V792 CVT_V792_BIT_CLEAR_2 write failed !\n");
			return FALSE;
		}
	}
	return TRUE;
}
////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_interrupt( cvt_V792_data* p_data, UINT8 level, UINT8 vector, UINT8 event_number)
{
	UINT16 reg_value;
	//
	// interrupt level
	reg_value= level& 0x07;
	if( level!= reg_value)
	{
		TRACE( "V792 bad interrupt level !\n");
		return FALSE;
	}
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_INT_LEVEL_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_INT_LEVEL write failed !\n");
		return FALSE;
	}
	//
	// interrupt vector
	reg_value= vector;
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_INT_VECTOR_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_INT_VECTOR write failed !\n");
		return FALSE;
	}
	//
	// event number (event trigger register)
	reg_value= event_number& 0x1f;
	if( event_number!= reg_value)
	{
		TRACE( "V792 bad interrupt event number !\n");
		return FALSE;
	}
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_EVENT_TRG_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_EVENT_TRG write failed !\n");
		return FALSE;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_readout_mode( cvt_V792_data* p_data, BOOL bus_error_enable, BOOL block_end_enable, BOOL align64_enable)
{
	UINT16 set_msk= 0;
	UINT16 clear_msk= 0;
	
	//
	// Bus error
	set_msk|= ( bus_error_enable)? CVT_V792_CTRL1_BERR_ENABLE_MSK: 0;
	clear_msk|= ( bus_error_enable)? 0: CVT_V792_CTRL1_BERR_ENABLE_MSK;
	//
	// Block end
	set_msk|= ( block_end_enable)? CVT_V792_CTRL1_BLKEND_MSK: 0;
	clear_msk|= ( block_end_enable)? 0: CVT_V792_CTRL1_BLKEND_MSK;
	//
	// Align 64
	set_msk|= ( align64_enable)? CVT_V792_CTRL1_ALIGN64_MSK: 0;
	clear_msk|= ( align64_enable)? 0: CVT_V792_CTRL1_ALIGN64_MSK;
	
	//
	// bitmask set
	if( set_msk)
	{
		if( !cvt_set_bitmask_reg( &p_data->m_common_data, CVT_V792_CONTROL_1_INDEX, &set_msk))
		{
			TRACE( "V792 CVT_V792_CONTROL_1 write failed !\n");
			return FALSE;
		}
	}
	//
	// bitmask clear
	if( clear_msk)
	{
		// clear the enable bit
		if( !cvt_clear_bitmask_reg( &p_data->m_common_data, CVT_V792_CONTROL_1_INDEX, &clear_msk))
		{
			TRACE( "V792 CVT_V792_CONTROL_1 write failed !\n");
			return FALSE;
		}
	}
	return TRUE;
}
////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_get_status( cvt_V792_data* p_data, BOOL *p_is_data_ready, BOOL *p_is_busy, BOOL *p_is_term_on, BOOL *p_is_buffer_full)
{
	UINT16 reg_value= 0;
	//
	// read status register 1
	if( !cvt_read_reg( &p_data->m_common_data, CVT_V792_STATUS_1_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_STATUS_1 read failed !\n");
		return FALSE;
	}
	*p_is_data_ready= ( reg_value& CVT_V792_STS1_DREADY_MSK)? TRUE: FALSE;
	*p_is_busy= ( reg_value& CVT_V792_STS1_BUSY_MSK)? TRUE: FALSE;
	*p_is_term_on= ( reg_value& CVT_V792_STS1_TERM_ON_MSK)? TRUE: FALSE;
	//
	// read status register 2
	if( !cvt_read_reg( &p_data->m_common_data, CVT_V792_STATUS_2_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_STATUS_2 read failed !\n");
		return FALSE;
	}
	*p_is_buffer_full= ( reg_value& CVT_V792_STS2_BUFFER_FULL_MSK)? TRUE: FALSE;

	return TRUE;
}


////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_crate_number( cvt_V792_data* p_data, UINT8 crate_number)
{
	UINT16 reg_value;
	//
	// Crate number
	reg_value= crate_number;
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_CRATE_SELECT_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_CRATE_SELECT write failed !\n");
		return FALSE;
	}

	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_pedestal( cvt_V792_data* p_data, UINT8 value)
{
	UINT16 reg_value= value;
	//
	// Iped register
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_IPED_INDEX, &reg_value))
	{
		TRACE( "V792 IPED write failed !\n");
		return FALSE;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_get_event_counter( cvt_V792_data* p_data, UINT32* p_counter)
{
	UINT16 reg_value;
	*p_counter= 0;
	// Counter L
	if( !cvt_read_reg( &p_data->m_common_data, CVT_V792_EVENT_COUNTER_LOW_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_EVENT_COUNTER_LOW_INDEX read failed !\n");
		return FALSE;
	}
	*p_counter= reg_value;
	// Counter H
	if( !cvt_read_reg( &p_data->m_common_data, CVT_V792_EVENT_COUNTER_HIGH_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_EVENT_COUNTER_HIGH_INDEX read failed !\n");
		return FALSE;
	}
	*p_counter|= ( ((UINT32)( reg_value& 0x00ff))<< 16);

	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_get_system_info( cvt_V792_data* p_data, UINT16 *p_firmware_rev, UINT8 *p_piggy_back_type, UINT16 *p_serial_number)
{
	UINT16 reg_value= 0;
	//
	// Firmware revision register
	if( !cvt_read_reg( &p_data->m_common_data, CVT_V792_FW_REV_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_FW_REV read failed !\n");
		return FALSE;
	}
	*p_firmware_rev= reg_value;
	//
	// Piggy back type
	if( !cvt_read_reg( &p_data->m_common_data, CVT_V792_STATUS_2_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_STATUS_2 read failed !\n");
		return FALSE;
	}
	*p_piggy_back_type= CVT_V792_GET_PIGGY_BACK_TYPE( reg_value);
	//
	// Serial number
	if( !cvt_read_reg( &p_data->m_common_data, CVT_V792_ROM_SERIAL_0_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_ROM_SERIAL_LSB read failed !\n");
		return FALSE;
	}
	*p_serial_number= reg_value& 0x00ff;
	if( !cvt_read_reg( &p_data->m_common_data, CVT_V792_ROM_SERIAL_1_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_ROM_SERIAL_MSB read failed !\n");
		return FALSE;
	}
	*p_serial_number|= (reg_value& 0x00ff)<< 8;

	return TRUE;

}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
const char * cvt_V792_get_piggy_back_str( cvt_V792_data* p_data, UINT8 piggy_back_type)
{
	const char *pback_str[]= {	"[CAEN TAC 32ch - A786]   ","[CAEN PK 32ch - A787]    ","[CAEN QDC 32ch - A796]   ",
								"[INFN TDC 32ch]          ","[INFN PK  32ch]          ","[UNKOWN TYPE]            ",
								"[CAEN QDC-IG 32ch - A863]","[UNKOWN TYPE]            ","[CAEN TAC 16ch - A986]   ",
								"[CAEN PK 16ch - A987]    ","[CAEN QDC 16ch - A855]   ","[UNKOWN TYPE]            ",
								"[UNKOWN TYPE]            ","[UNKOWN TYPE]            ","[CAEN DUAL-RANGE  8ch]   ",
								"[UNKOWN TYPE]            "};
	if( piggy_back_type>= sizeof( pback_str))
	{
		return "[UNKOWN TYPE]            ";
	}
	return pback_str[ piggy_back_type];

}
////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_software_reset( cvt_V792_data* p_data)
{
/* Pulses SOFT_RESET_BIT
	UINT16 reg_value= (UINT16)CVT_V792_BSC1_SOFT_RESET_MSK;
	// write bit set register
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_SET_1_INDEX, &reg_value))
	{
		TRACE( "V792 SW_RESET write failed !\n");
		return FALSE;
	}
	// write bit clear register
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_CLEAR_1_INDEX, &reg_value))
	{
		TRACE( "V792 SW_RESET write failed !\n");
		return FALSE;
	}
*/
	//
	// Write SINGLE_SHOT_RESET register
	UINT16 reg_value= 0;
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_SINGLE_SHOT_RESET_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_SINGLE_SHOT_RESET write failed !\n");
		return FALSE;
	}

	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_data_clear( cvt_V792_data* p_data)
{
	UINT16 reg_value= (UINT16)CVT_V792_BSC2_CLEAR_DATA_MSK;
	// write bit set register
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_SET_2_INDEX, &reg_value))
	{
		TRACE( "V792 SW_CLEAR write failed !\n");
		return FALSE;
	}
	// write bit clear register
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_BIT_CLEAR_2_INDEX, &reg_value))
	{
		TRACE( "V792 SW_CLEAR write failed !\n");
		return FALSE;
	}

	//
	// Event Counter reset
	reg_value= 0;
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_EVENT_COUNT_RESET_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_EVENT_COUNT_RESET write failed !\n");
		return FALSE;
	}
	return TRUE;
}


#ifdef CVT_V792_USE_DATA_QUEUE
////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_read_data( cvt_V792_data *p_data) 
{
	BOOL ret_val= TRUE;
	UINT32 size_bytes;
	UINT32 *buffer= NULL;
	
	size_bytes= cvt_V792_get_queue_free( p_data)* sizeof( UINT32);

	if( !size_bytes)
	{
		// Queue full
		return TRUE;
	}
	buffer= malloc( size_bytes);
	if( buffer== NULL)
	{
		TRACE( "cvt_V792_read_data: cannot allocate temporary buffer\n");
		return FALSE;
	}
	// Read MEB
	if( !cvt_V792_read_MEB( p_data, buffer, &size_bytes))
	{
		ret_val= FALSE;
		goto exit_point;
	}
	// Put into queue
	if( !cvt_V792_inqueue( p_data, buffer, size_bytes/ sizeof( UINT32))) 
	{
		ret_val= FALSE;
		goto exit_point;
	}

exit_point:
	if( buffer) 
	{
		free( buffer);
	}
	return ret_val;

}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_peek_event( cvt_V792_data *p_data, UINT32 *out_buff, long *p_out_buff_size, UINT32 *p_event_count, BOOL dequeue)
{
	long in_org_pos= p_data->m_queue_ini;
	int step= 0;
	UINT32 data;
	long queue_ptr= in_org_pos;
	long out_buff_size= *p_out_buff_size;
	long num_word_to_read= 0;

	// scan input stream
	while( queue_ptr!= p_data->m_queue_end)
	{
		data= p_data->m_queue[ queue_ptr];
		if( ++queue_ptr>= CVT_V792_QUEUE_SIZE_DWORD)
		{
			queue_ptr= 0;
		}
		switch( step)
		{
			case 0:
				//
				// wait for header
				if( !CVT_QTP_IS_HEADER( data))
				{
					// misplaced data ?!?!?
					// skip it
					cvt_V792_dequeue( p_data, NULL, 1);
					in_org_pos= queue_ptr;
					break;
				}
				// Header found:
				num_word_to_read= 1;
				++step;
				break;
			case 1:
				//
				// wait for EOB
				++num_word_to_read;
				// check if EOB
				if( !CVT_QTP_IS_EOB( data))
				{
					// channels data
					break;
				}
				// EOB found
				// check event count
				*p_event_count= CVT_QTP_GET_EOB_EVENT_COUNT( data);
				if( out_buff_size< num_word_to_read)
				{
					// p_data->m_queue_ini= in_org_pos;
					*p_out_buff_size= 0;
					return FALSE;
				}
				// event count matching
				if( queue_ptr>= in_org_pos)
				{
					memcpy( out_buff, p_data->m_queue+ in_org_pos, num_word_to_read*sizeof( UINT32));
					// p_data->m_queue_ini= queue_ptr;
				}
				else
				{
					// HACK: debug < or <= condition
					long first_block_word= CVT_V792_QUEUE_SIZE_DWORD- in_org_pos;
					long second_block_word= num_word_to_read- first_block_word;
					memcpy( out_buff, p_data->m_queue+ in_org_pos, first_block_word* sizeof(UINT32));
					memcpy( out_buff+ first_block_word, p_data->m_queue, second_block_word* sizeof(UINT32));
					// p_data->m_queue_ini= second_block_word;
				}
				if( dequeue) 
				{
					cvt_V792_dequeue( p_data, NULL, num_word_to_read);
				}

				*p_out_buff_size= num_word_to_read;
				// done
				return TRUE;
		}
	}
	// p_data->m_queue_ini= in_org_pos;
	*p_out_buff_size= 0;
	return FALSE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_peek_channel_data( cvt_V792_data *p_data, UINT8 channel_index, UINT16 *out_buff, long *p_out_buff_size, UINT32 *p_event_count, BOOL dequeue)
{
	long in_org_pos= p_data->m_queue_ini;
	int step= 0;
	UINT32 data;
	long queue_ptr= in_org_pos;
	long queue_data_ptr= in_org_pos;
	long out_buff_size= *p_out_buff_size;
	long num_word_to_read= 0;

	*p_out_buff_size= 0;
	// scan input stream
	while( queue_ptr!= p_data->m_queue_end)
	{
		data= p_data->m_queue[ queue_ptr];
		if( ++queue_ptr>= CVT_V792_QUEUE_SIZE_DWORD)
		{
			queue_ptr= 0;
		}
		switch( step)
		{
			case 0:
				//
				// wait for header
				if( !CVT_QTP_IS_HEADER( data))
				{
					// misplaced data ?!?!?
					// skip it
					cvt_V792_dequeue( p_data, NULL, 1);
					in_org_pos= queue_ptr;
					break;
				}
				// Header found:
				num_word_to_read= 1;
				queue_data_ptr= queue_ptr;
				++step;
				break;
			case 1:
				//
				// wait for EOB
				// check if EOB
				if( !CVT_QTP_IS_EOB( data))
				{
					++num_word_to_read;
					// channels data
					break;
				}
				// EOB found
				// check event count
				*p_event_count= CVT_QTP_GET_EOB_EVENT_COUNT( data);
				{
					int count_word; 
					*p_out_buff_size= 0;
					for( count_word= 0; count_word< num_word_to_read; count_word++)
					{
						data= p_data->m_queue[ queue_data_ptr];
						if( ++queue_data_ptr>= CVT_V792_QUEUE_SIZE_DWORD)
						{
							queue_data_ptr= 0;
						}
						
						if( ( data& CVT_QTP_DATA_TYPE_MSK)== CVT_QTP_DATUM)
						{
							// DATUM
							int ch= -1;
							UINT32 geo;
							UINT32 adc;
							UINT32 un;
							UINT32 ov;

							switch ( p_data->m_type) 
							{
							case CVT_V792_TYPE_A:
								ch= CVT_V792A_GET_DATUM_CH( data);
								geo= CVT_V792_GET_DATUM_GEO( data);
								adc= CVT_V792_GET_DATUM_ADC( data);
								un= CVT_V792_GET_DATUM_UN( data);
								ov= CVT_V792_GET_DATUM_OV( data);
								break;
							case CVT_V792_TYPE_N:
								ch= CVT_V792N_GET_DATUM_CH( data);
								geo= CVT_V792_GET_DATUM_GEO( data);
								adc= CVT_V792_GET_DATUM_ADC( data);
								un= CVT_V792_GET_DATUM_UN( data);
								ov= CVT_V792_GET_DATUM_OV( data);
								break;
							case CVT_V965:
								ch= CVT_V965_GET_DATUM_CH( data);
								geo= CVT_V965_GET_DATUM_GEO( data);
								adc= CVT_V965_GET_DATUM_ADC( data);
								un= CVT_V965_GET_DATUM_UN( data);
								ov= CVT_V965_GET_DATUM_OV( data);
								break;
							case CVT_V965_TYPE_A:
								ch= CVT_V965A_GET_DATUM_CH( data);
								geo= CVT_V965_GET_DATUM_GEO( data);
								adc= CVT_V965_GET_DATUM_ADC( data);
								un= CVT_V965_GET_DATUM_UN( data);
								ov= CVT_V965_GET_DATUM_OV( data);
								break;
							case CVT_V785_TYPE_A:
								ch= CVT_V785A_GET_DATUM_CH( data);
								geo= CVT_V785_GET_DATUM_GEO( data);
								adc= CVT_V785_GET_DATUM_ADC( data);
								un= CVT_V785_GET_DATUM_UN( data);
								ov= CVT_V785_GET_DATUM_OV( data);
								break;
							case CVT_V785_TYPE_N:
								ch= CVT_V785N_GET_DATUM_CH( data);
								geo= CVT_V785_GET_DATUM_GEO( data);
								adc= CVT_V785_GET_DATUM_ADC( data);
								un= CVT_V785_GET_DATUM_UN( data);
								ov= CVT_V785_GET_DATUM_OV( data);
								break;
							case CVT_V1785:
								ch= CVT_V1785_GET_DATUM_CH( data);
								geo= CVT_V1785_GET_DATUM_GEO( data);
								adc= CVT_V1785_GET_DATUM_ADC( data);
								un= CVT_V1785_GET_DATUM_UN( data);
								ov= CVT_V1785_GET_DATUM_OV( data);
								break;
							case CVT_V775_TYPE_A:
								ch= CVT_V775A_GET_DATUM_CH( data);
								geo= CVT_V775_GET_DATUM_GEO( data);
								adc= CVT_V775_GET_DATUM_ADC( data);
								un= CVT_V775_GET_DATUM_UN( data);
								ov= CVT_V775_GET_DATUM_OV( data);
								break;
							case CVT_V775_TYPE_N:
								ch= CVT_V775N_GET_DATUM_CH( data);
								geo= CVT_V775_GET_DATUM_GEO( data);
								adc= CVT_V775_GET_DATUM_ADC( data);
								un= CVT_V775_GET_DATUM_UN( data);
								ov= CVT_V775_GET_DATUM_OV( data);
								break;
							default:
								// Unmatched case !??!
								TRACE1( "cvt_V792_peek_channel_data: bad board type: %d\n", p_data->m_type);
								*p_out_buff_size= 0;
								return FALSE;
							}
							if( ch!= (int)channel_index)
							{
								continue;
							}
							// Store the adc value
							if( out_buff_size< ( *p_out_buff_size)+ 1)
							{
								// p_data->m_queue_ini= in_org_pos;
								*p_out_buff_size= 0;
								return FALSE;
							}
							out_buff[ (*p_out_buff_size)++]= (UINT16)adc;
						}
					}
				}
				if( dequeue) 
				{
					cvt_V792_dequeue( p_data, NULL, num_word_to_read);
				}
				// done
				return TRUE;
		}
	}
	// p_data->m_queue_ini= in_org_pos;
	*p_out_buff_size= 0;
	return FALSE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_inqueue( cvt_V792_data* p_data, const UINT32* in_buff, UINT32 in_buff_size)
{
	if( !in_buff_size)
		return TRUE;
	if( cvt_V792_get_queue_free( p_data)< (long)in_buff_size)
	{
		// No enough free space
		return FALSE;
	}
	if( p_data->m_queue_end+ in_buff_size< CVT_V792_QUEUE_SIZE_DWORD)
	{
		memcpy( p_data->m_queue+ p_data->m_queue_end, in_buff, in_buff_size*sizeof( UINT32));
		p_data->m_queue_end+= in_buff_size;
	}
	else
	{
		long first_block_word= CVT_V792_QUEUE_SIZE_DWORD- p_data->m_queue_end;
		long second_block_word= in_buff_size- first_block_word;
		memcpy( p_data->m_queue+ p_data->m_queue_end, in_buff, first_block_word* sizeof(UINT32));
		memcpy( p_data->m_queue, in_buff+ first_block_word, second_block_word* sizeof(UINT32));
		p_data->m_queue_end= second_block_word;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_dequeue( cvt_V792_data* p_data, UINT32 *out_buff, UINT32 out_buff_size)
{
	if( !out_buff_size)
		return TRUE;
	if( cvt_V792_get_queue_length( p_data)< (long)out_buff_size)
	{
		// No enough data
		return FALSE;
	}
	if( p_data->m_queue_ini+ out_buff_size< CVT_V792_QUEUE_SIZE_DWORD)
	{
		if( out_buff)
		{
			memcpy( out_buff, p_data->m_queue+ p_data->m_queue_ini, out_buff_size*sizeof( UINT32));
		}
		p_data->m_queue_ini+= out_buff_size;
	}
	else
	{
		long first_block_word= CVT_V792_QUEUE_SIZE_DWORD- p_data->m_queue_ini;
		long second_block_word= out_buff_size- first_block_word;
		if( out_buff)
		{
			memcpy( out_buff, p_data->m_queue+ p_data->m_queue_ini, first_block_word* sizeof(UINT32));
			memcpy( out_buff+ first_block_word, p_data->m_queue, second_block_word* sizeof(UINT32));
		}
		p_data->m_queue_ini= second_block_word;
	}
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
long cvt_V792_get_queue_free( cvt_V792_data* p_data)
{
	return (long)CVT_V792_QUEUE_SIZE_DWORD- (( p_data->m_queue_end- p_data->m_queue_ini+ CVT_V792_QUEUE_SIZE_DWORD)%CVT_V792_QUEUE_SIZE_DWORD);
}

////////////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////////////
long cvt_V792_get_queue_length( cvt_V792_data* p_data)
{
	return (long)(( p_data->m_queue_end- p_data->m_queue_ini+ CVT_V792_QUEUE_SIZE_DWORD)%CVT_V792_QUEUE_SIZE_DWORD);
}

#endif // CVT_V792_USE_DATA_QUEUE

////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////
BOOL cvt_V792_set_MCST_CBLT( cvt_V792_data* p_data, UINT8 address, MCST_CBLT_board_pos pos)
{
	UINT16 reg_value;
	//
	// MCST/CBLT address register
	reg_value= address;
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_MCST_CBLT_ADDRESS_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_MCST_CBLT_ADDRESS write failed !\n");
		return FALSE;
	}
	//
	// MCST/CBLT control register
	switch( pos)
	{
	case MCST_CBLT_board_pos_first:
		reg_value= CVT_V792_MCCTRL_FIRST_BOARD_MSK;
		break;
	case MCST_CBLT_board_pos_mid:
		reg_value= CVT_V792_MCCTRL_MID_BOARD_MSK;
		break;
	case MCST_CBLT_board_pos_last:
		reg_value= CVT_V792_MCCTRL_LAST_BOARD_MSK;
		break;
	case MCST_CBLT_board_disabled:
		reg_value= CVT_V792_MCCTRL_DISABLED_BOARD_MSK;
		break;
	default:
		TRACE1( "V792 cvt_V792_set_MCST_CBLT bad position indentifier '%d'!\n", pos);
		return FALSE;
	}
	if( !cvt_write_reg( &p_data->m_common_data, CVT_V792_MCST_CBLT_CTRL_INDEX, &reg_value))
	{
		TRACE( "V792 CVT_V792_MCST_CBLT_CTRL write failed !\n");
		return FALSE;
	}

	return TRUE;
}
////////////////////////////////////////////////////////////////////////////////////////////////
//
//     M I S C E L L A N E O U S   A P I s
//
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////////////////////
const char* cvt_V792_SW_rev( void)
{
	return "CAENqtpSDK Rev. 1.1";
}

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \file    cvt_V792.h
*   \brief   V792 VME board definitions
*   \author  NDA
*   \version 1.0
*   \date    03/2006
*            
*            Provides methods, properties and defines to handle V792 VME boards
*/
////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef __CVT_V792_DEF_H
#define __CVT_V792_DEF_H

////////////////////////////////////////////////////////////////////////////////////////////////
// File includes
////////////////////////////////////////////////////////////////////////////////////////////////
#include "cvt_common_defs.h"
#include "cvt_board_commons.h"
////////////////////////////////////////////////////////////////////////////////////////////////
// Global defines
////////////////////////////////////////////////////////////////////////////////////////////////

#define CVT_V792_THRESHOLD_NUM_A		32			/*!< \brief Number of threshold registers for V792AA or V792AC */
#define CVT_V792_THRESHOLD_NUM_N		16			/*!< \brief Number of threshold registers for V792NA or V792NC */

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \def     CVT_V792_USE_DATA_QUEUE
*   \brief   Enable/ disable V792 builtin data queue
*            
*            Define this to enable V792 builtin data queue; comment out to disable this feature
*/
////////////////////////////////////////////////////////////////////////////////////////////////
#define CVT_V792_USE_DATA_QUEUE	1

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \def     CVT_V792_QUEUE_SIZE_DWORD
*   \brief   V792 data queue size 
*            
*            This is the number of UINT32 words storage into the V792 builtin data queue
*/
////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef CVT_V792_USE_DATA_QUEUE
	#define CVT_V792_QUEUE_SIZE_DWORD	(1024*1024)	
#endif

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \enum    CVT_QTP_TYPES
*   \brief   The QTP board type
*            
*            Provides an entry foreach supported board type
*/
////////////////////////////////////////////////////////////////////////////////////////////////
typedef enum
{
	CVT_V792_TYPE_A,							/*!< \brief The board is V792A */
	CVT_V792_TYPE_N,							/*!< \brief The board is V792N */
	CVT_V862_TYPE_A,							/*!< \brief The board is V862A */
	CVT_V965,									/*!< \brief The board is V965  */
	CVT_V965_TYPE_A,							/*!< \brief The board is V965A */
	CVT_V785_TYPE_A,							/*!< \brief The board is V785A */
	CVT_V785_TYPE_N,							/*!< \brief The board is V785N */
	CVT_V1785,									/*!< \brief The board is V1785 */
	CVT_V775_TYPE_A,							/*!< \brief The board is V775A */
	CVT_V775_TYPE_N,							/*!< \brief The board is V775N */
} CVT_QTP_TYPES;

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \struct  cvt_V792_data
*   \brief   V792 board data structure
*            
*            This structure stores any board data, common and specific ones. 
*            It must be passed to each V792's API 
*/
////////////////////////////////////////////////////////////////////////////////////////////////
typedef struct
{
	cvt_board_data	m_common_data;				/*!< Boards common data structure: always let at offset 0 */
	//
	// Board specific data

	CVT_QTP_TYPES m_type;						/*!< The QTP board type */
#ifdef CVT_V792_USE_DATA_QUEUE
	UINT32 *m_queue;							/*!< Builtin read data queue: only if CVT_V792_USE_DATA_QUEUE is defined */
	long m_queue_ini;							/*!< Builtin read data queue bottom offset: valid data can be peeked from queue_ini to queue_end. (only if CVT_V792_USE_DATA_QUEUE is defined) */
	long m_queue_end;							/*!< Builtin read data queue top offset: new data can be stored from queue_end to queue_ini (only if CVT_V792_USE_DATA_QUEUE is defined) */
#endif

} cvt_V792_data;

////////////////////////////////////////////////////////////////////////////////////////////////
// Registers address
////////////////////////////////////////////////////////////////////////////////////////////////
#define CVT_V792_OUT_BUFFER_ADD						0x0000				/*!< \brief Output buffer relative address */
#define CVT_V792_FW_REV_ADD							0x1000				/*!< \brief Firmware Revision register relative address */
#define CVT_V792_GEO_ADDRESS_ADD					0x1002				/*!< \brief Geo Address register relative address */
#define CVT_V792_MCST_CBLT_ADDRESS_ADD				0x1004				/*!< \brief MCST/CBLT Address register relative address */
#define CVT_V792_BIT_SET_1_ADD						0x1006				/*!< \brief Bit Set 1 register relative address */
#define CVT_V792_BIT_CLEAR_1_ADD					0x1008				/*!< \brief Bit Clear 1 register relative address */
#define CVT_V792_INT_LEVEL_ADD						0x100A				/*!< \brief Interrupt level register relative address */
#define CVT_V792_INT_VECTOR_ADD						0x100C				/*!< \brief Interrupt vector register relative address */
#define CVT_V792_STATUS_1_ADD						0x100E				/*!< \brief Status 1 register relative address */
#define CVT_V792_CONTROL_1_ADD						0x1010				/*!< \brief Control 1 register relative address */
#define CVT_V792_ADER_HIGH_ADD						0x1012				/*!< \brief Address decoder HIGH register relative address */
#define CVT_V792_ADER_LOW_ADD						0x1014				/*!< \brief Address decoder LOW register relative address */
#define CVT_V792_SINGLE_SHOT_RESET_ADD				0x1016				/*!< \brief Single Shot Reset register relative address */
#define CVT_V792_MCST_CBLT_CTRL_ADD					0x101A				/*!< \brief MCST/CBLT control register relative address */
#define CVT_V792_EVENT_TRG_ADD						0x1020				/*!< \brief Event trigger register relative address */
#define CVT_V792_STATUS_2_ADD						0x1022				/*!< \brief Status 2 register relative address */
#define CVT_V792_EVENT_COUNTER_LOW_ADD				0x1024				/*!< \brief Event counter low register relative address */
#define CVT_V792_EVENT_COUNTER_HIGH_ADD				0x1026				/*!< \brief Event counter high register relative address */
#define CVT_V792_INC_EVENT_ADD						0x1028				/*!< \brief Increment Event register relative address */
#define CVT_V792_INC_OFFSET_ADD						0x102A				/*!< \brief Increment Offset register relative address */
#define CVT_V792_LOAD_TEST_ADD						0x102C				/*!< \brief Load Test register relative address */
#define CVT_V792_FCLR_WND_ADD						0x102E				/*!< \brief FCLR window register relative address */
#define CVT_V792_BIT_SET_2_ADD						0x1032				/*!< \brief Bit Set 2 register relative address */
#define CVT_V792_BIT_CLEAR_2_ADD					0x1034				/*!< \brief Bit Clear 2 register relative address */
#define CVT_V792_W_MEMORY_TEST_ADDRESS_ADD			0x1036				/*!< \brief W Memory Test Address register relative address */
#define CVT_V792_MEMORY_TEST_WORD_HIGH_ADD			0x1038				/*!< \brief Memory Test Word High register relative address */
#define CVT_V792_MEMORY_TEST_WORD_LOW_ADD			0x103A				/*!< \brief Memory Test Word Low register relative address */
#define CVT_V792_CRATE_SELECT_ADD					0x103C				/*!< \brief Crate select register relative address */
#define CVT_V792_TEST_EVENT_WRITE_ADD				0x103E				/*!< \brief Test Event register Write relative address */
#define CVT_V792_EVENT_COUNT_RESET_ADD				0x1040				/*!< \brief Event Counter Reset register relative address */
#define CVT_V792_IPED_ADD							0x1060				/*!< \brief Iped register relative address */
#define CVT_V792_R_TEST_ADDRESS_ADD					0x1064				/*!< \brief R Test Address register relative address */
#define CVT_V792_SW_COMM_ADD						0x1068				/*!< \brief SW Comm register relative address */
#define CVT_V792_SLIDE_CONSTANT_ADD					0x106A				/*!< \brief Sliding constant register relative address */
#define CVT_V792_AAD_ADD							0x1070				/*!< \brief AAD register relative address */
#define CVT_V792_BAD_ADD							0x1072				/*!< \brief BAD register relative address */
#define CVT_V792_THRESHOLD_0_ADD					0x1080				/*!< \brief Threshold 0 register relative address */
#define CVT_V792_THRESHOLD_1_ADD					0x1082				/*!< \brief Threshold 1 register relative address */
#define CVT_V792_THRESHOLD_2_ADD					0x1084				/*!< \brief Threshold 2 register relative address */
#define CVT_V792_THRESHOLD_3_ADD					0x1086				/*!< \brief Threshold 3 register relative address */
#define CVT_V792_THRESHOLD_4_ADD					0x1088				/*!< \brief Threshold 4 register relative address */
#define CVT_V792_THRESHOLD_5_ADD					0x108A				/*!< \brief Threshold 5 register relative address */
#define CVT_V792_THRESHOLD_6_ADD					0x108C				/*!< \brief Threshold 6 register relative address */
#define CVT_V792_THRESHOLD_7_ADD					0x108E				/*!< \brief Threshold 7 register relative address */
#define CVT_V792_THRESHOLD_8_ADD					0x1090				/*!< \brief Threshold 8 register relative address */
#define CVT_V792_THRESHOLD_9_ADD					0x1092				/*!< \brief Threshold 9 register relative address */
#define CVT_V792_THRESHOLD_10_ADD					0x1094				/*!< \brief Threshold 10 register relative address */
#define CVT_V792_THRESHOLD_11_ADD					0x1096				/*!< \brief Threshold 11 register relative address */
#define CVT_V792_THRESHOLD_12_ADD					0x1098				/*!< \brief Threshold 12 register relative address */
#define CVT_V792_THRESHOLD_13_ADD					0x109A				/*!< \brief Threshold 13 register relative address */
#define CVT_V792_THRESHOLD_14_ADD					0x109C				/*!< \brief Threshold 14 register relative address */
#define CVT_V792_THRESHOLD_15_ADD					0x109E				/*!< \brief Threshold 15 register relative address */
#define CVT_V792_THRESHOLD_16_ADD					0x10A0				/*!< \brief Threshold 16 register relative address */
#define CVT_V792_THRESHOLD_17_ADD					0x10A2				/*!< \brief Threshold 17 register relative address */
#define CVT_V792_THRESHOLD_18_ADD					0x10A4				/*!< \brief Threshold 18 register relative address */
#define CVT_V792_THRESHOLD_19_ADD					0x10A6				/*!< \brief Threshold 19 register relative address */
#define CVT_V792_THRESHOLD_20_ADD					0x10A8				/*!< \brief Threshold 20 register relative address */
#define CVT_V792_THRESHOLD_21_ADD					0x10AA				/*!< \brief Threshold 21 register relative address */
#define CVT_V792_THRESHOLD_22_ADD					0x10AC				/*!< \brief Threshold 22 register relative address */
#define CVT_V792_THRESHOLD_23_ADD					0x10AE				/*!< \brief Threshold 23 register relative address */
#define CVT_V792_THRESHOLD_24_ADD					0x10B0				/*!< \brief Threshold 24 register relative address */
#define CVT_V792_THRESHOLD_25_ADD					0x10B2				/*!< \brief Threshold 25 register relative address */
#define CVT_V792_THRESHOLD_26_ADD					0x10B4				/*!< \brief Threshold 26 register relative address */
#define CVT_V792_THRESHOLD_27_ADD					0x10B6				/*!< \brief Threshold 27 register relative address */
#define CVT_V792_THRESHOLD_28_ADD					0x10B8				/*!< \brief Threshold 28 register relative address */
#define CVT_V792_THRESHOLD_29_ADD					0x10BA				/*!< \brief Threshold 29 register relative address */
#define CVT_V792_THRESHOLD_30_ADD					0x10BC				/*!< \brief Threshold 30 register relative address */
#define CVT_V792_THRESHOLD_31_ADD					0x10BE				/*!< \brief Threshold 31 register relative address */
#define CVT_V792_ROM_OUI_2_ADD						0x8026				/*!< \brief Manufacturer identifier (IEEE OUI) (MSB) relative address */
#define CVT_V792_ROM_OUI_1_ADD						0x802A				/*!< \brief Manufacturer identifier (IEEE OUI) relative address */
#define CVT_V792_ROM_OUI_0_ADD						0x802E				/*!< \brief Manufacturer identifier (IEEE OUI) (LSB) relative address */
#define CVT_V792_ROM_VERSION_ADD					0x8032				/*!< \brief Purchased version of the Mod.V792 relative address */
#define CVT_V792_ROM_BOARD_ID_2_ADD					0x8036				/*!< \brief Board identifier (MSB) relative address */
#define CVT_V792_ROM_BOARD_ID_1_ADD					0x803A				/*!< \brief Board identifier relative address */
#define CVT_V792_ROM_BOARD_ID_0_ADD					0x803E				/*!< \brief Board identifier (LSB) relative address */
#define CVT_V792_ROM_REVISION_ADD					0x804E				/*!< \brief Hardware revision identifier relative address */
#define CVT_V792_ROM_SERIAL_1_ADD					0x8F02				/*!< \brief Serial number (MSB) relative address */
#define CVT_V792_ROM_SERIAL_0_ADD					0x8F06				/*!< \brief Serial number (LSB) relative address */

////////////////////////////////////////////////////////////////////////////////////////////////
// Registers data size
////////////////////////////////////////////////////////////////////////////////////////////////
#define CVT_V792_OUT_BUFFER_DATA_SIZE				cvD64				/*!< \brief Output buffer data size */
#define CVT_V792_FW_REV_DATA_SIZE					cvD16				/*!< \brief Firmware Revision register data size */
#define CVT_V792_GEO_ADDRESS_DATA_SIZE				cvD16				/*!< \brief Geo Address register data size */
#define CVT_V792_MCST_CBLT_ADDRESS_DATA_SIZE		cvD16				/*!< \brief MCST/CBLT Address register data size */
#define CVT_V792_BIT_SET_1_DATA_SIZE				cvD16				/*!< \brief Bit Set 1 register data size */
#define CVT_V792_BIT_CLEAR_1_DATA_SIZE				cvD16				/*!< \brief Bit Clear 1 register data size */
#define CVT_V792_INT_LEVEL_DATA_SIZE				cvD16				/*!< \brief Interrupt level register data size */
#define CVT_V792_INT_VECTOR_DATA_SIZE				cvD16				/*!< \brief Interrupt vector register data size */
#define CVT_V792_STATUS_1_DATA_SIZE					cvD16				/*!< \brief Status 1 register data size */
#define CVT_V792_CONTROL_1_DATA_SIZE				cvD16				/*!< \brief Control 1 register data size */
#define CVT_V792_ADER_HIGH_DATA_SIZE				cvD16				/*!< \brief Address decoder HIGH register data size */
#define CVT_V792_ADER_LOW_DATA_SIZE					cvD16				/*!< \brief Address decoder LOW register data size */
#define CVT_V792_SINGLE_SHOT_RESET_DATA_SIZE		cvD16				/*!< \brief Single Shot Reset register data size */
#define CVT_V792_MCST_CBLT_CTRL_DATA_SIZE			cvD16				/*!< \brief MCST/CBLT control register data size */
#define CVT_V792_EVENT_TRG_DATA_SIZE				cvD16				/*!< \brief Event trigger register data size */
#define CVT_V792_STATUS_2_DATA_SIZE					cvD16				/*!< \brief Status 2 register data size */
#define CVT_V792_EVENT_COUNTER_LOW_DATA_SIZE		cvD16				/*!< \brief Event counter low register data size */
#define CVT_V792_EVENT_COUNTER_HIGH_DATA_SIZE		cvD16				/*!< \brief Event counter high register data size */
#define CVT_V792_INC_EVENT_DATA_SIZE				cvD16				/*!< \brief Increment Event register data size */
#define CVT_V792_INC_OFFSET_DATA_SIZE				cvD16				/*!< \brief Increment Offset register data size */
#define CVT_V792_LOAD_TEST_DATA_SIZE				cvD16				/*!< \brief Load Test register data size */
#define CVT_V792_FCLR_WND_DATA_SIZE					cvD16				/*!< \brief FCLR window register data size */
#define CVT_V792_BIT_SET_2_DATA_SIZE				cvD16				/*!< \brief Bit Set 2 register data size */
#define CVT_V792_BIT_CLEAR_2_DATA_SIZE				cvD16				/*!< \brief Bit Clear 2 register data size */
#define CVT_V792_W_MEMORY_TEST_ADDRESS_DATA_SIZE	cvD16				/*!< \brief W Memory Test Address register data size */
#define CVT_V792_MEMORY_TEST_WORD_HIGH_DATA_SIZE	cvD16				/*!< \brief Memory Test Word High register data size */
#define CVT_V792_MEMORY_TEST_WORD_LOW_DATA_SIZE		cvD16				/*!< \brief Memory Test Word Low register data size */
#define CVT_V792_CRATE_SELECT_DATA_SIZE				cvD16				/*!< \brief Crate select register data size */
#define CVT_V792_TEST_EVENT_WRITE_DATA_SIZE			cvD16				/*!< \brief Test Event Write register data size */
#define CVT_V792_EVENT_COUNT_RESET_DATA_SIZE		cvD16				/*!< \brief Event Counter Reset register data size */
#define CVT_V792_IPED_DATA_SIZE						cvD16				/*!< \brief Iped register data size */
#define CVT_V792_R_TEST_ADDRESS_DATA_SIZE			cvD16				/*!< \brief R Test Address register data size */
#define CVT_V792_SW_COMM_DATA_SIZE					cvD16				/*!< \brief SW Comm register data size */
#define CVT_V792_SLIDE_CONSTANT_DATA_SIZE			cvD16				/*!< \brief Sliding constant register data size */
#define CVT_V792_AAD_DATA_SIZE						cvD16				/*!< \brief AAD register data size */
#define CVT_V792_BAD_DATA_SIZE						cvD16				/*!< \brief BAD register data size */
#define CVT_V792_THRESHOLD_DATA_SIZE				cvD16				/*!< \brief Threshold register data size */
#define CVT_V792_ROM_OUI_2_DATA_SIZE				cvD16				/*!< \brief Manufacturer identifier (IEEE OUI) (MSB) data size */
#define CVT_V792_ROM_OUI_1_DATA_SIZE				cvD16				/*!< \brief Manufacturer identifier (IEEE OUI) data size */
#define CVT_V792_ROM_OUI_0_DATA_SIZE				cvD16				/*!< \brief Manufacturer identifier (IEEE OUI) (LSB) data size */
#define CVT_V792_ROM_VERSION_DATA_SIZE				cvD16				/*!< \brief Purchased version of the Mod.V792 data size */
#define CVT_V792_ROM_BOARD_ID_2_DATA_SIZE			cvD16				/*!< \brief Board identifier (MSB) data size */
#define CVT_V792_ROM_BOARD_ID_1_DATA_SIZE			cvD16				/*!< \brief Board identifier data size */
#define CVT_V792_ROM_BOARD_ID_0_DATA_SIZE			cvD16				/*!< \brief Board identifier (LSB) data size */
#define CVT_V792_ROM_REVISION_DATA_SIZE				cvD16				/*!< \brief Hardware revision identifier data size */
#define CVT_V792_ROM_SERIAL_1_DATA_SIZE				cvD16				/*!< \brief Serial number (MSB) data size */
#define CVT_V792_ROM_SERIAL_0_DATA_SIZE				cvD16				/*!< \brief Serial number (LSB) data size */

////////////////////////////////////////////////////////////////////////////////////////////////
// Registers address modifiers
////////////////////////////////////////////////////////////////////////////////////////////////
#define CVT_V792_OUT_BUFFER_AM						cvA32_S_MBLT	/*!< \brief Output buffer address modifier */
#define CVT_V792_FW_REV_AM							cvA32_S_DATA	/*!< \brief Firmware Revision register address modifier */
#define CVT_V792_GEO_ADDRESS_AM						cvA32_S_DATA	/*!< \brief Geo Address register address modifier */
#define CVT_V792_MCST_CBLT_ADDRESS_AM				cvA32_S_DATA	/*!< \brief MCST/CBLT Address register address modifier */
#define CVT_V792_BIT_SET_1_AM						cvA32_S_DATA	/*!< \brief Bit Set 1 register address modifier */
#define CVT_V792_BIT_CLEAR_1_AM						cvA32_S_DATA	/*!< \brief Bit Clear 1 register address modifier */
#define CVT_V792_INT_LEVEL_AM						cvA32_S_DATA	/*!< \brief Interrupt level register address modifier */
#define CVT_V792_INT_VECTOR_AM						cvA32_S_DATA	/*!< \brief Interrupt vector register address modifier */
#define CVT_V792_STATUS_1_AM						cvA32_S_DATA	/*!< \brief Status 1 register address modifier */
#define CVT_V792_CONTROL_1_AM						cvA32_S_DATA	/*!< \brief Control 1 register address modifier */
#define CVT_V792_ADER_HIGH_AM						cvA32_S_DATA	/*!< \brief Address decoder HIGH register address modifier */
#define CVT_V792_ADER_LOW_AM						cvA32_S_DATA	/*!< \brief Address decoder LOW register address modifier */
#define CVT_V792_SINGLE_SHOT_RESET_AM				cvA32_S_DATA	/*!< \brief Single Shot Reset register address modifier */
#define CVT_V792_MCST_CBLT_CTRL_AM					cvA32_S_DATA	/*!< \brief MCST/CBLT control register address modifier */
#define CVT_V792_EVENT_TRG_AM						cvA32_S_DATA	/*!< \brief Event trigger register address modifier */
#define CVT_V792_STATUS_2_AM						cvA32_S_DATA	/*!< \brief Status 2 register address modifier */
#define CVT_V792_EVENT_COUNTER_LOW_AM				cvA32_S_DATA	/*!< \brief Event counter low register address modifier */
#define CVT_V792_EVENT_COUNTER_HIGH_AM				cvA32_S_DATA	/*!< \brief Event counter high register address modifier */
#define CVT_V792_INC_EVENT_AM						cvA32_S_DATA	/*!< \brief Increment Event register address modifier */
#define CVT_V792_INC_OFFSET_AM						cvA32_S_DATA	/*!< \brief Increment Offset register address modifier */
#define CVT_V792_LOAD_TEST_AM						cvA32_S_DATA	/*!< \brief Load Test register address modifier */
#define CVT_V792_FCLR_WND_AM						cvA32_S_DATA	/*!< \brief FCLR window register address modifier */
#define CVT_V792_BIT_SET_2_AM						cvA32_S_DATA	/*!< \brief Bit Set 2 register address modifier */
#define CVT_V792_BIT_CLEAR_2_AM						cvA32_S_DATA	/*!< \brief Bit Clear 2 register address modifier */
#define CVT_V792_W_MEMORY_TEST_ADDRESS_AM			cvA32_S_DATA	/*!< \brief W Memory Test Address register address modifier */
#define CVT_V792_MEMORY_TEST_WORD_HIGH_AM			cvA32_S_DATA	/*!< \brief Memory Test Word High register address modifier */
#define CVT_V792_MEMORY_TEST_WORD_LOW_AM			cvA32_S_DATA	/*!< \brief Memory Test Word Low register address modifier */
#define CVT_V792_CRATE_SELECT_AM					cvA32_S_DATA	/*!< \brief Crate select register address modifier */
#define CVT_V792_TEST_EVENT_WRITE_AM				cvA32_S_DATA	/*!< \brief Test Event Write register address modifier */
#define CVT_V792_EVENT_COUNT_RESET_AM				cvA32_S_DATA	/*!< \brief Event Counter Reset register address modifier */
#define CVT_V792_IPED_AM							cvA32_S_DATA	/*!< \brief Iped register address modifier */
#define CVT_V792_R_TEST_ADDRESS_AM					cvA32_S_DATA	/*!< \brief R Test Address register address modifier */
#define CVT_V792_SW_COMM_AM							cvA32_S_DATA	/*!< \brief SW Comm register address modifier */
#define CVT_V792_SLIDE_CONSTANT_AM					cvA32_S_DATA	/*!< \brief Sliding constant register address modifier */
#define CVT_V792_AAD_AM								cvA32_S_DATA	/*!< \brief AAD register address modifier */
#define CVT_V792_BAD_AM								cvA32_S_DATA	/*!< \brief BAD register address modifier */
#define CVT_V792_THRESHOLD_AM						cvA32_S_DATA	/*!< \brief Threshold register address modifier */
#define CVT_V792_ROM_OUI_2_AM						cvA32_S_DATA	/*!< \brief Manufacturer identifier (IEEE OUI) (MSB) address modifier */
#define CVT_V792_ROM_OUI_1_AM						cvA32_S_DATA	/*!< \brief Manufacturer identifier (IEEE OUI) address modifier */
#define CVT_V792_ROM_OUI_0_AM						cvA32_S_DATA	/*!< \brief Manufacturer identifier (IEEE OUI) (LSB) address modifier */
#define CVT_V792_ROM_VERSION_AM						cvA32_S_DATA	/*!< \brief Purchased version of the Mod.V792 address modifier */
#define CVT_V792_ROM_BOARD_ID_2_AM					cvA32_S_DATA	/*!< \brief Board identifier (MSB) address modifier */
#define CVT_V792_ROM_BOARD_ID_1_AM					cvA32_S_DATA	/*!< \brief Board identifier address modifier */
#define CVT_V792_ROM_BOARD_ID_0_AM					cvA32_S_DATA	/*!< \brief Board identifier (LSB) address modifier */
#define CVT_V792_ROM_REVISION_AM					cvA32_S_DATA	/*!< \brief Hardware revision identifier address modifier */
#define CVT_V792_ROM_SERIAL_1_AM					cvA32_S_DATA	/*!< \brief Serial number (MSB) address modifier */
#define CVT_V792_ROM_SERIAL_0_AM					cvA32_S_DATA	/*!< \brief Serial number (LSB) address modifier */

////////////////////////////////////////////////////////////////////////////////////////////////
// Registers indexes
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \enum    CVT_V792_REG_INDEX
*   \brief   The registers indexes
*            
*            Provides an entry for each register: This is the index into the CVT_V792_REG_TABLE board table
*/
////////////////////////////////////////////////////////////////////////////////////////////////
typedef enum
{
	CVT_V792_OUT_BUFFER_INDEX,					/*!< \brief Output buffer index */
	CVT_V792_FW_REV_INDEX,						/*!< \brief Firmware Revision register index */
	CVT_V792_GEO_ADDRESS_INDEX,					/*!< \brief Geo Address register index */
	CVT_V792_MCST_CBLT_ADDRESS_INDEX,			/*!< \brief MCST/CBLT Address register index */
	CVT_V792_BIT_SET_1_INDEX,					/*!< \brief Bit Set 1 register index */
	CVT_V792_BIT_CLEAR_1_INDEX,					/*!< \brief Bit Clear 1 register index */
	CVT_V792_INT_LEVEL_INDEX,					/*!< \brief Interrupt level register index */
	CVT_V792_INT_VECTOR_INDEX,					/*!< \brief Interrupt vector register index */
	CVT_V792_STATUS_1_INDEX,					/*!< \brief Status 1 register index */
	CVT_V792_CONTROL_1_INDEX,					/*!< \brief Control 1 register index */
	CVT_V792_ADER_HIGH_INDEX,					/*!< \brief Address decoder HIGH register index */
	CVT_V792_ADER_LOW_INDEX,					/*!< \brief Address decoder LOW register index */
	CVT_V792_SINGLE_SHOT_RESET_INDEX,			/*!< \brief Single Shot Reset register index */
	CVT_V792_MCST_CBLT_CTRL_INDEX,				/*!< \brief MCST/CBLT control register index */
	CVT_V792_EVENT_TRG_INDEX,					/*!< \brief Event trigger register index */
	CVT_V792_STATUS_2_INDEX,					/*!< \brief Status 2 register index */
	CVT_V792_EVENT_COUNTER_LOW_INDEX,			/*!< \brief Event counter low register index */
	CVT_V792_EVENT_COUNTER_HIGH_INDEX,			/*!< \brief Event counter high register index */
	CVT_V792_INC_EVENT_INDEX,					/*!< \brief Increment Event register index */
	CVT_V792_INC_OFFSET_INDEX,					/*!< \brief Increment Offset register index */
	CVT_V792_LOAD_TEST_INDEX,					/*!< \brief Load Test register index */
	CVT_V792_FCLR_WND_INDEX,					/*!< \brief FCLR window register index */
	CVT_V792_BIT_SET_2_INDEX,					/*!< \brief Bit Set 2 register index */
	CVT_V792_BIT_CLEAR_2_INDEX,					/*!< \brief Bit Clear 2 register index */
	CVT_V792_W_MEMORY_TEST_ADDRESS_INDEX,		/*!< \brief W Memory Test Address register index */
	CVT_V792_MEMORY_TEST_WORD_HIGH_INDEX,		/*!< \brief Memory Test Word High register index */
	CVT_V792_MEMORY_TEST_WORD_LOW_INDEX,		/*!< \brief Memory Test Word Low register index */
	CVT_V792_CRATE_SELECT_INDEX,				/*!< \brief Crate select register index */
	CVT_V792_TEST_EVENT_WRITE_INDEX,			/*!< \brief Test Event register Write index */
	CVT_V792_EVENT_COUNT_RESET_INDEX,			/*!< \brief Event Counter Reset register index */
	CVT_V792_IPED_INDEX,						/*!< \brief Iped register index */
	CVT_V792_R_TEST_ADDRESS_INDEX,				/*!< \brief R Test Address register index */
	CVT_V792_SW_COMM_INDEX,						/*!< \brief SW Comm register index */
	CVT_V792_SLIDE_CONSTANT_INDEX,				/*!< \brief Sliding constant register index */
	CVT_V792_AAD_INDEX,							/*!< \brief AAD register index */
	CVT_V792_BAD_INDEX,							/*!< \brief BAD register index */
	CVT_V792_THRESHOLD_0_INDEX,					/*!< \brief Threshold 0 register index  ( Ch  0 -> V792A Ch  0 -> V792N )*/
	CVT_V792_THRESHOLD_1_INDEX,					/*!< \brief Threshold 1 register index  ( Ch  1 -> V792A )*/
	CVT_V792_THRESHOLD_2_INDEX,					/*!< \brief Threshold 2 register index  ( Ch  2 -> V792A Ch  1 -> V792N )*/
	CVT_V792_THRESHOLD_3_INDEX,					/*!< \brief Threshold 3 register index  ( Ch  3 -> V792A )*/
	CVT_V792_THRESHOLD_4_INDEX,					/*!< \brief Threshold 4 register index  ( Ch  4 -> V792A Ch  2 -> V792N )*/
	CVT_V792_THRESHOLD_5_INDEX,					/*!< \brief Threshold 5 register index  ( Ch  5 -> V792A )*/
	CVT_V792_THRESHOLD_6_INDEX,					/*!< \brief Threshold 6 register index  ( Ch  6 -> V792A Ch  3 -> V792N )*/
	CVT_V792_THRESHOLD_7_INDEX,					/*!< \brief Threshold 7 register index  ( Ch  7 -> V792A )*/
	CVT_V792_THRESHOLD_8_INDEX,					/*!< \brief Threshold 8 register index  ( Ch  8 -> V792A Ch  4 -> V792N )*/
	CVT_V792_THRESHOLD_9_INDEX,					/*!< \brief Threshold 9 register index  ( Ch  9 -> V792A )*/
	CVT_V792_THRESHOLD_10_INDEX,				/*!< \brief Threshold 10 register index ( Ch 10 -> V792A Ch  5 -> V792N )*/
	CVT_V792_THRESHOLD_11_INDEX,				/*!< \brief Threshold 11 register index ( Ch 11 -> V792A )*/
	CVT_V792_THRESHOLD_12_INDEX,				/*!< \brief Threshold 12 register index ( Ch 12 -> V792A Ch  6 -> V792N )*/
	CVT_V792_THRESHOLD_13_INDEX,				/*!< \brief Threshold 13 register index ( Ch 13 -> V792A )*/
	CVT_V792_THRESHOLD_14_INDEX,				/*!< \brief Threshold 14 register index ( Ch 14 -> V792A Ch  7 -> V792N )*/
	CVT_V792_THRESHOLD_15_INDEX,				/*!< \brief Threshold 15 register index ( Ch 15 -> V792A )*/
	CVT_V792_THRESHOLD_16_INDEX,				/*!< \brief Threshold 16 register index ( Ch 16 -> V792A Ch  8 -> V792N )*/
	CVT_V792_THRESHOLD_17_INDEX,				/*!< \brief Threshold 17 register index ( Ch 17 -> V792A )*/
	CVT_V792_THRESHOLD_18_INDEX,				/*!< \brief Threshold 18 register index ( Ch 18 -> V792A Ch  9 -> V792N )*/
	CVT_V792_THRESHOLD_19_INDEX,				/*!< \brief Threshold 19 register index ( Ch 19 -> V792A )*/
	CVT_V792_THRESHOLD_20_INDEX,				/*!< \brief Threshold 20 register index ( Ch 20 -> V792A Ch 10 -> V792N )*/
	CVT_V792_THRESHOLD_21_INDEX,				/*!< \brief Threshold 21 register index ( Ch 21 -> V792A )*/
	CVT_V792_THRESHOLD_22_INDEX,				/*!< \brief Threshold 22 register index ( Ch 22 -> V792A Ch 11 -> V792N )*/
	CVT_V792_THRESHOLD_23_INDEX,				/*!< \brief Threshold 23 register index ( Ch 23 -> V792A )*/
	CVT_V792_THRESHOLD_24_INDEX,				/*!< \brief Threshold 24 register index ( Ch 24 -> V792A Ch 12 -> V792N )*/
	CVT_V792_THRESHOLD_25_INDEX,				/*!< \brief Threshold 25 register index ( Ch 25 -> V792A )*/
	CVT_V792_THRESHOLD_26_INDEX,				/*!< \brief Threshold 26 register index ( Ch 26 -> V792A Ch 13 -> V792N )*/
	CVT_V792_THRESHOLD_27_INDEX,				/*!< \brief Threshold 27 register index ( Ch 27 -> V792A )*/
	CVT_V792_THRESHOLD_28_INDEX,				/*!< \brief Threshold 28 register index ( Ch 28 -> V792A Ch 14 -> V792N )*/
	CVT_V792_THRESHOLD_29_INDEX,				/*!< \brief Threshold 29 register index ( Ch 29 -> V792A )*/
	CVT_V792_THRESHOLD_30_INDEX,				/*!< \brief Threshold 30 register index ( Ch 30 -> V792A Ch 15 -> V792N )*/
	CVT_V792_THRESHOLD_31_INDEX,				/*!< \brief Threshold 31 register index ( Ch 31 -> V792A )*/
	CVT_V792_ROM_OUI_2_INDEX,					/*!< \brief Manufacturer identifier (IEEE OUI) (MSB) index */
	CVT_V792_ROM_OUI_1_INDEX,					/*!< \brief Manufacturer identifier (IEEE OUI) index */
	CVT_V792_ROM_OUI_0_INDEX,					/*!< \brief Manufacturer identifier (IEEE OUI) (LSB) index */
	CVT_V792_ROM_VERSION_INDEX,					/*!< \brief Purchased version of the Mod.V792 index */
	CVT_V792_ROM_BOARD_ID_2_INDEX,				/*!< \brief Board identifier (MSB) index */
	CVT_V792_ROM_BOARD_ID_1_INDEX,				/*!< \brief Board identifier index */
	CVT_V792_ROM_BOARD_ID_0_INDEX,				/*!< \brief Board identifier (LSB) index */
	CVT_V792_ROM_REVISION_INDEX,				/*!< \brief Hardware revision identifier index */
	CVT_V792_ROM_SERIAL_1_INDEX,				/*!< \brief Serial number (MSB) index */
	CVT_V792_ROM_SERIAL_0_INDEX,				/*!< \brief Serial number (LSB) index */
} CVT_V792_REG_INDEX;

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \enum    cvt_V792_threshold_id 
*   \brief   V792 discriminator threshold IDs
*/
////////////////////////////////////////////////////////////////////////////////////////////////
typedef enum 
{
	CVT_V792_THRESHOLD_ALL= -1,						/*!< Catch all the thresholds */
	CVT_V792_THRESHOLD_0= 0,						/*!< Threshold 0 register */
	CVT_V792_THRESHOLD_1,							/*!< Threshold 1 register */
	CVT_V792_THRESHOLD_2,							/*!< Threshold 2 register */
	CVT_V792_THRESHOLD_3,							/*!< Threshold 3 register */
	CVT_V792_THRESHOLD_4,							/*!< Threshold 4 register */
	CVT_V792_THRESHOLD_5,							/*!< Threshold 5 register */
	CVT_V792_THRESHOLD_6,							/*!< Threshold 6 register */
	CVT_V792_THRESHOLD_7,							/*!< Threshold 7 register */
	CVT_V792_THRESHOLD_8,							/*!< Threshold 8 register */
	CVT_V792_THRESHOLD_9,							/*!< Threshold 9 register */
	CVT_V792_THRESHOLD_10,							/*!< Threshold 10 register */
	CVT_V792_THRESHOLD_11,							/*!< Threshold 11 register */
	CVT_V792_THRESHOLD_12,							/*!< Threshold 12 register */
	CVT_V792_THRESHOLD_13,							/*!< Threshold 13 register */
	CVT_V792_THRESHOLD_14,							/*!< Threshold 14 register */
	CVT_V792_THRESHOLD_15,							/*!< Threshold 15 register */
	CVT_V792_THRESHOLD_16,							/*!< Threshold 16 register */
	CVT_V792_THRESHOLD_17,							/*!< Threshold 17 register */
	CVT_V792_THRESHOLD_18,							/*!< Threshold 18 register */
	CVT_V792_THRESHOLD_19,							/*!< Threshold 19 register */
	CVT_V792_THRESHOLD_20,							/*!< Threshold 20 register */
	CVT_V792_THRESHOLD_21,							/*!< Threshold 21 register */
	CVT_V792_THRESHOLD_22,							/*!< Threshold 22 register */
	CVT_V792_THRESHOLD_23,							/*!< Threshold 23 register */
	CVT_V792_THRESHOLD_24,							/*!< Threshold 24 register */
	CVT_V792_THRESHOLD_25,							/*!< Threshold 25 register */
	CVT_V792_THRESHOLD_26,							/*!< Threshold 26 register */
	CVT_V792_THRESHOLD_27,							/*!< Threshold 27 register */
	CVT_V792_THRESHOLD_28,							/*!< Threshold 28 register */
	CVT_V792_THRESHOLD_29,							/*!< Threshold 29 register */
	CVT_V792_THRESHOLD_30,							/*!< Threshold 30 register */
	CVT_V792_THRESHOLD_31,							/*!< Threshold 31 register */
} cvt_V792_threshold_id;
												
////////////////////////////////////////////////////////////////////////////////////////////////
/*! \var     CVT_QTP_NUM_CHANNELS
*   \brief   The number of channels foreach specific board
*/
////////////////////////////////////////////////////////////////////////////////////////////////
extern const int CVT_QTP_NUM_CHANNELS[];

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \enum    CVT_V792_THRESHOLD_MSK
*   \brief   Threshold register bitmasks
*/
////////////////////////////////////////////////////////////////////////////////////////////////
typedef enum
{
	CVT_V792_THRESHOLD_VALUE_MSK    = 0x00FF,		/*!< Threshold data bitmask */
	CVT_V792_THRESHOLD_KILL_MSK     = 0x0100,		/*!< Threshold KILL bitmask */
} CVT_V792_THRESHOLD_MSK;

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \enum    CVT_V792_BIT_SET_CLEAR_1_MSK
*   \brief   V792 Bit Set\\Clear 1 register masks
*            
*/
////////////////////////////////////////////////////////////////////////////////////////////////
typedef enum
{
	CVT_V792_BSC1_BERR_FLAG_MSK		= 0x0008,		/*!< Bus Error Flag bit */
	CVT_V792_BSC1_SEL_ADDR_MSK		= 0x0010,		/*!< Select Address bit. */
	CVT_V792_BSC1_SOFT_RESET_MSK	= 0x0080,		/*!< Software RESET bit. */
} CVT_V792_BIT_SET_CLEAR_1_MSK;

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \enum    CVT_V792_STATUS_1_MSK
*   \brief   V792 Status 1 register bit masks
*            
*/
////////////////////////////////////////////////////////////////////////////////////////////////
typedef enum
{
	CVT_V792_STS1_DREADY_MSK		= 0x0001,		/*!< Data Ready */
	CVT_V792_STS1_GLOBAL_DREADY_MSK	= 0x0002,		/*!< At least one module has Data Ready. */
	CVT_V792_STS1_BUSY_MSK			= 0x0004,		/*!< Module Busy. */
	CVT_V792_STS1_GLOBAL_BUSY_MSK	= 0x0008,		/*!< At least one module is Busy. */
	CVT_V792_STS1_AMNESIA_MSK		= 0x0010,		/*!< GEO is not available from the JAUX. */
	CVT_V792_STS1_PURGED_MSK		= 0x0020,		/*!< The board is purged. */
	CVT_V792_STS1_TERM_ON_MSK		= 0x0040,		/*!< All Control Bus Terminations are ON. */
	CVT_V792_STS1_TERM_OFF_MSK		= 0x0080,		/*!< All Control Bus Terminations are OFF. */
	CVT_V792_STS1_EVREADY_MSK		= 0x0100,		/*!< Indicates that the number in the Event Trigger Register 
													     is greater than or equal to the number of events stored in the
														 memory and an interrupt request has been generated with interrupt
														 level different from 0 */
} CVT_V792_STATUS_1_MSK;

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \enum    CVT_V792_CONTROL_1_MSK
*   \brief   V792 Control 1 register bit masks
*            
*/
////////////////////////////////////////////////////////////////////////////////////////////////
typedef enum
{
	CVT_V792_CTRL1_BLKEND_MSK		= 0x0004,		/*!< End of Block bit */
	CVT_V792_CTRL1_PROG_RESET_MSK	= 0x0010,		/*!< Programmable reset mode bit */
	CVT_V792_CTRL1_BERR_ENABLE_MSK	= 0x0020,		/*!< Bus error enable bit */
	CVT_V792_CTRL1_ALIGN64_MSK		= 0x0040,		/*!< Align 64 bit */
} CVT_V792_CONTROL_1_MSK;

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \enum    CVT_V792_MCST_CBLT_CTRL_MSK
*   \brief   V792 CVT_V792/MCST Control register bit masks
*            
*/
////////////////////////////////////////////////////////////////////////////////////////////////
typedef enum
{
	CVT_V792_MCCTRL_DISABLED_BOARD_MSK	= 0x0000,		/*!< Disabled Board mask */
	CVT_V792_MCCTRL_LAST_BOARD_MSK		= 0x0001,		/*!< Last Board mask */
	CVT_V792_MCCTRL_FIRST_BOARD_MSK		= 0x0002,		/*!< First Board mask */
	CVT_V792_MCCTRL_MID_BOARD_MSK		= 0x0003,		/*!< Middle Board mask */
} CVT_V792_MCST_CBLT_CTRL_MSK;

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \enum    CVT_V792_STATUS_2_MSK
*   \brief   V792 Status 2 register bit masks
*            
*/
////////////////////////////////////////////////////////////////////////////////////////////////
typedef enum
{
	CVT_V792_STS2_BUFFER_EMPTY_MSK		= 0x0002,		/*!< Buffer empty. */
	CVT_V792_STS2_BUFFER_FULL_MSK		= 0x0004,		/*!< Buffer full. */
	CVT_V792_STS2_DSEL_0_MSK			= 0x0010,		/*!< Indicate the type of piggy-back plugged into the board (bit 0). */
	CVT_V792_STS2_DSEL_1_MSK			= 0x0020,		/*!< Indicate the type of piggy-back plugged into the board (bit 1). */
	CVT_V792_STS2_CSEL_0_MSK			= 0x0040,		/*!< Indicate the type of piggy-back plugged into the board (bit 2). */
	CVT_V792_STS2_CSEL_1_MSK			= 0x0080,		/*!< Indicate the type of piggy-back plugged into the board (bit 3). */
} CVT_V792_STATUS_2_MSK;

#define CVT_V792_STS2_PIGGY_BACK_TYPE_MSK		0x00F0																												/*!< CVT_V792_STATUS_2 register piggy-back type bitmasks. */
#define CVT_V792_GET_PIGGY_BACK_TYPE( reg)		((((UINT16)reg)& CVT_V792_STS2_PIGGY_BACK_TYPE_MSK)>> 4)															/*!< \brief Extract the piggy back type from UINT16 value */
#define CVT_V792_SET_PIGGY_BACK_TYPE( reg)		reg= (((UINT16)reg)& ~CVT_V792_STS2_PIGGY_BACK_TYPE_MSK)| ((((UINT16)value)<< 4)&CVT_V792_STS2_PIGGY_BACK_TYPE_MSK) /*!< \brief Sets the piggy back type into UINT16 value */

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \enum    CVT_V792_BIT_SET_CLEAR_2_MSK
*   \brief   V792 Bit Set\\Clear 2 register masks
*/
////////////////////////////////////////////////////////////////////////////////////////////////
typedef enum
{
	CVT_V792_BSC2_TEST_MEM_MSK      = 0x0001,		/*!< Random Memory Access Test Mode bit */
	CVT_V792_BSC2_OFFLINE_MSK       = 0x0002,		/*!< Offline bit. */
	CVT_V792_BSC2_CLEAR_DATA_MSK    = 0x0004,		/*!< Clear Data bit. */
	CVT_V792_BSC2_OVER_RANGE_DIS_MSK= 0x0008,		/*!< Overflow suppression disable bit. */
	CVT_V792_BSC2_LOW_THR_DIS_MSK   = 0x0010,		/*!< Zero suppression disable bit. */
	CVT_V792_BSC2_VALID_DIS_MSK     = 0x0020,		/*!< Valid suppression disable bit. */
	CVT_V792_BSC2_TEST_ACQ_MSK      = 0x0040,		/*!< Acquisition Test Mode bit. */
	CVT_V792_BSC2_SLIDE_EN_MSK      = 0x0080,		/*!< Sliding scale enale/disable bit. */
	CVT_V792_BSC2_STEP_TH_MSK       = 0x0100,		/*!< Zero suppression threshold resolution bit. */
	CVT_V792_BSC2_COMMON_STOP_MSK   = 0x0400,		/*!< enable/disable the common stop acquisition mode bit. */
	CVT_V792_BSC2_AUTO_INC_MSK      = 0x0800,		/*!< enable/disable the automatic increment of the readout pointer bit. */
	CVT_V792_BSC2_EMPTY_EN_MSK      = 0x1000,		/*!< Choose if writing the header and EOB when there are no accepted channels. */
	CVT_V792_BSC2_SLIDE_SUB_EN_MSK  = 0x2000,		/*!< Change operation mode for the sliding scale bit. */
	CVT_V792_BSC2_ALL_TRG_MSK       = 0x4000,		/*!< Choose how to increment the event counter. */
} CVT_V792_BIT_SET_CLEAR_2_MSK;

////////////////////////////////////////////////////////////////////////////////////////////////
// Output buffer macros
////////////////////////////////////////////////////////////////////////////////////////////////
#define CVT_QTP_DATA_TYPE_MSK				0x07000000		/*!< \brief Data type bit masks */

#define CVT_QTP_HEADER						0x02000000		/*!< \brief Header data type */
#define CVT_QTP_EOB							0x04000000		/*!< \brief End of block data type */
#define CVT_QTP_DATUM						0x00000000		/*!< \brief Datum data type */
#define CVT_QTP_NOT_VALID_DATUM				0x06000000		/*!< \brief Not valid datum data type */

#define CVT_QTP_HDR_GEO_MSK					0xf8000000		/*!< \brief The header GEO address msk */
#define CVT_QTP_HDR_CRATE_MSK				0x00ff0000		/*!< \brief The header crate msk */
#define CVT_QTP_HDR_CH_COUNT_MSK			0x00003f00		/*!< \brief The header channel count msk */

#define CVT_QTP_EOB_GEO_MSK					0xf8000000		/*!< \brief The EOB GEO address msk */
#define CVT_QTP_EOB_MAX_EVENT_COUNT			(0x00ffffff+ 1)	/*!< \brief The EOB Max Event count */
#define CVT_QTP_EOB_EVENT_COUNT_MSK			0x00ffffff		/*!< \brief The EOB Event count msk */

//
// V792 specific
#define CVT_V792_DATUM_GEO_MSK				0xf8000000		/*!< \brief The DATUM GEO address msk for V792 */
#define CVT_V792_DATUM_CH_MSK				0x001f0000		/*!< \brief The DATUM channel msk  for V792 ( bit 0x00010000 is always zero for V792N ) */
#define CVT_V792_DATUM_UN_MSK				0x00002000		/*!< \brief The DATUM Under threshold msk  for V792 */
#define CVT_V792_DATUM_OV_MSK				0x00001000		/*!< \brief The DATUM Overflow msk  for V792 */
#define CVT_V792_DATUM_ADC_MSK				0x00000fff		/*!< \brief The DATUM ADC msk  for V792 */
//
// V785 specific
#define CVT_V785_DATUM_GEO_MSK				0xf8000000		/*!< \brief The DATUM GEO address msk for V785 */
#define CVT_V785_DATUM_CH_MSK				0x001f0000		/*!< \brief The DATUM channel msk  for V785 ( bit 0x00010000 is always zero for V785N ) */
#define CVT_V785_DATUM_UN_MSK				0x00002000		/*!< \brief The DATUM Under threshold msk  for V785 */
#define CVT_V785_DATUM_OV_MSK				0x00001000		/*!< \brief The DATUM Overflow msk  for V785 */
#define CVT_V785_DATUM_ADC_MSK				0x00000fff		/*!< \brief The DATUM ADC msk  for V785 */
//
// V775 specific
#define CVT_V775_DATUM_GEO_MSK				0xf8000000		/*!< \brief The DATUM GEO address msk for V775 */
#define CVT_V775_DATUM_CH_MSK				0x001f0000		/*!< \brief The DATUM channel msk  for V775 ( bit 0x00010000 is always zero for V775N ) */
#define CVT_V775_DATUM_VD_MSK				0x00004000		/*!< \brief The DATUM valid bit msk for V775 */
#define CVT_V775_DATUM_UN_MSK				0x00002000		/*!< \brief The DATUM Under threshold msk  for V775 */
#define CVT_V775_DATUM_OV_MSK				0x00001000		/*!< \brief The DATUM Overflow msk  for V775 */
#define CVT_V775_DATUM_ADC_MSK				0x00000fff		/*!< \brief The DATUM ADC msk  for V775 */
//
// V862 specific
#define CVT_V862_DATUM_GEO_MSK				0xf8000000		/*!< \brief The DATUM GEO address msk for V862 */
#define CVT_V862_DATUM_CH_MSK				0x001f0000		/*!< \brief The DATUM channel msk  for V862 */
#define CVT_V862_DATUM_UN_MSK				0x00002000		/*!< \brief The DATUM Under threshold msk  for V862 */
#define CVT_V862_DATUM_OV_MSK				0x00001000		/*!< \brief The DATUM Overflow msk  for V862 */
#define CVT_V862_DATUM_ADC_MSK				0x00000fff		/*!< \brief The DATUM ADC msk  for V862 */
//
// V1785 specific
#define CVT_V1785_DATUM_GEO_MSK				0xf8000000		/*!< \brief The DATUM GEO address msk for V1785 */
#define CVT_V1785_DATUM_CH_MSK				0x001c0000		/*!< \brief The DATUM channel msk  for V1785 */
#define CVT_V1785_DATUM_RG_MSK				0x00020000		/*!< \brief The DATUM range msk  for V1785 */
#define CVT_V1785_DATUM_UN_MSK				0x00002000		/*!< \brief The DATUM Under threshold msk  for V1785 */
#define CVT_V1785_DATUM_OV_MSK				0x00001000		/*!< \brief The DATUM Overflow msk  for V1785 */
#define CVT_V1785_DATUM_ADC_MSK				0x00000fff		/*!< \brief The DATUM ADC msk  for V1785 */
//
// V965 specific
#define CVT_V965_DATUM_GEO_MSK				0xf8000000		/*!< \brief The DATUM GEO address msk for V965 */
#define CVT_V965_DATUM_CH_MSK				0x001e0000		/*!< \brief The DATUM channel msk  for V965 */
#define CVT_V965A_DATUM_CH_MSK				0x001c0000		/*!< \brief The DATUM channel msk  for V965A */
#define CVT_V965_DATUM_RG_MSK				0x00010000		/*!< \brief The DATUM range msk  for V965 */
#define CVT_V965A_DATUM_RG_MSK				0x00020000		/*!< \brief The DATUM range msk  for V965A */
#define CVT_V965_DATUM_UN_MSK				0x00002000		/*!< \brief The DATUM Under threshold msk  for V965 */
#define CVT_V965_DATUM_OV_MSK				0x00001000		/*!< \brief The DATUM Overflow msk  for V965 */
#define CVT_V965_DATUM_ADC_MSK				0x00000fff		/*!< \brief The DATUM ADC msk  for V965 */


#define CVT_QTP_IS_HEADER(data)					((data& CVT_QTP_DATA_TYPE_MSK)== CVT_QTP_HEADER)				/*!< \brief Check if Header data type */
#define CVT_QTP_IS_EOB(data)					((data& CVT_QTP_DATA_TYPE_MSK)== CVT_QTP_EOB)					/*!< \brief Check if End of block data type */
#define CVT_QTP_IS_DATUM(data)					((data& CVT_QTP_DATA_TYPE_MSK)== CVT_QTP_DATUM)					/*!< \brief Check if Datum data type */
#define CVT_QTP_IS_NOT_VALID_DATUM(data)		((data& CVT_QTP_DATA_TYPE_MSK)== CVT_QTP_NOT_VALID_DATUM)		/*!< \brief Check if Not valid datum data type */

#define CVT_QTP_GET_HDR_GEO(data)				((UINT32)((((UINT32)data)& CVT_QTP_HDR_GEO_MSK)>>27))			/*!< \brief Gets the header GEO address */
#define CVT_QTP_GET_HDR_CRATE(data)				((UINT32)((((UINT32)data)& CVT_QTP_HDR_CRATE_MSK)>>16))			/*!< \brief Gets the header crate */
#define CVT_QTP_GET_HDR_CH_COUNT(data)			((UINT32)((((UINT32)data)& CVT_QTP_HDR_CH_COUNT_MSK)>>8))		/*!< \brief Gets the header channel count */

#define CVT_QTP_GET_EOB_GEO(data)				((UINT32)((((UINT32)data)& CVT_QTP_EOB_GEO_MSK)>>27))			/*!< \brief Gets the EOB GEO address */
#define CVT_QTP_GET_EOB_EVENT_COUNT(data)		((UINT32)(((UINT32)data)& CVT_QTP_EOB_EVENT_COUNT_MSK))			/*!< \brief Gets the EOB Event count */

//
// V792 specific
#define CVT_V792_GET_DATUM_GEO(data)			((UINT32)((((UINT32)data)& CVT_V792_DATUM_GEO_MSK)>>27))		/*!< \brief Gets the DATUM GEO address for V792*/
#define CVT_V792_GET_DATUM_UN(data)				((UINT32)((((UINT32)data)& CVT_V792_DATUM_UN_MSK)>>13))			/*!< \brief Gets the DATUM Under threshold for V792*/
#define CVT_V792_GET_DATUM_OV(data)				((UINT32)((((UINT32)data)& CVT_V792_DATUM_OV_MSK)>>12))			/*!< \brief Gets the DATUM Overflow for V792*/
#define CVT_V792_GET_DATUM_ADC(data)			((UINT32)(((UINT32)data)& CVT_V792_DATUM_ADC_MSK))				/*!< \brief Gets the DATUM ADC for V792*/
#define CVT_V792N_GET_DATUM_CH(data)			((UINT32)((((UINT32)data)& CVT_V792_DATUM_CH_MSK)>>17))			/*!< \brief Gets the DATUM channel for V792N */
#define CVT_V792A_GET_DATUM_CH(data)			((UINT32)((((UINT32)data)& CVT_V792_DATUM_CH_MSK)>>16))			/*!< \brief Gets the DATUM channel for V792A */
//
// V785 specific
#define CVT_V785_GET_DATUM_GEO(data)			((UINT32)((((UINT32)data)& CVT_V785_DATUM_GEO_MSK)>>27))		/*!< \brief Gets the DATUM GEO address for V785*/
#define CVT_V785_GET_DATUM_UN(data)				((UINT32)((((UINT32)data)& CVT_V785_DATUM_UN_MSK)>>13))			/*!< \brief Gets the DATUM Under threshold for V785*/
#define CVT_V785_GET_DATUM_OV(data)				((UINT32)((((UINT32)data)& CVT_V785_DATUM_OV_MSK)>>12))			/*!< \brief Gets the DATUM Overflow for V785*/
#define CVT_V785_GET_DATUM_ADC(data)			((UINT32)(((UINT32)data)& CVT_V785_DATUM_ADC_MSK))				/*!< \brief Gets the DATUM ADC for V785*/
#define CVT_V785N_GET_DATUM_CH(data)			((UINT32)((((UINT32)data)& CVT_V785_DATUM_CH_MSK)>>17))			/*!< \brief Gets the DATUM channel for V785N */
#define CVT_V785A_GET_DATUM_CH(data)			((UINT32)((((UINT32)data)& CVT_V785_DATUM_CH_MSK)>>16))			/*!< \brief Gets the DATUM channel for V785A */
//
// V775 specific
#define CVT_V775_GET_DATUM_GEO(data)			((UINT32)((((UINT32)data)& CVT_V775_DATUM_GEO_MSK)>>27))		/*!< \brief Gets the DATUM GEO address for V775*/
#define CVT_V775_GET_DATUM_VD(data)				((UINT32)((((UINT32)data)& CVT_V775_DATUM_VD_MSK)>>14))			/*!< \brief Gets the DATUM Valid bit for V775*/
#define CVT_V775_GET_DATUM_UN(data)				((UINT32)((((UINT32)data)& CVT_V775_DATUM_UN_MSK)>>13))			/*!< \brief Gets the DATUM Under threshold for V775*/
#define CVT_V775_GET_DATUM_OV(data)				((UINT32)((((UINT32)data)& CVT_V775_DATUM_OV_MSK)>>12))			/*!< \brief Gets the DATUM Overflow for V775*/
#define CVT_V775_GET_DATUM_ADC(data)			((UINT32)(((UINT32)data)& CVT_V775_DATUM_ADC_MSK))				/*!< \brief Gets the DATUM ADC for V775*/
#define CVT_V775N_GET_DATUM_CH(data)			((UINT32)((((UINT32)data)& CVT_V775_DATUM_CH_MSK)>>17))			/*!< \brief Gets the DATUM channel for V775N */
#define CVT_V775A_GET_DATUM_CH(data)			((UINT32)((((UINT32)data)& CVT_V775_DATUM_CH_MSK)>>16))			/*!< \brief Gets the DATUM channel for V775A */

//
// V862 specific
#define CVT_V862_GET_DATUM_GEO(data)			((UINT32)((((UINT32)data)& CVT_V862_DATUM_GEO_MSK)>>27))		/*!< \brief Gets the DATUM GEO address for V862*/
#define CVT_V862_GET_DATUM_UN(data)				((UINT32)((((UINT32)data)& CVT_V862_DATUM_UN_MSK)>>13))			/*!< \brief Gets the DATUM Under threshold for V862*/
#define CVT_V862_GET_DATUM_OV(data)				((UINT32)((((UINT32)data)& CVT_V862_DATUM_OV_MSK)>>12))			/*!< \brief Gets the DATUM Overflow for V862*/
#define CVT_V862_GET_DATUM_ADC(data)			((UINT32)(((UINT32)data)& CVT_V862_DATUM_ADC_MSK))				/*!< \brief Gets the DATUM ADC for V862*/
#define CVT_V862_GET_DATUM_CH(data)				((UINT32)((((UINT32)data)& CVT_V862_DATUM_CH_MSK)>>16))			/*!< \brief Gets the DATUM channel for V862 */
//
// V1785 specific
#define CVT_V1785_GET_DATUM_GEO(data)			((UINT32)((((UINT32)data)& CVT_V1785_DATUM_GEO_MSK)>>27))		/*!< \brief Gets the DATUM GEO address for V1785*/
#define CVT_V1785_GET_DATUM_RG(data)			((UINT32)((((UINT32)data)& CVT_V1785_DATUM_RG_MSK)>>17))		/*!< \brief Gets the DATUM range for V1785*/
#define CVT_V1785_GET_DATUM_UN(data)			((UINT32)((((UINT32)data)& CVT_V1785_DATUM_UN_MSK)>>13))		/*!< \brief Gets the DATUM Under threshold for V1785*/
#define CVT_V1785_GET_DATUM_OV(data)			((UINT32)((((UINT32)data)& CVT_V1785_DATUM_OV_MSK)>>12))		/*!< \brief Gets the DATUM Overflow for V1785*/
#define CVT_V1785_GET_DATUM_ADC(data)			((UINT32)(((UINT32)data)& CVT_V1785_DATUM_ADC_MSK))				/*!< \brief Gets the DATUM ADC for V1785*/
#define CVT_V1785_GET_DATUM_CH(data)			((UINT32)((((UINT32)data)& CVT_V1785_DATUM_CH_MSK)>>18))		/*!< \brief Gets the DATUM channel for V1785 */
//
// V965 specific
#define CVT_V965_GET_DATUM_GEO(data)			((UINT32)((((UINT32)data)& CVT_V965_DATUM_GEO_MSK)>>27))		/*!< \brief Gets the DATUM GEO address for V965*/
#define CVT_V965_GET_DATUM_RG(data)				((UINT32)((((UINT32)data)& CVT_V965_DATUM_RG_MSK)>>16))			/*!< \brief Gets the DATUM range for V965 */
#define CVT_V965A_GET_DATUM_RG(data)			((UINT32)((((UINT32)data)& CVT_V965A_DATUM_RG_MSK)>>17))		/*!< \brief Gets the DATUM range for V965A */
#define CVT_V965_GET_DATUM_UN(data)				((UINT32)((((UINT32)data)& CVT_V965_DATUM_UN_MSK)>>13))			/*!< \brief Gets the DATUM Under threshold for V965 */
#define CVT_V965_GET_DATUM_OV(data)				((UINT32)((((UINT32)data)& CVT_V965_DATUM_OV_MSK)>>12))			/*!< \brief Gets the DATUM Overflow for V965*/
#define CVT_V965_GET_DATUM_ADC(data)			((UINT32)(((UINT32)data)& CVT_V965_DATUM_ADC_MSK))				/*!< \brief Gets the DATUM ADC for V965*/
#define CVT_V965_GET_DATUM_CH(data)				((UINT32)((((UINT32)data)& CVT_V965_DATUM_CH_MSK)>>17))			/*!< \brief Gets the DATUM channel for V965 */
#define CVT_V965A_GET_DATUM_CH(data)			((UINT32)((((UINT32)data)& CVT_V965A_DATUM_CH_MSK)>>18))		/*!< \brief Gets the DATUM channel for V965A */


////////////////////////////////////////////////////////////////////////////////////////////////
// Global variables declaration
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
// Global methods declaration
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
//
//     B O A R D S   H A N D L I N G
//
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_open( cvt_V792_data* p_data, UINT16 base_address, long vme_handle, CVT_QTP_TYPES type);
*   \brief   V792 VME boards data initialization
*            
*            Provides specific handling for V792 boards opening.
*   \param   p_data Pointer to board data
*   \param   base_address The board base address (MSW)
*   \param   vme_handle The VME handle
*   \param   type The board type : it must be a \ref CVT_QTP_TYPES value
*   \return  TRUE: board successfully opened
*   \note    Must be called before any other board specific API.
*   \sa      CVT_QTP_TYPES 
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_open( cvt_V792_data* p_data, UINT16 base_address, long vme_handle, CVT_QTP_TYPES type);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_close( cvt_V792_data* p_data);
*   \brief   V792 VME boards closing and resource free
*            
*            Provides specific handling for V792 boards closing.
*   \param   p_data Pointer to board data
*   \return  TRUE: board successfully closed
*   \note    Must be called when done with any other board specific API.
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_close( cvt_V792_data* p_data);

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
/*! \fn      BOOL cvt_V792_set_geo_address( cvt_V792_data* p_data, UINT8 value);
*   \brief   Sets the geographical address.
*            
*            Writes the geographical address register.
*   \param   p_data Pointer to board data
*   \param   value The value to write
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_geo_address( cvt_V792_data* p_data, UINT8 value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_get_geo_address( cvt_V792_data* p_data, UINT8 *p_value);
*   \brief   Gets the geographical address.
*            
*            Reads the geographical address register.
*   \param   p_data Pointer to board data
*   \param   p_value The value read
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_get_geo_address( cvt_V792_data* p_data, UINT8 *p_value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_MCST_CBLT_address( cvt_V792_data* p_data, UINT8 value);
*   \brief   Sets the MCST\\CBLT address.
*            
*            Writes the geographical address register. 
*            This register contains the most significant bits of the MCST/CBLT address of the module set via VME.
*   \param   p_data Pointer to board data
*   \param   value The value to write
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_MCST_CBLT_address( cvt_V792_data* p_data, UINT8 value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_get_MCST_CBLT_address( cvt_V792_data* p_data, UINT8 *p_value);
*   \brief   Gets the MCST\\CBLT address.
*            
*            Reads the MCST\\CBLT address register.
*            This register contains the most significant bits of the MCST/CBLT address of the module set via VME.
*   \param   p_data Pointer to board data
*   \param   p_value The value read
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_get_MCST_CBLT_address( cvt_V792_data* p_data, UINT8 *p_value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_bit_set_1( cvt_V792_data* p_data, UINT16 value);
*   \brief   Set the bit set 1 register.
*            
*            Writes the specified value to the bit set 1 register.
*   \param   p_data Pointer to board data
*   \param   value The value to set
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_bit_set_1( cvt_V792_data* p_data, UINT16 value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_bit_clear_1( cvt_V792_data* p_data, UINT16 value);
*   \brief   Set the bit clear 1 register.
*            
*            Writes the specified value to the bit clear 1 register.
*   \param   p_data Pointer to board data
*   \param   value The value to set
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_bit_clear_1( cvt_V792_data* p_data, UINT16 value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_interrupt_level( cvt_V792_data* p_data, UINT8 value);
*   \brief   Set the interrupt level register.
*            
*            Writes the specified value to the interrupt level register.
*   \param   p_data Pointer to board data
*   \param   value The value to set
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_interrupt_level( cvt_V792_data* p_data, UINT8 value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_get_interrupt_level( cvt_V792_data* p_data, UINT8 *p_value);
*   \brief   Set the interrupt level register.
*            
*            Reads the interrupt level register.
*   \param   p_data Pointer to board data
*   \param   p_value The value read
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_get_interrupt_level( cvt_V792_data* p_data, UINT8 *p_value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_interrupt_vector( cvt_V792_data* p_data, UINT8 value);
*   \brief   Set the interrupt level register.
*            
*            Writes the specified value to the interrupt vector register.
*   \param   p_data Pointer to board data
*   \param   value The value to set
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_interrupt_vector( cvt_V792_data* p_data, UINT8 value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_get_interrupt_vector( cvt_V792_data* p_data, UINT8 *p_value);
*   \brief   Get the interrupt vector register.
*            
*            Reads the interrupt vector register.
*   \param   p_data Pointer to board data
*   \param   p_value The value read
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_get_interrupt_vector( cvt_V792_data* p_data, UINT8 *p_value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_get_status_1( cvt_V792_data* p_data, UINT16 *p_value);
*   \brief   Get the status 1 register.
*            
*            Reads the status 1 register.
*   \param   p_data Pointer to board data
*   \param   p_value The value read
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_get_status_1( cvt_V792_data* p_data, UINT16 *p_value);
	
////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_control_1( cvt_V792_data* p_data, UINT16 value);
*   \brief   Set to control 1 register.
*            
*            Writes the specified value to the control 1 register.
*   \param   p_data Pointer to board data
*   \param   value The value to write
*   \return  TRUE: Procedure successfully executed
*   \sa      CVT_V792_CONTROL_1_MSK
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_control_1( cvt_V792_data* p_data, UINT16 value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_get_control_1( cvt_V792_data* p_data, UINT16 *p_value);
*   \brief   Get the control 1 register.
*            
*            Reads the control 1 register.
*   \param   p_data Pointer to board data
*   \param   p_value The value read
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_get_control_1( cvt_V792_data* p_data, UINT16 *p_value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_clear_bitmask_control_1( cvt_V792_data* p_data, CVT_V792_CONTROL_1_MSK value);
*   \brief   Clear a bitmask to control 1 register.
*            
*            Clears the bits of the bitmask into the register. The inverted bitmask is ANDed to the actual regiter content
*   \param   p_data Pointer to board data
*   \param   value The bitmask to clear: may be a combination of \ref CVT_V792_CONTROL_1_MSK bitmask
*   \return  TRUE: Procedure successfully executed
*   \sa      CVT_V792_CONTROL_1_MSK
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_clear_bitmask_control_1( cvt_V792_data* p_data, CVT_V792_CONTROL_1_MSK value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_bitmask_control_1( cvt_V792_data* p_data, CVT_V792_CONTROL_1_MSK value);
*   \brief   Set a bitmask to control 1 register.
*            
*            Sets the bits of the bitmask into the register. The bitmask is ORed to the actual regiter content
*   \param   p_data Pointer to board data
*   \param   value The bitmask to set : may be a combination of \ref CVT_V792_CONTROL_1_MSK bitmask
*   \return  TRUE: Procedure successfully executed
*   \sa      CVT_V792_CONTROL_1_MSK
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_bitmask_control_1( cvt_V792_data* p_data, CVT_V792_CONTROL_1_MSK value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_clear_bitmask_control_1( cvt_V792_data* p_data, CVT_V792_CONTROL_1_MSK value);
*   \brief   Clear a bitmask to control 1 register.
*            
*            Clears the bits of the bitmask into the register. The inverted bitmask is ANDed to the actual regiter content
*   \param   p_data Pointer to board data
*   \param   value The bitmask to clear: may be a combination of \ref CVT_V792_CONTROL_1_MSK bitmask
*   \return  TRUE: Procedure successfully executed
*   \sa      CVT_V792_CONTROL_1_MSK
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_clear_bitmask_control_1( cvt_V792_data* p_data, CVT_V792_CONTROL_1_MSK value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_ader_high( cvt_V792_data* p_data, UINT8 value);
*   \brief   Set the address decoder high register.
*            
*            Writes the specified value to the address decoder high register.
*   \param   p_data Pointer to board data
*   \param   value The value to set
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_ader_high( cvt_V792_data* p_data, UINT8 value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_get_ader_high( cvt_V792_data* p_data, UINT8 *p_value);
*   \brief   Get the address decoder high register.
*            
*            Reads the address decoder high register.
*   \param   p_data Pointer to board data
*   \param   p_value The value read
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_get_ader_high( cvt_V792_data* p_data, UINT8 *p_value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_ader_low( cvt_V792_data* p_data, UINT8 value);
*   \brief   Set the address decoder low register.
*            
*            Writes the specified value to the address decoder low register.
*   \param   p_data Pointer to board data
*   \param   value The value to set
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_ader_high( cvt_V792_data* p_data, UINT8 value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_get_ader_low( cvt_V792_data* p_data, UINT8 *p_value);
*   \brief   Get the address decoder high register.
*            
*            Reads the address decoder low register.
*   \param   p_data Pointer to board data
*   \param   p_value The value read
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_get_ader_low( cvt_V792_data* p_data, UINT8 *p_value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_single_shot_reset( cvt_V792_data* p_data);
*   \brief   Perfoms a single shot reset.
*            
*            Writes a dummy value to the single shot reset register.
*   \param   p_data Pointer to board data
*   \return  TRUE: Procedure successfully executed
*   \warning This register must be used very carefully and for debugging purposes only. 
*            In order to reset the board, it is recommended to use the Bit Set 1 Register
*   \sa      cvt_V792_set_bitmask_control_1
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_single_shot_reset( cvt_V792_data* p_data);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_MCST_CBLT_control( cvt_V792_data* p_data, UINT8 value);
*   \brief   Set the MCST/CBLT Control register.
*            
*            Writes the specified value to the MCST/CBLT Control register.
*   \param   p_data Pointer to board data
*   \param   value The value to set
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_MCST_CBLT_control( cvt_V792_data* p_data, UINT8 value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_get_MCST_CBLT_control( cvt_V792_data* p_data, UINT8 *p_value);
*   \brief   Get the MCST/CBLT Control register.
*            
*            Reads the MCST/CBLT Control register.
*   \param   p_data Pointer to board data
*   \param   p_value The value read
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_get_MCST_CBLT_control( cvt_V792_data* p_data, UINT8 *p_value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_event_trigger( cvt_V792_data* p_data, UINT8 value);
*   \brief   Set the Event trigger register.
*            
*            Writes the specified value to the Event trigger register.
*   \param   p_data Pointer to board data
*   \param   value The value to set
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_event_trigger( cvt_V792_data* p_data, UINT8 value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_get_event_trigger( cvt_V792_data* p_data, UINT8 *p_value);
*   \brief   Get the Event trigger register.
*            
*            Reads the event trigger register.
*   \param   p_data Pointer to board data
*   \param   p_value The value read
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_get_event_trigger( cvt_V792_data* p_data, UINT8 *p_value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_get_status_2( cvt_V792_data* p_data, UINT16 *p_value);
*   \brief   Get the status 1 register.
*            
*            Reads the status 2 register.
*   \param   p_data Pointer to board data
*   \param   p_value The value read
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_get_status_2( cvt_V792_data* p_data, UINT16 *p_value);
	
////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_threshold( cvt_V792_data* p_data, cvt_V792_threshold_id threshold_id, UINT16 value);
*   \brief   Set the specified threshold register.
*            
*            Writes the specified value to the specified threshold register. If threshold_id is CVT_V792_THRESHOLD_ALL
*            every threshold will be setted to the same value.
*   \param   p_data Pointer to board data
*   \param   threshold_id The threshold to set. It must be a \ref cvt_V792_threshold_id valid id.
*   \param   value The value to set : Threshold value and Kill bit (\ref CVT_V792_THRESHOLD_KILL_MSK)
*                  If the kill bit be isn't setted the channel will be disabled
*   \return  TRUE: Procedure successfully executed
*   \sa      cvt_V792_threshold_id 
*   \sa      CVT_V792_THRESHOLD_KILL_MSK
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_threshold( cvt_V792_data* p_data, cvt_V792_threshold_id threshold_id, UINT16 value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_dual_range_threshold( cvt_V792_data* p_data, cvt_V792_threshold_id threshold_id, UINT16 value_high, UINT16 value_low);
*   \brief   Set the specified threshold register for dual range board (V965 and V1785).
*            
*            Writes the specified value to the specified threshold register as high and low part. 
*            If threshold_id is CVT_V792_THRESHOLD_ALL every threshold will be setted to the same value.
*   \param   p_data Pointer to board data
*   \param   threshold_id The threshold to set. It must be a \ref cvt_V792_threshold_id valid id.
*   \param   value_high The high value to set : Threshold value and Kill bit (\ref CVT_V792_THRESHOLD_KILL_MSK)
*                  If the kill bit be isn't setted the channel will be disabled
*   \param   value_low The low value to set : Threshold value and Kill bit (\ref CVT_V792_THRESHOLD_KILL_MSK)
*                  If the kill bit be isn't setted the channel will be disabled
*   \return  TRUE: Procedure successfully executed
*   \sa      cvt_V792_threshold_id 
*   \sa      CVT_V792_THRESHOLD_KILL_MSK
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_dual_range_threshold( cvt_V792_data* p_data, cvt_V792_threshold_id threshold_id, UINT16 value_high, UINT16 value_low);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_bit_set_2( cvt_V792_data* p_data, UINT16 value);
*   \brief   Set the bit set 2 register.
*            
*            Writes the specified value to the bit set 2 register.
*   \param   p_data Pointer to board data
*   \param   value The value to set
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_bit_set_2( cvt_V792_data* p_data, UINT16 value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_bit_clear_2( cvt_V792_data* p_data, UINT16 value);
*   \brief   Set the bit clear 2 register.
*            
*            Writes the specified value to the bit clear 2 register.
*   \param   p_data Pointer to board data
*   \param   value The value to set
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_bit_clear_2( cvt_V792_data* p_data, UINT16 value);
	
////////////////////////////////////////////////////////////////////////////////////////////////
//
//     L E V E L   2   A P I s
//
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_sliding_scale( cvt_V792_data* p_data, BOOL enable);
*   \brief   Setups the sliding scale feature.
*            
*            Enable/ Disable sliding scale.
*   \param   p_data Pointer to board data
*   \param   enable Enable/disable the sliding scale
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_sliding_scale( cvt_V792_data* p_data, BOOL enable);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_zero_suppression( cvt_V792_data* p_data, BOOL enable, BOOL step_threshold, const UINT16* thresholds_buff);
*   \brief   Setups the zero suppression feature.
*            
*            Enable/ Disable zero suppression and setups other relevant parameters.
*   \param   p_data Pointer to board data
*   \param   enable Enable/disable the zero suppression (low threshold check).
*   \param   step_threshold If setted the comparison is ADC CONVERTED VALUE < THRESHOLD VALUE x 2, otherwise ADC CONVERTED VALUE < THRESHOLD VALUE x 16
*   \param   thresholds_buff The thresholds buffer: It must contain as many values as specified for the board type ( \ref CVT_QTP_NUM_CHANNELS).
*                            It may be NULL (e.g. when disabling this feature or the thresholds are already setted)
*                            The value consists of the threshold level and the kill bit (\ref CVT_V792_THRESHOLD_KILL_MSK)
*                            If the kill bit be isn't setted the channel will be disabled
*   \return  TRUE: Procedure successfully executed
*   \sa      CVT_V792_THRESHOLD_KILL_MSK
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_zero_suppression( cvt_V792_data* p_data, BOOL enable, BOOL step_threshold, const UINT16* thresholds_buff);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_dual_range_zero_suppression( cvt_V792_data* p_data, BOOL enable, BOOL step_threshold, const UINT16* thresholds_buff_high, const UINT16* thresholds_buff_low)
*   \brief   Setups the zero suppression feature for dual range boards ( V965 or V1785)
*            
*            Enable/ Disable zero suppression and setups other relevant parameters.
*   \param   p_data Pointer to board data
*   \param   enable Enable/disable the zero suppression (low threshold check).
*   \param   step_threshold If setted the comparison is ADC CONVERTED VALUE < THRESHOLD VALUE x 2, otherwise ADC CONVERTED VALUE < THRESHOLD VALUE x 16
*   \param   thresholds_buff_high The high thresholds buffer: It must contain as many values as specified for the board type ( \ref CVT_QTP_NUM_CHANNELS).
*                            It may be NULL (e.g. when disabling this feature or the thresholds are already setted)
*                            The value consists of the threshold level and the kill bit (\ref CVT_V792_THRESHOLD_KILL_MSK)
*                            If the kill bit be isn't setted the channel will be disabled
*   \param   thresholds_buff_low The low thresholds buffer: It must contain as many values as specified for the board type ( \ref CVT_QTP_NUM_CHANNELS).
*                            It may be NULL (e.g. when disabling this feature or the thresholds are already setted)
*                            The value consists of the threshold level and the kill bit (\ref CVT_V792_THRESHOLD_KILL_MSK)
*                            If the kill bit be isn't setted the channel will be disabled
*   \return  TRUE: Procedure successfully executed
*   \sa      CVT_V792_THRESHOLD_KILL_MSK
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_dual_range_zero_suppression( cvt_V792_data* p_data, BOOL enable, BOOL step_threshold, const UINT16* thresholds_buff_high, const UINT16* thresholds_buff_low);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_overflow_suppression( cvt_V792_data* p_data, BOOL enable);
*   \brief   Setups the overflow suppression feature.
*            
*            Enable/ Disable overflow suppression and setups other relevant parameters.
*   \param   p_data Pointer to board data
*   \param   enable Enable/disable the overflow suppression.
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_overflow_suppression( cvt_V792_data* p_data, BOOL enable);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_MEB( cvt_V792_data* p_data, BOOL empty_enable, BOOL auto_incr_enable);
*   \brief   Setups the Multiple Event Buffer.
*            
*            Setups the relevant parameters for the Multiple Event Buffer.
*   \param   p_data Pointer to board data
*   \param   empty_enable Enable/disable storing events with no data (only header and EOB).
*   \param   auto_incr_enable Enable/disable the event counter autoincrement feature.
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_MEB( cvt_V792_data* p_data, BOOL empty_enable, BOOL auto_incr_enable);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_event_counter( cvt_V792_data* p_data, BOOL count_all_events, BOOL reset_count);
*   \brief   Setups the event counter.
*            
*            Setups the relevant parameters for the event counter.
*   \param   p_data Pointer to board data
*   \param   count_all_events Enable counting all events or accepted ones only.
*   \param   reset_count If TRUE resets the event counter.
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_event_counter( cvt_V792_data* p_data, BOOL count_all_events, BOOL reset_count);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_read_MEB( cvt_V792_data* p_data, void* p_buff, UINT32* p_buff_size);
*   \brief   Reads data from Multiple event buffer and stores to user buffer
*            
*            Call cvt_FIFO_BLT_read for Multiple Event Buffer and stores data into user buffer.
*            If the returned data is just a V792 not valid datum, this is discarded.
*   \param   p_data Pointer to board data
*   \param   p_buff The target buffer: caller allocated
*   \param   p_buff_size The buffer size (bytes). On exit will holds the number of bytes really read.
*   \return  TRUE: Procedure successfully executed
*   \sa      cvt_FIFO_BLT_read
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_read_MEB( cvt_V792_data* p_data, void* p_buff, UINT32* p_buff_size);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_fast_clear_window( cvt_V792_data* p_data, UINT16 value);
*   \brief   Set the fast clear window register.
*            
*            Writes the specified value to the fast clear window register.
*   \param   p_data Pointer to board data
*   \param   value The value to set
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_fast_clear_window( cvt_V792_data* p_data, UINT16 value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_acquisition_mode( cvt_V792_data* p_data, BOOL sliding_scale_enable, BOOL zero_suppression_enable, BOOL overflow_suppression_enable, BOOL valid_suppression_enable, BOOL common_stop_enable, BOOL empty_enable, BOOL count_all_events);
*   \brief   Setups the acquisition mode parameters.
*            
*            Setups the relevant parameters for the acquisition mode.
*   \param   p_data Pointer to board data
*   \param   sliding_scale_enable Enable/Disable the sliding scale feature.
*   \param   zero_suppression_enable Enable/Disable the zero suppression feature.
*   \param   overflow_suppression_enable Enable/Disable the overflow suppression feature.
*   \param   valid_suppression_enable Enable/Disable the valid suppression feature (V775 ONLY).
*   \param   common_stop_enable Enable/Disable the common stop acquisition mode (V775 ONLY).
*   \param   empty_enable Enable/Disable empty event storing feature.
*   \param   count_all_events Enable counting all triggers or accepted one only.
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_acquisition_mode( cvt_V792_data* p_data, BOOL sliding_scale_enable, BOOL zero_suppression_enable, BOOL overflow_suppression_enable, BOOL valid_suppression_enable, BOOL common_stop_enable, BOOL empty_enable, BOOL count_all_events);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_interrupt( cvt_V792_data* p_data, UINT8 level, UINT8 vector, UINT8 event_number);
*   \brief   Setups interrupt parameters.
*            
*            Setups the relevant parameters for interrupt usage.
*   \param   p_data Pointer to board data
*   \param   level The interrupt level.
*   \param   vector The interrupt vector.
*   \param   event_number The number of events to get an interrupt (setting 0 interrupt feature is disabled).
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_interrupt( cvt_V792_data* p_data, UINT8 level, UINT8 vector, UINT8 event_number);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_readout_mode( cvt_V792_data* p_data, BOOL bus_error_enable, BOOL block_end_enable, BOOL align64_enable);
*   \brief   Setups data readout mode parameters.
*            
*            Setups the relevant parameters for data readout.
*   \param   p_data Pointer to board data
*   \param   bus_error_enable Enable bus error: the module is enabled to generate a Bus error to finish a block transfer. 
*   \param   block_end_enable Enable block end: the module sends all data to the CPU until the first EOB word
*                             (end of first event) is reached;
*   \param   align64_enable Enable align 64: the module is enabled to add dummy word when the number of words is odd.
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_readout_mode( cvt_V792_data* p_data, BOOL bus_error_enable, BOOL block_end_enable, BOOL align64_enable);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_get_status( cvt_V792_data* p_data, BOOL *p_is_data_ready, BOOL *p_is_busy, BOOL *p_is_term_on, BOOL *p_is_buffer_full);
*   \brief   Gets information about board status.
*            
*   \param   p_data Pointer to board data
*   \param   p_is_data_ready Returns TRUE if data is ready. 
*   \param   p_is_busy Returns TRUE if board is busy. 
*   \param   p_is_term_on Returns TRUE if term bit is on. 
*   \param   p_is_buffer_full Returns TRUE if buffer is full. 
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_get_status( cvt_V792_data* p_data, BOOL *p_is_data_ready, BOOL *p_is_busy, BOOL *p_is_term_on, BOOL *p_is_buffer_full);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_crate_number( cvt_V792_data* p_data, UINT8 crate_number);
*   \brief   Setups the crate number.
*            
*   \param   p_data Pointer to board data
*   \param   crate_number The crate number. 
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_crate_number( cvt_V792_data* p_data, UINT8 crate_number);
	
////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_pedestal( cvt_V792_data* p_data, UINT8 value);
*   \brief   Set the pedestal value.
*            
*            Writes the specified value to the Iped register
*   \param   p_data Pointer to board data
*   \param   value The value to set
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_pedestal( cvt_V792_data* p_data, UINT8 value);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_get_event_counter( cvt_V792_data* p_data, UINT32* p_counter);
*   \brief   Gets the event counter value.
*            
*            Gets the event counter actual value.
*   \param   p_data Pointer to board data
*   \param   p_counter Pointer to the counter value.
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_get_event_counter( cvt_V792_data* p_data, UINT32* p_counter);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_get_system_info( cvt_V792_data* p_data, UINT16 *p_firmware_rev, UINT8 *p_piggy_back_type, UINT16 *p_serial_number);
*   \brief   Gets board's system information.
*            
*            Reads the firmware revision register, the board type and the serial number.
*   \param   p_data Pointer to board data
*   \param   p_firmware_rev The firmare release (MSB major release, LSB minor release).
*   \param   p_piggy_back_type The mounted piggy back type.
*   \param   p_serial_number The serial number.
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_get_system_info( cvt_V792_data* p_data, UINT16 *p_firmware_rev, UINT8 *p_piggy_back_type, UINT16 *p_serial_number);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      const char * cvt_V792_get_piggy_back_str( cvt_V792_data* p_data, UINT8 piggy_back_type);
*   \brief   Returns ths mounted piggy back type string.
*            
*   \param   p_data Pointer to board data
*   \param   piggy_back_type The piggy back type numerical value (as returned from \ref BOOL cvt_V792_get_system_info( cvt_V792_data* p_data, UINT16 *p_firmware_rev, UINT8 *p_piggy_back_type, UINT16 *p_serial_number).
*   \return  A user friendly string describing the piggy back type 
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API const char * cvt_V792_get_piggy_back_str( cvt_V792_data* p_data, UINT8 piggy_back_type);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_software_reset( cvt_V792_data* p_data);
*   \brief   Performs a software reset.
*            
*            Writes a dummy value into CVT_V792_SINGLE_SHOT_RESET register.
*   \param   p_data Pointer to board data
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_software_reset( cvt_V792_data* p_data);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_data_clear( cvt_V792_data* p_data);
*   \brief   Performs a data clear.
*            
*            Pulses the \ref CVT_V792_BSC2_CLEAR_DATA_MSK bit and clears the event counter.
*   \param   p_data Pointer to board data
*   \return  TRUE: Procedure successfully executed
*   \sa      CVT_V792_BSC2_CLEAR_DATA_MSK
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_data_clear( cvt_V792_data* p_data);

#ifdef CVT_V792_USE_DATA_QUEUE

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_read_data( cvt_V792_data *p_data);
*   \brief   Reads new data from MEB and stores into internal queue.
*            
*            Reads the largest allowed buffer (queue free size) from MEB and stores into queue
*   \param   p_data Pointer to board data
*   \return  TRUE: Procedure successfully executed
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_read_data( cvt_V792_data *p_data);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_peek_event( cvt_V792_data *p_data, UINT32 *out_buff, long *p_out_buff_size, UINT32 *p_event_count, BOOL dequeue);
*   \brief   Tries to peek (i.e. copy but not removing) an event from queue and save data to output buffer
*            
*            Scans the queue beginning from m_queue_ini, searching for an event block, and stores data output stream. 
*   \param   p_data Pointer to board data
*   \param   out_buff Pointer to output buffer
*   \param   p_out_buff_size  Size of out_buff (UINT32). The number of data written into out_buff (UINT32) on exit.
*   \param   p_event_count The found event count
*   \param   dequeue If TRUE removes the found events from queue
*   \return  TRUE: some event found FALSE elsewhere
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_peek_event( cvt_V792_data *p_data, UINT32 *out_buff, long *p_out_buff_size, UINT32 *p_event_count, BOOL dequeue);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      cvt_V792_peek_channel_data( cvt_V792_data *p_data, UINT8 channel_index, UINT16 *out_buff, long *p_out_buff_size, UINT32 *p_event_count, BOOL dequeue)
*   \brief   Tries to peek (i.e. copy but not removing) an event from queue and save ADC data for the specific channel to output buffer
*            
*            Scans the queue beginning from m_queue_ini, searching for an event block, and stores data output stream. 
*   \param   p_data Pointer to board data
*   \param   channel_index The channel to get data from
*   \param   out_buff Pointer to output buffer
*   \param   p_out_buff_size  Size of out_buff (UINT32). The number of data written into out_buff (UINT32) on exit.
*   \param   p_event_count Requested event count. The found event count, on exit
*   \param   dequeue If TRUE removes the found events from queue. !!! WARNING !!! the whole event will be deleted (i.e. other channels' data will be lost too)
*   \return  TRUE: some event found, FALSE elsewhere
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_peek_channel_data( cvt_V792_data *p_data, UINT8 channel_index, UINT16 *out_buff, long *p_out_buff_size, UINT32 *p_event_count, BOOL dequeue);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_inqueue( cvt_V792_data* p_data, const UINT32* in_buff, UINT32 in_buff_size);
*   \brief   Tries to put in_buff_size data word to queue
*            
*            Checks for the queue available free words and stores data beginning from m_queue_end. 
*            m_queue_end will be update
*   \param   p_data Pointer to board data
*   \param   in_buff Pointer to input buffer
*   \param   in_buff_size  Size of in_buff (UINT32)
*   \return  TRUE: data stored correctly FALSE elsewhere (e.g. queue overflow)
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_inqueue( cvt_V792_data* p_data, const UINT32* in_buff, UINT32 in_buff_size);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_dequeue( cvt_V792_data* p_data, UINT32 *out_buff, UINT32 out_buff_size);
*   \brief   Tries to get out_buff_size data word from the queue
*            
*            Checks for queue available data words and gets data beginning from m_queue_ini. 
*            m_queue_ini will be update: if out_buff!= NULL the data will be stored into it.
*   \param   p_data Pointer to board data
*   \param   out_buff Pointer to output buffer (may be NULL)
*   \param   out_buff_size  Size of out_buff (UINT32)
*   \return  TRUE: data removed correctly, FALSE elsewhere (e.g. queue underflow)
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_dequeue( cvt_V792_data* p_data, UINT32 *out_buff, UINT32 out_buff_size);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      long cvt_V792_get_queue_free( cvt_V792_data* p_data);
*   \brief   Gets the number of free queue space (words)
*            
*   \param   p_data Pointer to board data
*   \return  The number of free queue space (UINT32 words)
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API long cvt_V792_get_queue_free( cvt_V792_data* p_data);

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      long cvt_V792_get_queue_length( cvt_V792_data* p_data);
*   \brief   Gets the queue length
*            
*   \param   p_data Pointer to board data
*   \return  The queue length (UINT32 words)
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API long cvt_V792_get_queue_length( cvt_V792_data* p_data);

#endif  // CVT_V792_USE_DATA_QUEUE

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      BOOL cvt_V792_set_MCST_CBLT( cvt_V792_data* p_data, UINT8 address, MCST_CBLT_board_pos pos);
*   \brief   Setups MCST/CBLT parameters this board.
*            
*            Setups the relevant parameters for MCST/CBLT usage.
*   \param   p_data Pointer to board data.
*   \param   address The MCST/CBLT address.
*   \param   pos The board position into the MCST / CBLT chain: it must be a \ref MCST_CBLT_board_pos identifier
*   \return  TRUE: Procedure successfully executed
*   \sa      MCST_CBLT_board_pos
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API BOOL cvt_V792_set_MCST_CBLT( cvt_V792_data* p_data, UINT8 address, MCST_CBLT_board_pos pos);

////////////////////////////////////////////////////////////////////////////////////////////////
//
//     M I S C E L L A N E O U S   A P I s
//
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
/*! \fn      const char* cvt_V792_SW_rev( void);
*   \brief   Returns the SDK software release
*            
*   \return  const char*: The SDK revision string
*/
////////////////////////////////////////////////////////////////////////////////////////////////
CVT_DLL_API const char* cvt_V792_SW_rev( void);

#endif  // __CVT_V792_DEF_H

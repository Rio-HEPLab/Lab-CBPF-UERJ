########################################################################
#                                                                      
#              --- CAEN SpA - Computing Division ---                   
#                                                                      
#   V1190Demo Software Project                                        
#                                                                      
#   Created  :  Apr 2007      (Rel. 1.0)                                             
#                                                                      
#   Auth: NDA
#                                                                      
########################################################################

EXE	=	CAENqtpDemo

CC	=	gcc

#COPTS	=	-fPIC -DLINUX -Wall 
COPTS	=	-g -fPIC -DLINUX -Wall 

FLAGS	=	-Wall -s
#FLAGS	=	-Wall

DEPLIBS	=       -lCAENVME -lCAENqtpSDK

LIBS	=	

#INCLUDEDIR =	-I. -I../include
INCLUDEDIR =	-I. -I../include -I/opt/CAEN/CAENqtp-1.2.1/SDK -I/opt/CAEN/CAENqtp-1.2.1/SDK/common/include/CAENVMELib

OBJS		= main.o user_settings.o keyb.o
OBJS_CAENQDC	= CAENQDC.o user_settings.o keyb.o

#INCLUDES =	../include/CAENVMELib/CAENVMElib.h ../include/CAENVMELib/CAENVMEoslib.h ../include/CAENVMELib/CAENVMEtypes.h ../include/CAENVMETool/cvt_board_commons.h ../include/CAENVMETool/cvt_common_defs.h ../include/CAENVMETool/cvt_V792.h keyb.h user_settings.h

#########################################################################

all	:	$(EXE) CAENQDC

clean	:
		/bin/rm -f $(OBJS) $(EXE) CAENQDC

$(EXE)	:	$(OBJS)
		/bin/rm -f $(EXE)
		$(CC) $(FLAGS) -o $(EXE) $(OBJS) $(DEPLIBS)

CAENQDC	:	$(OBJS_CAENQDC)
		/bin/rm -f CAENQDC
		$(CC) $(FLAGS) -o CAENQDC $(OBJS_CAENQDC) $(DEPLIBS)

$(OBJS)	:	makefile

%.o	:	%.c
		$(CC) $(COPTS) $(INCLUDEDIR) -c -o $@ $<


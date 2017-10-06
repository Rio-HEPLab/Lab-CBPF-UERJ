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

EXE	=	CAENCommon

CC	=	gcc

COPTS	=	-fPIC -DLINUX -Wall 
#COPTS	=	-g -fPIC -DLINUX -Wall 

FLAGS	=	-Wall -s
#FLAGS	=	-Wall

DEPLIBS	=       -lCAENVME -lcurses -lCAENV1x90SDK -lCAENqtpSDK

LIBS	=	

#INCLUDEDIR =	-I./ -I../../../SDK/CAENV1x90SDK/include
INCLUDEDIR =	-I. -I../include -I/opt/CAEN/CAENV1x90_1_1_linux/CAENVMEToolBox/SDK/CAENV1x90SDK/include -I/opt/CAEN/CAENqtp-1.2.1/SDK -I/opt/CAEN/CAENqtp-1.2.1/SDK/common/include/CAENVMELib

OBJS	=	console.o keyb.o user_settings_TDC.o user_settings_QDC.o main.o

#########################################################################

all	:	$(EXE)

clean	:
		/bin/rm -f $(OBJS) $(EXE)

$(EXE)	:	$(OBJS)
		/bin/rm -f $(EXE)
		$(CC) $(FLAGS) -o $(EXE) $(OBJS) $(DEPLIBS)

$(OBJS)	:	Makefile

%.o	:	%.c
		$(CC) $(COPTS) $(INCLUDEDIR) -c -o $@ $<


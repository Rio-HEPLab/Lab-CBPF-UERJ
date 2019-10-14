
LIB	=	libCAENqtpSDK.so.1.2

CC	=	gcc

COPTS	=	-fPIC -DLINUX -Wall 
#COPTS	=	-g -fPIC -DLINUX -Wall 

FLAGS	=	-shared -Wall,-soname -nostartfiles -s
#FLAGS	=	-shared -Wall,-soname -nostartfiles

DEPLIBS	=       -lCAENVME

LIBS	=	

INCLUDEDIR =	-I./ -I./common/include/CAENVMELib/ -I./common/src/

INSTALLDIR =	/usr/lib

OBJS	=	./common/src/cvt_board_commons.o cvt_V792.o 

INCLUDES =	./common/include/CAENVMELib/CAENVMElib.h ./common/include/CAENVMELib/CAENVMEtypes.h ./common/include/CAENVMELib/CAENVMEoslib.h

#########################################################################

all	:	$(LIB)

clean	:
		/bin/rm -f $(OBJS) $(LIB)

install	:	$(LIB)
		/bin/rm -f $(INSTALLDIR)/libCAENqtpSDK.so
		install $(LIB) $(INSTALLDIR)
		ln -sf $(INSTALLDIR)/$(LIB) $(INSTALLDIR)/libCAENqtpSDK.so
		ldconfig

$(LIB)	:	$(OBJS)
		/bin/rm -f $(BINLIB)
		$(CC) $(FLAGS) -o $(LIB) $(OBJS) $(DEPLIBS)

$(OBJS)	:	$(INCLUDES) makefile

%.o	:	%.c
		$(CC) $(COPTS) $(INCLUDEDIR) -c -o $@ $<


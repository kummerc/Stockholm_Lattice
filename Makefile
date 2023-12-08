#
# Copyright (c) 2013 Forschungszentrum Juelich
#
# Author(s): Dirk Pleiter
#
# This software is available to you under a choice of one of two
# licenses.  You may choose to be licensed under the terms of the GNU
# General Public License (GPL) Version 2, available from the file
# COPYING in the main directory of this source tree, or the
# OpenIB.org BSD license below:
#
#     Redistribution and use in source and binary forms, with or
#     without modification, are permitted provided that the following
#     conditions are met:
#
#      - Redistributions of source code must retain the above
#        copyright notice, this list of conditions and the following
#        disclaimer.
#
#      - Redistributions in binary form must reproduce the above
#        copyright notice, this list of conditions and the following
#        disclaimer in the documentation and/or other materials
#        provided with the distribution.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#---------------------------------------------------------------------------------------------------

TARGET   = uedu
SRC      = main.cpp common.cpp geom.cpp rng.cpp u.cpp
CC       = acpp
CFLAGS   = -O3 -DU_SIMPLE

.PHONY: clean veryclean dep

all:	${TARGET}

uedu:	${SRC:%.c=%.o}
	${CC} -o $@ $^ -lm

%.o:	%.c
	${CC} ${CFLAGS} -c $<

clean:
	-rm -f ${SRC:%.c=%.o}
	-rm -f depend.mk

veryclean: clean
	-rm -f ${TARGET}
	-rm -rf doc

doc:	${SRC} ${wildcard *.h} doxy.par
	doxygen doxy.par

dep:	${SRC}
	@${CC} ${CFLAGS} -M $^ > depend.mk
depend.mk:
	touch $@
include depend.mk

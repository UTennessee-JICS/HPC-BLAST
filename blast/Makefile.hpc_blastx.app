#  Copyright 2016 UTK JICS AACE
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#           http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
#
#  ===========================================================================
# $Id: Makefile.blastp.app 421234 2013-12-10 19:27:57Z camacho $

WATCHERS = camacho madden fongah2

APP = hpc_blastx
SRC = hpc_blastx_app
LIB_ = $(BLAST_INPUT_LIBS) $(BLAST_LIBS) $(OBJMGR_LIBS)
LIB = blast_app_util $(LIB_:%=%$(STATIC))

# De-universalize Mac builds to work around a PPC toolchain limitation
CFLAGS 	 = $(FAST_CXXFLAGS:ppc=i386) 
CXXFLAGS = $(FAST_CXXFLAGS:ppc=i386) 
LDFLAGS  = $(FAST_LDFLAGS:ppc=i386) 

CPPFLAGS = -DNCBI_MODULE=BLAST $(ORIG_CPPFLAGS)
LIBS = $(CMPRS_LIBS) $(DL_LIBS) $(NETWORK_LIBS) $(ORIG_LIBS)

REQUIRES = objects -Cygwin

PROJ_TAG = gbench

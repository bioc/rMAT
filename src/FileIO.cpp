/////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 2.1 of the License,
// or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
/////////////////////////////////////////////////////////////////


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "FileIO.h"

//////////////////////////////////////////////////////////////////////

#ifdef WIN32
#pragma warning(disable: 4996) // don't show deprecated warnings.
#include <winsock2.h>
#else
#include <sys/types.h>
#include <netinet/in.h>
#include <inttypes.h>
#endif

// Some machines (sparc) dont support unaligned memory access
// The mmaped files are chock full of unaligned accesses.
// When AFFY_NOUNALIGNED_MEM is defined we will do the 
// alignment in software.
// This feature can be enabled on intel for testing.

// This is here to enable this feature automaticly on the sparc and ppc.
// they cant do unaligned loads.
#ifdef __sparc__
#define AFFY_UNALIGNED_IN_SW
#endif
#ifdef __POWERPC__
#define AFFY_UNALIGNED_IN_SW
#endif

using namespace std;

void
ReadUInt32_I(IFSTREAM& instr, uint32_t& val) 
{
  uint32_t v=0;
  instr.read((char*)&v,sizeof(v));
  val=itohl(v);
}
void
ReadInt32_I(IFSTREAM& instr,int32_t& val) 
{
  ReadUInt32_I(instr,(uint32_t&)val);
}
void
ReadFloat_I(IFSTREAM& instr,float& val) 
{
  ReadUInt32_I(instr,(uint32_t&)val);
}
//
void
ReadUInt16_I(IFSTREAM& instr, uint16_t& val) 
{
  uint16_t v=0;
  instr.read((char*)&v,sizeof(v));
  val=itohs(v);
}
void
ReadInt16_I(IFSTREAM& instr,int16_t& val) 
{
  ReadUInt16_I(instr,(uint16_t&)val);
}

// No byte swapping needed.
void
ReadUInt8(IFSTREAM& instr, uint8_t& val) 
{
  instr.read((char*)&val,sizeof(val));
}
void
ReadInt8(IFSTREAM& instr,int8_t& val) 
{
  ReadUInt8(instr,(uint8_t&)val);
}

//====================
//
// Network byte order
//

void
ReadUInt32_N(IFSTREAM& instr, uint32_t& val) 
{
  uint32_t v=0;
  instr.read((char*)&v,sizeof(v));
  val=ntohl(v);
}
void
ReadInt32_N(IFSTREAM& instr,int32_t& val) 
{
  ReadUInt32_N(instr,(uint32_t&)val);
}
void
ReadFloat_N(IFSTREAM& instr,float& val) 
{
  ReadUInt32_N(instr,(uint32_t&)val);
}
//
void
ReadUInt16_N(IFSTREAM& instr, uint16_t& val) 
{
  uint16_t v=0;
  instr.read((char*)&v,sizeof(v));
  val=ntohs(v);
}
void
ReadInt16_N(IFSTREAM& instr, int16_t& val) 
{
  ReadUInt16_N(instr,(uint16_t&)val);
}

void
ReadFloatFromOldBPMAP_N(IFSTREAM &instr, float &fval)
{
#ifndef IS_BIG_ENDIAN
	instr.read((char *)&fval, FLOAT_SIZE);
	fval = (float) ntohl((uint32_t)fval);
#else
	int i1=0;
	instr.read((char *)&i1, INT32_SIZE);
 	i1=affy_swap32(i1);
	fval=(float)affy_swap32((uint32_t)(*(float *)&i1));
#endif
}

//==============================

// c char*
void
ReadFixedCString(IFSTREAM& instr, char* str, uint32_t len)
{
  instr.read(str,len);
}

void
ReadFixedUCString(IFSTREAM& instr, unsigned char* str, uint32_t len)
{
  instr.read((char *)str,len);
}

void
ReadCString_I(IFSTREAM& instr, char*& str)
{
  uint32_t slen;
  ReadUInt32_I(instr,slen);
  str=new char[slen+1];
  instr.read(str,slen);
  // ensure null -- this replaces the last char read.
  str[slen]=0; 
}

void
ReadCharacterArray(IFSTREAM& instr, char* str, uint32_t len)
{
  instr.read(str,len);
}

// The data file may or may not have the null terminator as
// part of the length. Using a STD::STRING function as the buffer to
// read directly into may result in the string with a length one more
// than is intended. The CHAR * buffer will contain an extra null that
// will get removed when copying to the output STD::STRING.
void
ReadFixedString(IFSTREAM& instr, string& str, uint32_t len)
{
  char *s = new char[len+1];
  instr.read(s,len);
  s[len]=0;
  str = s;
  delete[] s;
}

void
ReadString_I(IFSTREAM& instr, string& str)
{
  uint32_t len;
  ReadUInt32_I(instr,len);
  ReadFixedString(instr,str,len);
}

void
ReadUIntLenString_I(IFSTREAM& instr, std::string &s)
{
	uint32_t len;
	ReadUInt32_I(instr, len);
	ReadFixedString(instr, s, len);
}

void
ReadCString_N(IFSTREAM& instr, char*& str)
{
  uint32_t slen;
  ReadUInt32_N(instr,slen);
  str=new char[slen+1];
  instr.read(str,slen);
  // ensure null -- this replaces the last char read.
  str[slen]=0; 
}

void
ReadString_N(IFSTREAM& instr, string& str)
{
  uint32_t len;
  ReadUInt32_N(instr,len);
  ReadFixedString(instr,str,len);
}

void
ReadUIntLenString_N(IFSTREAM& instr, std::string &s)
{
	uint32_t len;
	ReadUInt32_N(instr, len);
	ReadFixedString(instr, s, len);
}

void 
ReadNextLine(IFSTREAM& instr,char* line,int len)
{
	strcpy(line,"");
	while (!instr.eof())
	{
		instr.getline(line,len);
		if (strlen(line)>0)
		{
			if (line[strlen(line)-1]=='\r')
				line[strlen(line)-1]='\0';
			if (strlen(line)>0)
				return;
		}
	}
}


//==============================
//

// There are two issues here: alignment and byte order.  
// Byte order can be handled with the "htoi" functions.
// Alignment has to be handled with functins.

// The "Get*" functions are for reading data from a mem-mapped file.
// OLD: Callers of these really should use "itoh*" directly.
//      Opps! They cant as "itoh" does not handle aligment.

// The signed versions
int32_t 
MmGetInt32_I(int32_t* ptr)
{
  return MmGetUInt32_I((uint32_t*)ptr);
}
int16_t 
MmGetInt16_I(int16_t* ptr)
{
  return MmGetUInt16_I((uint16_t*)ptr);
}
int32_t 
MmGetInt32_N(int32_t* ptr)
{
  return MmGetUInt32_N((uint32_t*)ptr);
}
int16_t 
MmGetInt16_N(int16_t* ptr)
{
  return MmGetUInt16_N((uint16_t*)ptr);
}
int8_t 
MmGetInt8(int8_t* ptr)
{
  return MmGetUInt8((uint8_t*)ptr);
}

//==========
// If unaligned accesses to memory are allowed, (intel)
// use these functions. 

#ifndef AFFY_UNALIGNED_IN_SW

//// 32
uint32_t
MmGetUInt32_I(uint32_t* ptr)
{
  return itohl(*ptr);
}
void
MmSetUInt32_I(uint32_t* ptr,uint32_t val)
{
  *(uint32_t*)ptr=htoil(val);
}
//
uint32_t
MmGetUInt32_N(uint32_t* ptr)
{
  return ntohl(*ptr);
}
void
MmSetUInt32_N(uint32_t* ptr,uint32_t val)
{
  *(uint32_t*)ptr=htonl(val);
}
//// 16
uint16_t
MmGetUInt16_I(uint16_t* ptr)
{
  return itohs(*ptr);
}
void
MmSetUInt16_I(uint16_t* ptr,uint16_t val)
{
  *(uint16_t*)ptr=htois(val);
}
//
uint16_t
MmGetUInt16_N(uint16_t* ptr)
{
  return ntohs(*ptr);
}
void
MmSetUInt16_N(uint16_t* ptr,uint16_t val)
{
  *(uint16_t*)ptr=htons(val);
}
//
uint8_t
MmGetUInt8(uint8_t* ptr)
{
  return *ptr; // no byte-swapping needed.
}
void
MmSetUInt8(uint8_t* ptr,uint8_t val)
{
  *ptr=val;
}
//
float
MmGetFloat_I(float* ptr)
{
  uint32_t v=itohl(*(int*)ptr);
  return *((float*)&v);
}
void
MmSetFloat_I(float* ptr,float val)
{
  *(uint32_t*)ptr=htoil(*(uint32_t*)&val);
}
float
MmGetFloat_N(float* ptr)
{
  uint32_t v=ntohl(*(int*)ptr);
  return *((float*)&v);
}
void
MmSetFloat_N(float* ptr,float val)
{
  *(uint32_t*)ptr=htonl(*(uint32_t*)&val);
}

float
MmGetFloatFromOldBPMAP_N(float *ptr)
{
	float fval;
#ifndef IS_BIG_ENDIAN
	fval = (float) ntohl((uint32_t)*ptr);
#else
	int i1=*(int *)ptr;
 	i1=affy_swap32(i1);
	fval=(float)affy_swap32((uint32_t)(*(float *)&i1));
#endif
	return fval;
}

#endif

// We dont have unaligned access to memory, fake it
// The conversion from and to little endian is done as part of
// the unaligned mem access

#ifdef AFFY_UNALIGNED_IN_SW

//// 32
uint32_t
MmGetUInt32_I(uint32_t* ptr)
{
  uint8_t* cptr=(uint8_t*)ptr;
  uint32_t val=0;
  val|=(*cptr++);
  val|=(*cptr++)<< 8;
  val|=(*cptr++)<<16;
  val|=(*cptr++)<<24;
  return val;
}
void
MmSetUInt32_I(uint32_t* ptr,uint32_t val)
{
  uint8_t* cptr=(uint8_t*)ptr;
  *cptr++=((val    )&0xFF);
  *cptr++=((val>> 8)&0xFF);
  *cptr++=((val>>16)&0xFF);
  *cptr++=((val>>24)&0xFF);
}

uint32_t
MmGetUInt32_N(uint32_t* ptr)
{
  uint8_t* cptr=(uint8_t*)ptr;
  uint32_t val=0;
  val|=(*cptr++)<<24;
  val|=(*cptr++)<<16;
  val|=(*cptr++)<< 8;
  val|=(*cptr++);
  return val;
}
void
MmSetUInt32_N(uint32_t* ptr,uint32_t val)
{
  uint8_t* cptr=(uint8_t*)ptr;
  // val=htonl(val); // with it in big order
  *cptr++=((val>>24)&0xFF);
  *cptr++=((val>>16)&0xFF);
  *cptr++=((val>> 8)&0xFF);
  *cptr++=((val    )&0xFF);
}

//// 16
uint16_t
MmGetUInt16_I(uint16_t* ptr)
{
  uint8_t* cptr=(uint8_t*)ptr;
  uint16_t val=0;
  val|=(*cptr++);
  val|=(*cptr++)<<8;
  return val;
}
void
MmSetUInt16_I(uint16_t* ptr,uint16_t val)
{
  uint8_t* cptr=(uint8_t*)ptr;
  *cptr++=((val    )&0xFF);
  *cptr++=((val>> 8)&0xFF);
}
uint16_t
MmGetUInt16_N(uint16_t* ptr)
{
  uint8_t* cptr=(uint8_t*)ptr;
  uint16_t val=0;
  val|=(*cptr++)<< 8;
  val|=(*cptr++);
  return val;
}
void
MmSetUInt16_N(uint16_t* ptr,uint16_t val)
{
  uint8_t* cptr=(uint8_t*)ptr;
  *cptr++=((val>>8)&0xFF);
  *cptr++=((val   )&0xFF);
}

//// 8
uint8_t
MmGetUInt8(uint8_t* ptr)
{
  return *(uint8_t*)ptr;
}
void
MmSetUInt8(uint8_t* ptr,uint8_t val)
{
  *(uint8_t*)ptr=val;
}
//
float
MmGetFloat_I(float* ptr)
{
  uint32_t val=MmGetUInt32_I((uint32_t*)ptr);
  return *((float*)&val);
}
void
MmSetFloat_I(float* ptr,float val)
{
  MmSetUInt32_I((uint32_t*)ptr,*(uint32_t*)&val);
}
float
MmGetFloat_N(float* ptr)
{
  uint32_t val=MmGetUInt32_N((uint32_t*)ptr);
  return *((float*)&val);
}

float
MmGetFloatFromOldBPMAP_N(float *ptr)
{
	float fval;
#ifndef IS_BIG_ENDIAN
	fval = (float) ntohl((uint32_t)MmGetFloat_I(ptr));
#else
	int i1=*(int *)ptr;
 	i1=affy_swap32(i1);
	fval=(float)affy_swap32((uint32_t)(*(float *)&i1)); // cast to float, then deref
#endif
	return fval;
}

void
MmSetFloat_N(float* ptr,float val)
{
  MmSetUInt32_N((uint32_t*)ptr,*(uint32_t*)&val);
}
#endif

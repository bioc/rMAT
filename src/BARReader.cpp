
/**Because there were some memory errors when running the code which ported Affymetrix Fusion SDK BPMAPFileData's code directly into R, I reimplemented it using its core code from ReadHeaderSection() and ReadDataSection().*/

//THE TERMS FROM AFFYMETRIX:

////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
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

#include <iostream>
#include <stdio.h>
#include "BARFileData.h"
#include "affy-base-types.h"
#include "FileIO.h"

#include <stdlib.h>

using namespace affxbar;
using namespace std;

extern "C"
{
#include <R.h>
#include <Rdefines.h>
#include <wchar.h>
#include <wctype.h>

/**Because there were some memory errors when running the code which ported Affymetrix Fusion SDK BPMAPFileData's code directly into R, I reimplemented it using its core code from ReadHeaderSection() and ReadDataSection().*/
  SEXP Parser (SEXP fileName)
  {
    const char *fname = CHAR (STRING_ELT (fileName, 0));
    // Open the file.
      std::ifstream instr;
     instr.open (fname, std::ios::in | std::ios::binary);
    
     // Check if open
    if (!instr)   
      {	
	cout << "Unable to open the file." << endl;
	return R_NilValue;     
	}
    
      // Magic number
      std::string magic;
      ReadFixedString (instr, magic, 8);
    
      // Version
    float m_Version;
    
   ReadFloat_N (instr, m_Version);
    
      // Number of sequendes
      int32_t cType;
    int m_NumberSequences;
    
   ReadInt32_N (instr, cType);
    
   m_NumberSequences = cType;
    
      // Columns
    int i = 0;
    
	ReadInt32_N (instr, cType);
    int m_NumberColumns;
    
	m_NumberColumns = cType;
    GDACFILES_BAR_DATA_TYPE_VECTOR m_ColumnTypes;
    
	m_ColumnTypes.resize (m_NumberColumns);
    
for (i = 0; i < m_NumberColumns; i++)     
      {
	
	ReadInt32_N (instr, cType);
	m_ColumnTypes[i] = (GDACFILES_BAR_DATA_TYPE) cType;      
}
    
      // Parameter
      std::string str;
    
	ReadInt32_N (instr, cType);
	    int m_NumberParameters;
    
	m_NumberParameters = cType;
	    TagValuePairTypeVector m_Parameters;
    
	m_Parameters.resize (m_NumberParameters);
    
	TagValuePairType param;
    
for (i = 0; i < m_NumberParameters; i++) 
      {
	
ReadString_N (instr, str);
m_Parameters[i].Tag = str;
	
ReadString_N (instr, str);
m_Parameters[i].Value = str;

//	cout << m_Parameters[i].Tag << endl;
//	cout << m_Parameters[i].Value << endl;
      
}
    
      // Determine the position of the start of the data
    int m_DataStartPosition = instr.tellg ();
    
      // Skip to the data section
      instr.seekg (m_DataStartPosition);
   
 
int totalDataPoints = 0;

    string *m_Name = new string[m_NumberSequences];
    string *m_GroupName = new string[m_NumberSequences];
    string *m_VersionSeq = new string[m_NumberSequences];
    int *m_NumberDataPoints = Calloc (m_NumberSequences, int);

    BarSequenceResultData ***tempStorage;
    tempStorage = new BarSequenceResultData **[m_NumberSequences];

//    cout << "m_NumberSequences is " << m_NumberSequences << endl;
    int curProbePos = 0;

    for (int i = 0; i < m_NumberSequences; i++)
      {
//cout<<"PROCESSING " <<i <<endl;
	ReadString_N (instr, m_Name[i]);
//	cout << m_Name[i] << endl;
	bool bVersion2 = ((int) (m_Version + 0.1) == 2);
	if (bVersion2)
	  {
	    
	  ReadString_N (instr, m_GroupName[i]);
//	    cout << m_GroupName[i] << endl;
	 }
	
ReadString_N (instr, m_VersionSeq[i]);

	
	if (bVersion2)
	   {
		int32_t nParams = 0;	    
		ReadInt32_N (instr, nParams);
//	    cout << "nParams is " << nParams << endl;
	    string tag;
	    string val;
	    
	for (int iParam = 0; iParam < nParams; iParam++)     
	      {
		ReadString_N (instr, tag);		
		ReadString_N (instr, val);
//                              m_Results[i].AddParameter(tag, val);
//		cout << tag << ", " << val << endl;  
		} 
	}
	
	int32_t cType;	
	ReadInt32_N (instr, cType);
	m_NumberDataPoints[i] = cType;
//	cout << "m_NumberDataPoints[i] is " << m_NumberDataPoints[i] << endl;
	
//              m_Results[i].m_NumberColumns = m_NumberColumns;
//              m_Results[i].m_pColumnTypes = &m_ColumnTypes;
//              m_Results[i].m_DataStartPosition = instr.tellg();
	  totalDataPoints += m_NumberDataPoints[i];
//	cout << "totalDataPoints is " << totalDataPoints << endl;

	tempStorage[i] = new BarSequenceResultData *[m_NumberDataPoints[i]];
	for (int k = 0; k < m_NumberDataPoints[i]; k++)
	  {
	    tempStorage[i][k] = new BarSequenceResultData[m_NumberColumns];
	  }


	string chrName = m_Name[i];

          int curChr = atoi (chrName.substr (3, 2).c_str ());

//	cout << "curChr is " << curChr << endl;
//
//         cout<<"Curchr is " <<curChr <<endl;
	for (int ipt = 0; ipt < m_NumberDataPoints[i]; ipt++)
	  {

	    tempStorage[i][ipt][0].iValue = curChr;

	    int32_t integerValue;
	    ReadInt32_N (instr, integerValue);
//	    cout << "integerValue " << integerValue << endl;
	    tempStorage[i][ipt][1].iValue = integerValue;

	    for(int k=2; k<=m_NumberColumns; k++){
		    float floatValue;
		    ReadFloat_N (instr, floatValue);
		    tempStorage[i][ipt][k].fValue = floatValue;
//		    cout << "floatValue " << floatValue << endl;
	     }

	  }

      }



//    cout << "DONE" << endl;




	SEXP hit_names, hitList;
	PROTECT (hit_names = allocVector (STRSXP, m_NumberColumns+1));
	PROTECT (hitList = allocVector (VECSXP, m_NumberColumns+1));

	SEXP chromosomeR, positionR;
	SEXP *signalR = new SEXP[m_NumberColumns-1];
	int *p_chromosome = NULL, *p_position = NULL;
	double **p_signal = Calloc(m_NumberColumns-1, double *);


	PROTECT (chromosomeR = NEW_INTEGER (totalDataPoints));
	p_chromosome = INTEGER_POINTER (chromosomeR);
	SET_STRING_ELT (hit_names, 0, mkChar ("chr"));
	SET_VECTOR_ELT (hitList, 0, chromosomeR);


	PROTECT (positionR = NEW_INTEGER (totalDataPoints));
	p_position = INTEGER_POINTER (positionR);
	SET_STRING_ELT (hit_names, 1, mkChar ("pos"));
	SET_VECTOR_ELT (hitList, 1, positionR);


	for(int k=0; k< (m_NumberColumns-1); k++){
		PROTECT (signalR[k] = NEW_NUMERIC (totalDataPoints));
		p_signal[k] = NUMERIC_POINTER (signalR[k]);
		char str[10];
//		char str2[4];
//		strcpy (str,"signal");
//		strcat (str,"signal");

const int MAX_SIZE = 255;
char buf[MAX_SIZE];
snprintf( str, MAX_SIZE, "signal%d", (k+1) );
//		itoa( (k+1), str2, 10);	
		//	char *  itoa ( int value, char * str, int base );
		//		snprintf( buf, MAX_SIZE, "FPS: %d", myFps );
//cout<<str <<endl;
		SET_STRING_ELT (hit_names, k+2, mkChar ( str));
		SET_VECTOR_ELT (hitList, k+2, signalR[k]);
	}


	curProbePos=0;
	for (int curSeq = 0; curSeq < m_NumberSequences; curSeq++)
	  {
	    for (int ipt = 0; ipt < m_NumberDataPoints[curSeq]; ipt++)
	      {
		p_chromosome[curProbePos] =  tempStorage[curSeq][ipt][0].iValue;
		p_position[curProbePos] = tempStorage[curSeq][ipt][1].iValue;
		for(int k=0; k<(m_NumberColumns-1); k++){
			p_signal[k][curProbePos] = tempStorage[curSeq][ipt][2+k].fValue;
		}
			curProbePos++;
	      }
	  }



	setAttrib (hitList, R_NamesSymbol, hit_names);
	UNPROTECT (3+m_NumberColumns);

    // Close the file
      instr.close ();
  
	return hitList;
  }


  /**ParsedMATBar**/
  SEXP ParseBar (SEXP fileName)
  {

    const char *fname = CHAR (STRING_ELT (fileName, 0));
    TagValuePairType param;
    CBARFileData bar;
    bar.SetFileName (fname);
    if (bar.Exists ())
      {
	bar.GetFileName ();
	bar.ReadHeader ();
	bar.Read ();

	int NUMSEQ = bar.GetNumberSequences ();
	int *seqLength = new int[NUMSEQ];
	int totalDataPoints = 0;

	for (int i = 0; i < NUMSEQ; i++)
	  {
	    CGDACSequenceResultItem res;
	    bar.GetResults (i, res);
	    seqLength[i] = res.GetNumberDataPoints ();
	    totalDataPoints += seqLength[i];
	  }


	SEXP hit_names, hitList;
	PROTECT (hit_names = allocVector (STRSXP, 4));
	PROTECT (hitList = allocVector (VECSXP, 4));

	SEXP chromosomeR, positionR, matScoreR, pValueR;
	int *p_chromosome = NULL, *p_position = NULL;
	double *p_matScore = NULL, *p_pValue = NULL;


	PROTECT (chromosomeR = NEW_INTEGER (totalDataPoints));
	p_chromosome = INTEGER_POINTER (chromosomeR);
	SET_STRING_ELT (hit_names, 0, mkChar ("chr"));
	SET_VECTOR_ELT (hitList, 0, chromosomeR);


	PROTECT (positionR = NEW_INTEGER (totalDataPoints));
	p_position = INTEGER_POINTER (positionR);
	SET_STRING_ELT (hit_names, 1, mkChar ("pos"));
	SET_VECTOR_ELT (hitList, 1, positionR);

	PROTECT (matScoreR = NEW_NUMERIC (totalDataPoints));
	p_matScore = NUMERIC_POINTER (matScoreR);
	SET_STRING_ELT (hit_names, 2, mkChar ("MATScore"));
	SET_VECTOR_ELT (hitList, 2, matScoreR);

	PROTECT (pValueR = NEW_NUMERIC (totalDataPoints));
	p_pValue = NUMERIC_POINTER (pValueR);
	SET_STRING_ELT (hit_names, 3, mkChar ("pValue"));
	SET_VECTOR_ELT (hitList, 3, pValueR);

	BarSequenceResultData ***tempStorage;
	tempStorage = new BarSequenceResultData **[NUMSEQ];
	int curProbePos = 0;

	for (int curSeq = 0; curSeq < NUMSEQ; curSeq++)
	  {
	    CGDACSequenceResultItem res;
	    bar.GetResults (curSeq, res);

	    tempStorage[curSeq] =
	      new BarSequenceResultData *[res.GetNumberDataPoints ()];
	    for (int i = 0; i < res.GetNumberDataPoints (); i++)
	      {
		tempStorage[curSeq][i] = new BarSequenceResultData[4];
	      }

	    string chrName = res.GetName ();
	    int curChr = atoi (chrName.substr (3, 2).c_str ());
	    /*cout<<"Curchr is " <<curChr <<endl; */
	    for (int ipt = 0; ipt < res.GetNumberDataPoints (); ipt++)
	      {
		tempStorage[curSeq][ipt][0].iValue = curChr;
		p_chromosome[curProbePos] =
		  tempStorage[curSeq][ipt][0].iValue;
		res.GetData (ipt, 0, tempStorage[curSeq][ipt][1]);
		p_position[curProbePos] = tempStorage[curSeq][ipt][1].iValue;
		res.GetData (ipt, 1, tempStorage[curSeq][ipt][2]);
		p_matScore[curProbePos] = tempStorage[curSeq][ipt][2].fValue;
		res.GetData (ipt, 2, tempStorage[curSeq][ipt][3]);
		p_pValue[curProbePos] = tempStorage[curSeq][ipt][3].fValue;
		curProbePos++;
	      }
	  }
	setAttrib (hitList, R_NamesSymbol, hit_names);
	UNPROTECT (6);
	return hitList;
      }
    else
      {
	cout << "File does not exists..." << endl;
      }
    return R_NilValue;
  }













































/**** OBSOLETE ****/
  /**ParsedMATBar**/
  SEXP ParseMATBar (SEXP fileName)
  {

    const char *fname = CHAR (STRING_ELT (fileName, 0));
    TagValuePairType param;
    CBARFileData bar;
    bar.SetFileName (fname);
    if (bar.Exists ())
      {
	bar.GetFileName ();
	bar.ReadHeader ();
	bar.Read ();

	int NUMSEQ = bar.GetNumberSequences ();
	int *seqLength = new int[NUMSEQ];
	int totalDataPoints = 0;

	for (int i = 0; i < NUMSEQ; i++)
	  {
	    CGDACSequenceResultItem res;
	    bar.GetResults (i, res);
	    seqLength[i] = res.GetNumberDataPoints ();
	    totalDataPoints += seqLength[i];
	  }


	SEXP hit_names, hitList;
	PROTECT (hit_names = allocVector (STRSXP, 5));
	PROTECT (hitList = allocVector (VECSXP, 5));

	SEXP chromosomeR, positionR, matScoreR, regionR, pValueR;
	int *p_chromosome = NULL, *p_position = NULL, *p_region = NULL;
	double *p_matScore = NULL, *p_pValue = NULL;


	PROTECT (chromosomeR = NEW_INTEGER (totalDataPoints));
	p_chromosome = INTEGER_POINTER (chromosomeR);
	SET_STRING_ELT (hit_names, 0, mkChar ("chr"));
	SET_VECTOR_ELT (hitList, 0, chromosomeR);


	PROTECT (positionR = NEW_INTEGER (totalDataPoints));
	p_position = INTEGER_POINTER (positionR);
	SET_STRING_ELT (hit_names, 1, mkChar ("pos"));
	SET_VECTOR_ELT (hitList, 1, positionR);

	PROTECT (matScoreR = NEW_NUMERIC (totalDataPoints));
	p_matScore = NUMERIC_POINTER (matScoreR);
	SET_STRING_ELT (hit_names, 2, mkChar ("MATScore"));
	SET_VECTOR_ELT (hitList, 2, matScoreR);

	PROTECT (pValueR = NEW_NUMERIC (totalDataPoints));
	p_pValue = NUMERIC_POINTER (pValueR);
	SET_STRING_ELT (hit_names, 3, mkChar ("pValue"));
	SET_VECTOR_ELT (hitList, 3, pValueR);

	BarSequenceResultData ***tempStorage;
	tempStorage = new BarSequenceResultData **[NUMSEQ];
	int curProbePos = 0;

	for (int curSeq = 0; curSeq < NUMSEQ; curSeq++)
	  {
	    CGDACSequenceResultItem res;
	    bar.GetResults (curSeq, res);

	    tempStorage[curSeq] =
	      new BarSequenceResultData *[res.GetNumberDataPoints ()];
	    for (int i = 0; i < res.GetNumberDataPoints (); i++)
	      {
		tempStorage[curSeq][i] = new BarSequenceResultData[4];
	      }

	    string chrName = res.GetName ();
	    int curChr = atoi (chrName.substr (3, 2).c_str ());
	    /*cout<<"Curchr is " <<curChr <<endl; */
	    for (int ipt = 0; ipt < res.GetNumberDataPoints (); ipt++)
	      {
		tempStorage[curSeq][ipt][0].iValue = curChr;
		p_chromosome[curProbePos] =
		  tempStorage[curSeq][ipt][0].iValue;
		res.GetData (ipt, 0, tempStorage[curSeq][ipt][1]);
		p_position[curProbePos] = tempStorage[curSeq][ipt][1].iValue;
		res.GetData (ipt, 1, tempStorage[curSeq][ipt][2]);
		p_matScore[curProbePos] = tempStorage[curSeq][ipt][2].fValue;
		res.GetData (ipt, 2, tempStorage[curSeq][ipt][3]);
		p_pValue[curProbePos] = tempStorage[curSeq][ipt][3].fValue;
		curProbePos++;
	      }
	  }
	setAttrib (hitList, R_NamesSymbol, hit_names);
	UNPROTECT (6);
	return hitList;
      }
    else
      {
	cout << "File does not exists..." << endl;
      }
    return R_NilValue;
  }

/**** OBSOLETE ****/
  /**ParseNormalizeBar**/
  SEXP ParseNormalizeBar (SEXP fileName)
  {

    const char *fname = CHAR (STRING_ELT (fileName, 0));
    TagValuePairType param;
    CBARFileData bar;
    bar.SetFileName (fname);
    if (bar.Exists ())
      {
	bar.GetFileName ();
	bar.ReadHeader ();
	bar.Read ();

	int NUMSEQ = bar.GetNumberSequences ();
	int *seqLength = new int[NUMSEQ];
	int totalDataPoints = 0;

	for (int i = 0; i < NUMSEQ; i++)
	  {
	    CGDACSequenceResultItem res;
	    bar.GetResults (i, res);
	    seqLength[i] = res.GetNumberDataPoints ();
	    totalDataPoints += seqLength[i];
	  }

	/*
	   cout<<"Num data point is " <<res.GetNumberDataPoints() <<endl;
	   cout<<"Number of columns is " <<res.GetNumberColumns() <<endl;
	   cout <<"Number of parameters is " << res.GetNumberParameters() <<endl;
	   IGNORING THE SEQUENCE SPECIFIC PARAMETER.. USING THE GLOBAL
	   for(int i=0; i<bar.GetNumberParameters(); i++){
	   cout<<  bar.GetParameter(i).Tag <<"\t";
	   }
	 */


	SEXP hit_names, hitList;
	PROTECT (hit_names = allocVector (STRSXP, 3));
	PROTECT (hitList = allocVector (VECSXP, 3));

	SEXP chromosomeR, positionR, signalR;
	int *p_chromosome = NULL, *p_position = NULL;
	double *p_signal = NULL;

	PROTECT (chromosomeR = NEW_INTEGER (totalDataPoints));
	p_chromosome = INTEGER_POINTER (chromosomeR);
	SET_STRING_ELT (hit_names, 0, mkChar ("chr"));
	SET_VECTOR_ELT (hitList, 0, chromosomeR);

	PROTECT (positionR = NEW_INTEGER (totalDataPoints));
	p_position = INTEGER_POINTER (positionR);
	SET_STRING_ELT (hit_names, 1, mkChar ("pos"));
	SET_VECTOR_ELT (hitList, 1, positionR);

	PROTECT (signalR = NEW_NUMERIC (totalDataPoints));
	p_signal = NUMERIC_POINTER (signalR);
	SET_STRING_ELT (hit_names, 2, mkChar ("signal"));
	SET_VECTOR_ELT (hitList, 2, signalR);


	BarSequenceResultData ***tempStorage;
	tempStorage = new BarSequenceResultData **[NUMSEQ];
	int curProbePos = 0;
	cout << "STILL OKAY !!\n";
	cout << "NUMSEQ IS " << NUMSEQ << endl;
	for (int curSeq = 0; curSeq < NUMSEQ; curSeq++)
	  {
	    CGDACSequenceResultItem res;
	    bar.GetResults (curSeq, res);

	    tempStorage[curSeq] =
	      new BarSequenceResultData *[res.GetNumberDataPoints ()];
	    cout << "NUM DATA POINTS IS " << res.
	      GetNumberDataPoints () << endl;
	    for (int i = 0; i < res.GetNumberDataPoints (); i++)
	      {
		tempStorage[curSeq][i] = new BarSequenceResultData[2];
	      }

	    string chrName = res.GetName ();
	    int curChr = atoi (chrName.substr (3, 2).c_str ());
	    /*cout<<"Curchr is " <<curChr <<endl; */
	    for (int ipt = 0; ipt < res.GetNumberDataPoints (); ipt++)
	      {
		tempStorage[curSeq][ipt][0].iValue = curChr;
		p_chromosome[curProbePos] =
		  tempStorage[curSeq][ipt][0].iValue;
		res.GetData (ipt, 0, tempStorage[curSeq][ipt][1]);
		p_position[curProbePos] = tempStorage[curSeq][ipt][1].iValue;
//              res.GetData (ipt, 1, tempStorage[curSeq][ipt][2]);
//              p_signal[curProbePos] = tempStorage[curSeq][ipt][2].fValue;
		curProbePos++;
	      }

	    cout << "FINISH WITHOUT ERROR " << endl;
	  }
	setAttrib (hitList, R_NamesSymbol, hit_names);
	UNPROTECT (5);
	return hitList;
      }
    else
      {
	cout << "File does not exists..." << endl;
      }


    return R_NilValue;
  }

}

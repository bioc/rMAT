/**Because there were some memory errors when running the code which ported Affymetrix Fusion SDK BPMAPFileWriter's code directly into R, I reimplemented it using its core code from Save().*/

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

//
#include <iostream>
#include<fstream>
#include <stdio.h>

#include "FileWriter.h"
#include "BARFileWriter.h"
#include "BARFileData.h"


using namespace affxbarwriter;
using namespace affxbar;
using namespace std;

extern "C"
{
#include <R.h>
#include <Rdefines.h>
#include <wchar.h>
#include <wctype.h>
#include <string>




/**assume sorted by chromosome number, then by position **/
  SEXP testBAR (SEXP fileName, SEXP barNameString, SEXP MATScores, SEXP pValue, SEXP position, SEXP chNos, SEXP numSeq)
  {

    const char *fname = CHAR (STRING_ELT (fileName, 0));
  ofstream m_NewBarFile; 
	m_NewBarFile.open(fname, std::ios::out | std::ios::binary);

	// Magic number
	char magic[9];
	snprintf(magic,sizeof(magic), "barr%c%c%c%c", '\r', '\n', '\032', '\n');


WriteFixedString(m_NewBarFile, std::string(magic), 8);


float BAR_VERSION = 2.0f;
WriteFloat_N(m_NewBarFile, BAR_VERSION);	

int NUM_SEQ = 101;
WriteInt32_N(m_NewBarFile, NUM_SEQ);


const char *genomeName = CHAR (STRING_ELT (barNameString, 0));
WriteString_N(m_NewBarFile, genomeName);

	m_NewBarFile.close();
  

  }




      /**assume sorted by chromosome number, then by position **/
  SEXP WriteBAR (SEXP fileName, SEXP barNameString, SEXP MATScores, SEXP pValue, SEXP position, SEXP chNos, SEXP numSeq)
  {

    int signal2Length = length(pValue);

    int NUM_SEQ = INTEGER(numSeq)[0];
    int *lengthList = new int[NUM_SEQ];
    int curListPtr = 0, startListPtr = 0;
    int curChromosome = INTEGER(chNos)[0];
    int i;
 
	
    return R_NilValue;

  }




/**OBSOLETE --- MEMORY PROBLEM ON SOME COMPUTERS ***/
/**assume sorted by chromosome number, then by position **/
  SEXP WriteMATBAR (SEXP fileName, SEXP barNameString, SEXP MATScores,
		    SEXP pValue, SEXP position, SEXP chNos, SEXP numSeq)
  {

    int NUM_SEQ = INTEGER (numSeq)[0];
    int *lengthList = new int[NUM_SEQ];
    int curListPtr = 0, startListPtr = 0;
    int curChromosome = INTEGER (chNos)[0];
    int i;
    for (i = 0; i < length (MATScores); i++)
      {
	if (curChromosome != INTEGER (chNos)[i])
	  {
	    lengthList[curListPtr] = i - startListPtr;
	    startListPtr = i;
	    curListPtr++;
	    curChromosome = INTEGER (chNos)[i];
	  }
      }
    lengthList[curListPtr] = i - startListPtr;

    /*
    for(i=0; i < NUM_SEQ; i++){
	cout<< lengthList[i] << endl;
    }
    */

    const char *fname = CHAR (STRING_ELT (fileName, 0));
    CBARFileWriter outFile;

    outFile.SetFileName (fname);
    outFile.SetNumberSequences (NUM_SEQ);	/* number of results */
    outFile.AddColumn (BAR_DATA_INTEGER);
    outFile.AddColumn (BAR_DATA_FLOAT);
    outFile.AddColumn (BAR_DATA_FLOAT);
    /*outFile.AddColumn (BAR_DATA_FLOAT); // cannot have more than 2 columns */

    outFile.AddAlgorithmParameter ("position",
				   "the position of the probe along the chromosome");
    outFile.AddAlgorithmParameter ("Standardized MATScores",
				   "the computed MATScore of probe i / max MATScore of all probes x 100%");
    outFile.AddAlgorithmParameter ("pValue", "");
    /*outFile.AddAlgorithmParameter ("regions",
				   "0 denoted unenriched, nonzero denoted enriched"); */
    BarSequenceResultData ***tempStorage;
    tempStorage = Calloc (NUM_SEQ, BarSequenceResultData **);
    int curProbePos = 0;
    for (int curSeq = 0; curSeq < NUM_SEQ; curSeq++)
      {
	CGDACSequenceResultItem *result = outFile.GetResultsPtr (curSeq);

	char *chrStr = Calloc (20, char);
	char chNoStr[10] = "";
	sprintf (chNoStr, "%d", INTEGER (chNos)[curProbePos]);
	strcpy (chrStr, "chr");
	strcat (chrStr, chNoStr);

	result->SetName (chrStr);
	const char *BARNAME = CHAR (STRING_ELT (barNameString, 0));
	result->SetGroupName (BARNAME);
	/*result->SetVersion ("1.0");*/
	result->SetNumberDataPoints (lengthList[curSeq]);

	/*tempStorage[curSeq] = new BarSequenceResultData *[lengthList[curSeq]]; */
	tempStorage[curSeq] =  Calloc (lengthList[curSeq], BarSequenceResultData *);
	for (i = 0; i < lengthList[curSeq]; i++)
	  {
	    /*tempStorage[curSeq][i] = new BarSequenceResultData[4];*/
	    tempStorage[curSeq][i] = Calloc (3, BarSequenceResultData);
	  }

	for (i = 0; i < lengthList[curSeq]; i++)
	  {

	    tempStorage[curSeq][i][0].iValue =
	      INTEGER (position)[curProbePos];
	    tempStorage[curSeq][i][1].fValue = REAL (MATScores)[curProbePos];
	    tempStorage[curSeq][i][2].fValue = REAL (pValue)[curProbePos];
	    /*tempStorage[curSeq][i][3].fValue = REAL (region)[curProbePos];*/

	    result->SetDataPoint (i, 0, tempStorage[curSeq][i][0]);
	    result->SetDataPoint (i, 1, tempStorage[curSeq][i][1]);
	    result->SetDataPoint (i, 2, tempStorage[curSeq][i][2]);
	    /*result->SetDataPoint (i, 3, tempStorage[curSeq][i][3]);*/

	    curProbePos++;
	  }

      }		/* end of for */
    outFile.CreateNewFile ();
    outFile.Save ();
    return R_NilValue;

  }



/**OBSOLETE --- MEMORY PROBLEM ON SOME COMPUTERS ***/
/**assume sorted by chromosome number, then by position **/
  SEXP WriteNormalizedBAR (SEXP fileName, SEXP barNameString, SEXP normData,
			   SEXP position, SEXP chNos, SEXP numSeq)
  {

    int NUM_SEQ = INTEGER (numSeq)[0];
    int *lengthList = new int[NUM_SEQ];
    int curListPtr = 0, startListPtr = 0;
    int curChromosome = INTEGER (chNos)[0];
    int i;
    for (i = 0; i < length (normData); i++)
      {
	if (curChromosome != INTEGER (chNos)[i])
	  {
	    lengthList[curListPtr] = i - startListPtr;
	    startListPtr = i;
	    curListPtr++;
	    curChromosome = INTEGER (chNos)[i];
	  }
      }
    lengthList[curListPtr] = i - startListPtr;

   /*
    for(i=0; i < NUM_SEQ; i++){
      cout<< lengthList[i] << endl;
    }
   */

    const char *fname = CHAR (STRING_ELT (fileName, 0));

    CBARFileWriter outFile;
    outFile.SetFileName (fname);
    outFile.SetNumberSequences (NUM_SEQ);	/* number of results */
    outFile.AddColumn (BAR_DATA_INTEGER);
    outFile.AddColumn (BAR_DATA_FLOAT);
    outFile.AddAlgorithmParameter ("position",
				   "the position of the probe along the chromosome");
    outFile.AddAlgorithmParameter ("signal", "normalized signal");

    BarSequenceResultData ***tempStorage;
    tempStorage = new BarSequenceResultData **[NUM_SEQ];
    int curProbePos = 0;
    for (int curSeq = 0; curSeq < NUM_SEQ; curSeq++)
      {
	CGDACSequenceResultItem *result = outFile.GetResultsPtr (curSeq);

	char *chrStr = Calloc (20, char);
	char chNoStr[10] = "";
	sprintf (chNoStr, "%d", INTEGER (chNos)[curProbePos]);
	strcpy (chrStr, "chr");
	strcat (chrStr, chNoStr);

	result->SetName (chrStr);
	const char *BARNAME = CHAR (STRING_ELT (barNameString, 0));
	result->SetGroupName (BARNAME);
	/*result->SetVersion ("1.0");*/
	result->SetNumberDataPoints (lengthList[curSeq]);
	/*
	cout << "length list is " << lengthList[curSeq] << endl;
	int i = 0;
	*/

	tempStorage[curSeq] = new BarSequenceResultData *[lengthList[curSeq]];
	for (i = 0; i < lengthList[curSeq]; i++)
	  {
	    tempStorage[curSeq][i] = new BarSequenceResultData[2];
	  }

	int addPos = 0;
	for (i = 0; i < lengthList[curSeq]; i++)
	  {

	    tempStorage[curSeq][i][0].iValue =
	      INTEGER (position)[curProbePos];
	    tempStorage[curSeq][i][1].fValue = REAL (normData)[curProbePos];

	    result->SetDataPoint (i, 0, tempStorage[curSeq][i][0]);
	    result->SetDataPoint (i, 1, tempStorage[curSeq][i][1]);

	    curProbePos++;
	  }

      } /*end of for*/
    outFile.CreateNewFile ();
    outFile.Save ();
    return R_NilValue;

  }

}

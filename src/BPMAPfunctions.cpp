/*
*Last Edited: Aug 21, 2007
*
*
* Author: Charles Cheung and Raphael Gottardo
*
*/

#include "BPMAPFileData.h"
#include <iostream>
#include <string>

using namespace std;
using namespace affxbpmap;


#include <R.h>
#include <Rdefines.h>
#include <wchar.h>
#include <wctype.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern "C" {



SEXP convertSeqToChNo(SEXP seqNum, SEXP cTableSeq, SEXP cTableChNo)
{

  int nProbes = length(seqNum);
  int cTableRows = length(cTableSeq);

  SEXP chNo;
  int *p_chNo = NULL;
  PROTECT(chNo = NEW_INTEGER(nProbes));
  p_chNo = INTEGER_POINTER(chNo);

  int cacheSeq = -1;
  int cacheChNo = -1;
  int foundIndex;
  int seqNum_i;

for(int i=0; i<nProbes; i++)
{
	seqNum_i = INTEGER(seqNum)[i];
  //The seqNum of the current probe is the same as that of the last probe, so we don't need to query the table for the chNo because they are the same
	if(seqNum_i== cacheSeq)
		{
			p_chNo[i] = cacheChNo;
		}
	//Perform a search on the table for the chromosome number for that sequence number.
	else{
			foundIndex = -1;
			for(int k=0; k< cTableRows; k++)
			{
			  if(INTEGER(cTableSeq)[k] == seqNum_i)
			  {
				foundIndex = k;
			  }
			}
			if(foundIndex == -1)
			{
				cout<<"Out of range.\n" <<endl;
				cout<<"Sequence number: " <<seqNum_i <<" is not found.\n" <<endl;
				return R_NilValue;
			}
 			  cacheSeq = seqNum_i;
			  cacheChNo = INTEGER(cTableChNo)[foundIndex];
			  p_chNo[i] = cacheChNo;
		}
}


 UNPROTECT(1);
   return chNo;

}

/*Last Modified on June 28, 2007 */
/*Added on November 3, 2006 to parse all BPMAP sequences in the BPMAP file using the Affymetrix Fusion SDK */
/*Essentially a modification of readAllBPMAPSeq */
/*Modifications: 1. read the list of sequences*/

SEXP readBPMAPSeq(SEXP fileName, SEXP rangelist,  SEXP readPM, SEXP readMM, SEXP readProbeLength, SEXP readPMProbeSeq, SEXP readMatchScore, SEXP readPosition, SEXP readTopStrand, SEXP verbose, SEXP maxRange, SEXP minRange)
{

/*checking inputs are of the right types:*/
if(!isString(fileName)){ cout << "Filename is not a string. Parsing aborted.\n"; return R_NilValue;}
if( ! (isNumeric(rangelist) ||rangelist ==R_NilValue) ){ cout << "rangelist is not a vector of integers. Parsing aborted.\n"; return R_NilValue;}
if(!isLogical(readPM)){ cout << "readPM is not a boolean (logical). Parsing aborted.\n"; return R_NilValue;}
if(!isLogical(readMM)){ cout << "readMM is not a boolean (logical). Parsing aborted.\n"; return R_NilValue;}
if(!isLogical(readProbeLength)){ cout << "readProbeLength is not a boolean (logical). Parsing aborted.\n"; return R_NilValue;}
if(!isLogical(readPMProbeSeq)){ cout << "readPMProbeSeq is not a boolean (logical). Parsing aborted.\n"; return R_NilValue;}
if(!isLogical(readMatchScore)){ cout << "readMatchScore is not a boolean (logical). Parsing aborted.\n"; return R_NilValue;}
if(!isLogical(readPosition)){ cout << "readPosition is not a boolean (logical). Parsing aborted.\n"; return R_NilValue;}

if(!isLogical(readTopStrand)){ cout << "readTopStrand is not a boolean (logical). Parsing aborted.\n"; return R_NilValue;}
if(!isLogical(verbose)){ cout << "verbose is not a boolean (logical). Parsing aborted.\n"; return R_NilValue;}

/*saving out the data to pointers-------------------------------*/
int *isVerbose = LOGICAL(verbose);
int isPM = *LOGICAL(readPM);
int isMM = *LOGICAL(readMM);
int isProbeLength = *LOGICAL(readProbeLength);
int isPMProbeSeq = *LOGICAL(readPMProbeSeq);
int isMatchScore = *LOGICAL(readMatchScore);
int isPosition = *LOGICAL(readPosition);
int isTopStrand = *LOGICAL(readTopStrand);
/*---------------------------------------------------------------------*/

CBPMAPFileData data;
data.SetFileName( CHAR(STRING_ELT(fileName, 0)) );

/*Check whether we can read the header, thereby confirming whether it is a BPMAP file. */
if( data.ReadHeader()==false){
       cout<<"Fail to read header\n" <<endl;
       return R_NilValue;
}

int totalSeqNumber = data.GetNumberSequences();
/*Check whether the range is correct*/
if(( *INTEGER(maxRange) > totalSeqNumber)||(*INTEGER(minRange) < 0))
{
	cout<<"Range is out of bound\n" <<endl;
	return R_NilValue;
}

/*Loop through to determine the TOTALHITS*/
int TOTALHITS=0;
int lengthList;
int *rangeListArray;
if(rangelist!=R_NilValue){
	lengthList = length(rangelist);
	rangeListArray = new int[lengthList];
	for(int i=0; i<lengthList; i++){
		rangeListArray[i] = INTEGER(rangelist)[i]; /*current sequence number to read ;*/
	}
}
else{  /*we have to read all */
	lengthList =  data.GetNumberSequences(); /*get the total number of sequences stored in the file*/
	rangeListArray = new int[lengthList];
	for(int i=0; i<lengthList; i++){
		rangeListArray[i] = i;
	}
}


if(*isVerbose){
	cout<<"Number of sequences to parse is " <<lengthList <<endl;
}

/*read the content of the data*/
if(data.Read()==false){
    cout<<"Fail to read\n"  <<endl;
    return R_NilValue;
}

int cur=-1;
int numSeqToParse = 0; /*a counter from the first sequence to the last*/
while(numSeqToParse<lengthList)
{
/*	cur = (int)REAL(rangelist)[numSeqToParse]; //current sequence number to read */
	cur = rangeListArray[numSeqToParse]; /*current sequence number to read*/
	  /*sequence's description            */
	CGDACSequenceItem seq;
	data.GetSequenceItem(cur, seq);  /*get sequence information*/
	TOTALHITS=TOTALHITS+ seq.GetNumberHits();
	numSeqToParse++;
}

if(*isVerbose){
	cout<<"total number of sequences to read is " <<TOTALHITS<<endl;
	cout <<endl;
	cout<<"Start reading..."<<endl;
}


  /*the data structure (a list) and the accompanying header for the list to return in R*/
  SEXP hit_names, hitList;

  /*vectors stored into R */
  SEXP PMX, PMY, MMX, MMY, ProbeLength, MatchScore, Position, TopStrand, SeqNum;
  SEXP PMProbe = NULL;

  int *p_pmx = NULL;
  int *p_pmy = NULL;
  int *p_mmx = NULL;
  int *p_mmy = NULL;
  int *p_ProbeLength = NULL;
  double *p_MatchScore = NULL;
  int *p_Position =NULL;
  int *p_TopStrand =NULL;
  int *p_SeqNum =NULL;

  /*Allocating storage of R variables:*/

/*1. count the number of outputs:*/
int count=1; /*for the sequence number*/
if(isPM) count=count+2;
if(isMM) count = count+2;
if(isProbeLength) count++;
if(isPMProbeSeq) count++;
if(isMatchScore) count++;
if(isPosition) count++;
if(isTopStrand) count++;

/*Declaring the name of the vector*/
PROTECT(hit_names = allocVector(STRSXP,count));    /*Declaring the name of the vector*/
/* Creating a list with appropriate vector elements:  //added 1 for SeqNum */
PROTECT(hitList = allocVector(VECSXP, count));

int eleCount=0;

if(isPM)
{
  PROTECT(PMX = NEW_INTEGER(TOTALHITS));
  p_pmx = INTEGER_POINTER(PMX);
SET_STRING_ELT(hit_names, eleCount,mkChar("PMX"));/*Declaring the name of the vector*/
   SET_VECTOR_ELT(hitList, eleCount, PMX);
eleCount++;

  PROTECT(PMY = NEW_INTEGER(TOTALHITS));
  p_pmy = INTEGER_POINTER(PMY);
SET_STRING_ELT(hit_names,eleCount,mkChar("PMY"));/*Declaring the name of the vector*/
   SET_VECTOR_ELT(hitList, eleCount, PMY);
eleCount++;
}
if(isMM)
{
  PROTECT(MMX = NEW_INTEGER(TOTALHITS));
  p_mmx = INTEGER_POINTER(MMX);
SET_STRING_ELT(hit_names,eleCount,mkChar("MMX"));/*Declaring the name of the vector*/
   SET_VECTOR_ELT(hitList, eleCount, MMX);
 eleCount++;

PROTECT(MMY = NEW_INTEGER(TOTALHITS));
  p_mmy = INTEGER_POINTER(MMY);
SET_STRING_ELT(hit_names, eleCount,mkChar("MMY"));/*Declaring the name of the vector*/
   SET_VECTOR_ELT(hitList,  eleCount, MMY);
 eleCount++;
}
if(isProbeLength)
{
  PROTECT(ProbeLength = NEW_INTEGER(TOTALHITS));
  p_ProbeLength = INTEGER_POINTER(ProbeLength);

SET_STRING_ELT(hit_names, eleCount,mkChar("ProbeLength"));/*Declaring the name of the vector*/

   SET_VECTOR_ELT(hitList,  eleCount, ProbeLength);
 eleCount++;
}
if(isPMProbeSeq)
{

	PROTECT(PMProbe = allocVector(STRSXP,TOTALHITS));
SET_STRING_ELT(hit_names, eleCount,mkChar("PMProbe"));/*Declaring the name of the vector*/

   SET_VECTOR_ELT(hitList,  eleCount, PMProbe);
 eleCount++;
}
if(isMatchScore)
{
  PROTECT(MatchScore = NEW_NUMERIC(TOTALHITS));
  p_MatchScore = NUMERIC_POINTER(MatchScore);

SET_STRING_ELT(hit_names, eleCount,mkChar("MatchScore"));/*Declaring the name of the vector*/
   SET_VECTOR_ELT(hitList,  eleCount, MatchScore);
 eleCount++;
}
if(isPosition)
{
  PROTECT(Position = NEW_INTEGER(TOTALHITS));
  p_Position = INTEGER_POINTER(Position);

SET_STRING_ELT(hit_names, eleCount, mkChar("Position"));/*Declaring the name of the vector*/
   SET_VECTOR_ELT(hitList,  eleCount, Position);
 eleCount++;
}
if(isTopStrand)
{
  PROTECT(TopStrand = NEW_INTEGER(TOTALHITS));
  p_TopStrand = INTEGER_POINTER(TopStrand);

SET_STRING_ELT(hit_names, eleCount,mkChar("TopStrand"));/*Declaring the name of the vector*/
   SET_VECTOR_ELT(hitList,  eleCount, TopStrand);
 eleCount++;

}
/*added - to identify the sequence number*/
  PROTECT(SeqNum = NEW_INTEGER(TOTALHITS));
  p_SeqNum = INTEGER_POINTER(SeqNum);
	SET_STRING_ELT(hit_names,  eleCount, mkChar("SeqNum"));/*Declaring the name of the vector*/
  SET_VECTOR_ELT(hitList,  eleCount, SeqNum);
   eleCount++;

  /*next, store hitStorage for specific sequence of interest*/
  GDACSequenceHitItemType *hitStorage=NULL;
  hitStorage = new GDACSequenceHitItemType[TOTALHITS];




int curStartIndex=0; /*for the starting index of each sequence*/
numSeqToParse=0;
while(numSeqToParse<lengthList)
{
	int cur = rangeListArray[numSeqToParse];
	CGDACSequenceItem seq;
	data.GetSequenceItem(cur, seq);  /*get sequence information*/
	if(*isVerbose)
	{
		cout<<"Sequence " <<cur <<", " <<seq.FullName() <<" containing: " << seq.GetNumberHits() <<" number of hits" <<endl;
	}
 	int NUMREAD=seq.GetNumberHits();
	for (int i=0; i< NUMREAD; i++) {
	     seq.GetHitItem(i,hitStorage[curStartIndex+i], true);  /* Initialize all elements to zero.*/
	}
	/*Reading and storing information to the allocated memory of the R variables*/
	for(int pos=curStartIndex; pos< curStartIndex+NUMREAD; pos++)
	{  /*change to reflect assignment of each sequence*/
		if(isPM)
		{
		  p_pmx[pos]=hitStorage[pos].PMX;
	        		p_pmy[pos]=hitStorage[pos].PMY;
	        		}
	        if(isMM){
	          p_mmx[pos]=hitStorage[pos].MMX;
			        p_mmy[pos]=hitStorage[pos].MMY;
			        }
        	if(isProbeLength){	p_ProbeLength[pos]= (int) hitStorage[pos].ProbeLength;}
		if(isPMProbeSeq){	SET_STRING_ELT(PMProbe, pos ,  mkChar( hitStorage[pos].PMProbe.c_str()  )   );  }
	        if(isMatchScore){	p_MatchScore[pos] = hitStorage[pos].MatchScore;}
	        if(isPosition){	p_Position[pos]=hitStorage[pos].Position;}
	        if(isTopStrand){	p_TopStrand[pos]= (int) hitStorage[pos].TopStrand; }
		/*added - to identify the origin of the hit*/
	        p_SeqNum[pos]= (int) cur;
	  }

	/*Update curStartIndex to prepare for storing the next set of sequence*/
	curStartIndex = curStartIndex + NUMREAD;
	numSeqToParse++;
} /*closing of the loop for each sequence*/

   /*and attaching the vector names:*/
   setAttrib(hitList, R_NamesSymbol, hit_names);

   /*finally delete the temperary variables*/
   delete[] hitStorage;
   hitStorage=NULL;

   UNPROTECT(count+2);
   return hitList;
}




/*Return the header information of the BPMAP file*/
/*including the BPMAP version number and the number of sequences in the BPMAP file*/
SEXP readBPMAPFileHeader(SEXP fileName)
{
     CBPMAPFileData data;
     const char *fname = CHAR(STRING_ELT(fileName,0));
     data.SetFileName(fname);
     /* data.SetFileName( CHAR(STRING_ELT(fileName, 0)) );*/
     if( data.ReadHeader()==false){
       cout<<"Fail to read header\n" <<endl;
       return R_NilValue;
     }

     SEXP  myVersion, numSeq, list_names, list;

     int *p_numSeq;
     /*Allocating storage:*/
     PROTECT(numSeq = NEW_INTEGER(1));
     p_numSeq = INTEGER_POINTER(numSeq);
     p_numSeq[0]=data.GetNumberSequences(); /*get the number of sequences stored in the file*/


     double *p_myVersion;
     PROTECT(myVersion = NEW_NUMERIC(1));
     p_myVersion = NUMERIC_POINTER(myVersion);
     p_myVersion[0]=data.GetVersion();  /*get the version number of the file*/

     /*the header of the returning list*/
     PROTECT(list_names = allocVector(STRSXP,2));
     SET_STRING_ELT(list_names,0,mkChar("Version"));
     SET_STRING_ELT(list_names,1,mkChar("NumSeq"));

     /* Creating a list with 2 vector elements: version number and number of sequences*/
     PROTECT(list = allocVector(VECSXP, 2));
     SET_VECTOR_ELT(list, 0, myVersion);
     SET_VECTOR_ELT(list, 1, numSeq);
     /*      and attaching the vector names:*/
     setAttrib(list, R_NamesSymbol, list_names);

     UNPROTECT(4);
     return list; /*returning the list*/

}

/*
* Loading Sequence Information
* fileName: R string specifying the full path of the file
* seqToRead: the sequence to read: a number from 0 to (number of sequences in file -1)
* Returns Nothing.
*/
SEXP readBPMAPSeqHeader(SEXP fileName, SEXP seqToRead)
{

  int cur = INTEGER_VALUE(seqToRead); /*the sequence to load*/
  CBPMAPFileData data;
  data.SetFileName( CHAR(STRING_ELT(fileName, 0)) );
  if(data.Read()==false)  {  /*reading the file information*/
     cout<<"Fail to read\n"  <<endl;
     return R_NilValue;
  }

   cout<< "Reading Sequence Information from " <<data.GetFileName() <<endl;
   /*sequence's description            */
   CGDACSequenceItem seq; /*data structure for the sequence*/
   data.GetSequenceItem(cur, seq);  /*get sequence information*/
   cout<<seq.FullName() <<" with " <<endl;
   cout<<"Containing: " << seq.GetNumberHits() <<" number of hits" <<endl;

   return R_NilValue;
}

/*New function
Output each sequence name and number of hit for the sequence into two vectors
*/
SEXP readBPMAPAllSeqHeader(SEXP fileName)
{

  CBPMAPFileData data;
  data.SetFileName( CHAR(STRING_ELT(fileName, 0)) );

/*Obtain the number of sequences to read*/
     if( data.ReadHeader()==false){
       cout<<"Fail to read header\n" <<endl;
       return R_NilValue;
     }
     int TOTALNUMSEQ=data.GetNumberSequences(); /*get the number of sequences stored in the file*/


  if(data.Read()==false)  {  /*reading the file information*/
     cout<<"Fail to read\n"  <<endl;
     return R_NilValue;
  }

   cout<< "Reading Sequence Information from " <<data.GetFileName() <<endl;
   /*sequence's description*/

     SEXP  seqName, groupName, version, probeMapping, seqNumber, seqNumHits, list_names, list;

     /*Allocating storage:*/
     PROTECT(seqName = allocVector(STRSXP, TOTALNUMSEQ));
     PROTECT(groupName = allocVector(STRSXP, TOTALNUMSEQ));
     PROTECT(version = allocVector(STRSXP, TOTALNUMSEQ));

     int *p_seqNumHits;
     PROTECT(seqNumHits = NEW_INTEGER(TOTALNUMSEQ));
     p_seqNumHits = INTEGER_POINTER(seqNumHits);

     int *p_probeMapping;
     PROTECT(probeMapping = NEW_INTEGER(TOTALNUMSEQ));
     p_probeMapping = INTEGER_POINTER(probeMapping);

     int *p_seqNumber;
     PROTECT(seqNumber = NEW_INTEGER(TOTALNUMSEQ));
     p_seqNumber = INTEGER_POINTER(seqNumber);


for(int cur=0; cur<TOTALNUMSEQ ; cur++){
   CGDACSequenceItem seq; /*data structure for the sequence*/
   data.GetSequenceItem(cur, seq);  /*get sequence information*/
/*   cout<<seq.FullName() <<" with " <<endl;*/
   /*cout<<"Containing: " << seq.GetNumberHits() <<" number of hits" <<endl;*/
   SET_STRING_ELT(seqName, cur , mkChar(seq.GetName().c_str()));
   SET_STRING_ELT(groupName, cur , mkChar(seq.GroupName().c_str()));
   SET_STRING_ELT(version, cur , mkChar(seq.GetSeqVersion().c_str()));
   p_seqNumHits[cur]=seq.GetNumberHits();
   p_probeMapping[cur]=seq.GetProbeMapping();
   p_seqNumber[cur]=seq.GetNumber();

   /*
   cout<<"Sequence name: " <<seq.GetName() <<", Group name: " <<seq.GroupName() <<", Seq Version: " << seq.GetSeqVersion()
   <<", Probe Mapping: " <<seq.GetProbeMapping() <<", Get Number: " <<seq.GetNumber() <<", Get NumberHits: " <<seq.GetNumberHits()
   <<", Num Parameters: " <<seq.GetNumberParameters() <<endl;
   */
}



     /*the header of the returning list*/
     PROTECT(list_names = allocVector(STRSXP,6));
     SET_STRING_ELT(list_names,0,mkChar("SeqName"));
     SET_STRING_ELT(list_names,1,mkChar("GroupName"));
     SET_STRING_ELT(list_names,2,mkChar("version"));
     SET_STRING_ELT(list_names,3,mkChar("probeMapping"));
     SET_STRING_ELT(list_names,4,mkChar("seqNum"));
     SET_STRING_ELT(list_names,5,mkChar("NumHits"));

     /* Creating a list with 2 vector elements: version number and number of sequences*/
     PROTECT(list = allocVector(VECSXP, 6));
     SET_VECTOR_ELT(list, 0, seqName);
     SET_VECTOR_ELT(list, 1, groupName);
     SET_VECTOR_ELT(list, 2, version);
     SET_VECTOR_ELT(list, 3, probeMapping);
     SET_VECTOR_ELT(list, 4, seqNumber);
     SET_VECTOR_ELT(list, 5, seqNumHits);

     /*      and attaching the vector names:*/
     setAttrib(list, R_NamesSymbol, list_names);

     UNPROTECT(8);
     return list; /*returning the list*/

}




/*
 *Merging BPMAP data frame and CEL data frame.
 *Each data frame must contain columns X and Y to specify probe positions on the array.
 *returns named list of merged columns
*/

SEXP BPMAPCelMerger(SEXP list, SEXP celList)
{

  int nProtect=0;
  SEXP BPMAPx = NULL, BPMAPy=NULL, names = getAttrib(list, R_NamesSymbol);
  for (int i=0; i<length(list); i++)
    {
      if(strcmp(CHAR(STRING_ELT(names, i)), "X") == 0)
      {
	BPMAPx=VECTOR_ELT(list,i);
      }

      else if(strcmp(CHAR(STRING_ELT(names, i)), "Y") == 0)
      {
	BPMAPy=VECTOR_ELT(list,i);
      }
   }

  if ((BPMAPx==NULL)||(BPMAPy==NULL)){
    cout << "BPMAP file does not contain variable X or Y. Read incorrectly" <<endl;
    return R_NilValue;
  }

/*########################### Declare output intensities ###############################*/
  SEXP *inten = new SEXP[length(celList)-2];
  SEXP *intenOut=new SEXP[length(celList)-2];  /*contain X,Y, and all the intensities*/
  double **p_intenOut;
  p_intenOut = new double*[length(celList)-2];

  for(int i=0; i< length(celList)-2; i++){
    PROTECT(intenOut[i] = NEW_NUMERIC(length(VECTOR_ELT(list,0))));
    nProtect++;
    p_intenOut[i] = NUMERIC_POINTER(intenOut[i]);
  }
/*########################################################################################*/


  SEXP x=NULL, y=NULL, celnames = getAttrib(celList, R_NamesSymbol);
  int intenCount=0;
  for (int i=0; i<length(celList); i++){
    if(strcmp(CHAR(STRING_ELT(celnames, i)), "X") == 0)
    {
	x=VECTOR_ELT(celList,i);
    }
    else if(strcmp(CHAR(STRING_ELT(celnames, i)), "Y") == 0)
    {
	y=VECTOR_ELT(celList,i);
    }

    else{
	/*cout<<"intensity in column " <<i <<endl;*/
	/*Save the variables out*/
	inten[intenCount] = VECTOR_ELT(celList,i);
	intenCount++;
	}
  }


  if ((x==NULL)||(y==NULL))
  {
	cout << "Cel file does not contain variable X or Y. Read incorrectly" <<endl;
	return R_NilValue;
  }

  int nrecord = length(x);
  int nBPMAP = length(BPMAPx);
  int curA=0, curB=0;
  int stop=0;
  while(curA<nBPMAP &&curB<nrecord)
  {

    if (INTEGER(BPMAPy)[curA]==INTEGER(y)[curB])
    {
      if (INTEGER(BPMAPx)[curA]==INTEGER(x)[curB])
      {
	stop++;
	for (int col=0; col<length(celList) -2; col++)
	{
	  p_intenOut[col][curA]=REAL(inten[col])[curB];
	}
	curA++;
      }
      else if(INTEGER(BPMAPx)[curA] < INTEGER(x)[curB])
      {
	cout<<"LEFT OVER READ... ERROR" <<endl;
	return R_NilValue;
	}
      else if(INTEGER(BPMAPx)[curA] > INTEGER(x)[curB])
      {
	curB++;
      }
    }

    else if (INTEGER(BPMAPy)[curA]>INTEGER(y)[curB])
    {
        curB++;
    }
    else
    {
	cout<<"FAILED.. lists not sorted"<<endl;
	cout<<curA <<", " <<curB <<endl;
	return R_NilValue;
    }
  }


  /*Now creating a list to output all of them*/
  /*############################### creating list names ####################*/
  SEXP hitList_names;
  PROTECT(hitList_names = allocVector(STRSXP, length(list)+length(celList)-2));
  nProtect++;

  /*############################### creating a list to output all of them #######################*/
  SEXP hitList;
  PROTECT(hitList = allocVector(VECSXP, length(list)-2 + length(celList)));
  nProtect++;

  /* attaching BPMAP vectors to list:*/
  for (int i=0; i<length(list); i++)
  {
    /*Save the variables out*/
    SET_VECTOR_ELT(hitList, i, VECTOR_ELT(list,i));
    SET_STRING_ELT(hitList_names, i, mkChar(CHAR(STRING_ELT(names, i))));
  }

  /*assume we have X and Y as the first and second columns of intensities*/
  int excludeXYcount=0;  /*because X and Y in the combined cel files might not happen to be on a specific column, this variable is needed*/
  for(int i=0; i< length(celList); i++)
  {

    if(!((strcmp(CHAR(STRING_ELT(celnames, i)), "X") == 0)||(strcmp(CHAR(STRING_ELT(celnames, i)), "Y") == 0)))
    {
    SET_VECTOR_ELT(hitList, excludeXYcount+length(list), intenOut[excludeXYcount]); /*append to hitList*/
    SET_STRING_ELT(hitList_names, length(list)+excludeXYcount, mkChar(CHAR(STRING_ELT(celnames, i))));
    excludeXYcount++;
    }

  }
  setAttrib(hitList, R_NamesSymbol, hitList_names);

  UNPROTECT(nProtect);
  /*cout<<stop <<endl;*/
  return hitList;
}


/*
 *Convert a vector of sequences of same length to an numerical matrix
 */

SEXP matrixSeq(SEXP seq, SEXP sizeSeq)
{
 int size = INTEGER_VALUE(sizeSeq);
 SEXP mat;
 int seqNum = LENGTH(seq);
 char *tmp=new char[size];

 PROTECT(seq = AS_CHARACTER(seq));
 PROTECT(mat = allocMatrix(INTSXP, seqNum, size));

 for(int i=0; i< seqNum; i++)
 {
   for(int j=0; j<size; j++)
   {
     strcpy(tmp,CHAR(STRING_ELT(seq,i)));
     if( tmp[j]=='A')
     {
       INTEGER(mat)[j*seqNum + i] = 1;   /* not j+size*i*/
     }
     else if (tmp[j]=='G')
     {
       INTEGER(mat)[j*seqNum + i] = 2;
     }
     else if (tmp[j]=='C')
     {
       INTEGER(mat)[j*seqNum + i] = 3;
     }
     else
     {
       INTEGER(mat)[j*seqNum + i] = 4;
     }
   }
 }
 UNPROTECT(2);
 return mat;
}

/*
* Convert a vector of sequences into a numerical matrix of adjacent pair codes of sequences.
*
*/
SEXP seqPair(SEXP seq, SEXP sizeSeq)
{
  int size = INTEGER_VALUE(sizeSeq);
  SEXP mat;
  int seqNum = LENGTH(seq);
  char *tmp=new char[size];

  PROTECT(seq = AS_CHARACTER(seq));
  PROTECT(mat = allocMatrix(INTSXP, seqNum, size));

  for(int i=0; i< seqNum; i++){
    for(int j=0; j<size; j++){
      strcpy(tmp,CHAR(STRING_ELT(seq,i)));
      if(( tmp[j]=='A')&&( tmp[j+1]=='A'))
      {
        INTEGER(mat)[j*seqNum + i] = 11;   /* not j+size*i*/
      }
      else if (( tmp[j]=='A')&&(tmp[j+1]=='G'))
      {
        INTEGER(mat)[j*seqNum + i] = 12;
      }
      else if (( tmp[j]=='A')&&(tmp[j+1]=='C'))
      {
        INTEGER(mat)[j*seqNum + i] = 13;
      }
      else if (( tmp[j]=='A')&&(tmp[j+1]=='T'))
      {
        INTEGER(mat)[j*seqNum + i] = 14;
      }
      else if (( tmp[j]=='G')&&(tmp[j+1]=='A'))
      {
        INTEGER(mat)[j*seqNum + i] = 21;
      }
      else if (( tmp[j]=='G')&&(tmp[j+1]=='G'))
      {
        INTEGER(mat)[j*seqNum + i] = 22;
      }
      else if (( tmp[j]=='G')&&(tmp[j+1]=='C'))
      {
        INTEGER(mat)[j*seqNum + i] = 23;
      }
      else if (( tmp[j]=='G')&&(tmp[j+1]=='T'))
      {
        INTEGER(mat)[j*seqNum + i] = 24;
      }
      else if (( tmp[j]=='C')&&(tmp[j+1]=='A'))
      {
        INTEGER(mat)[j*seqNum + i] = 31;
      }
      else if (( tmp[j]=='C')&&(tmp[j+1]=='G'))
      {
        INTEGER(mat)[j*seqNum + i] = 32;
      }
      else if (( tmp[j]=='C')&&(tmp[j+1]=='C'))
      {
        INTEGER(mat)[j*seqNum + i] = 33;
      }
      else if (( tmp[j]=='C')&&(tmp[j+1]=='T'))
      {
        INTEGER(mat)[j*seqNum + i] = 34;
      }
      else if (( tmp[j]=='T')&&(tmp[j+1]=='A'))
      {
        INTEGER(mat)[j*seqNum + i] = 41;
      }
      else if (( tmp[j]=='T')&&(tmp[j+1]=='G'))
      {
        INTEGER(mat)[j*seqNum + i] = 42;
      }
      else if (( tmp[j]=='T')&&(tmp[j+1]=='C'))
      {
        INTEGER(mat)[j*seqNum + i] = 43;
      }
      else if (( tmp[j]=='T')&&(tmp[j+1]=='T'))
      {
        INTEGER(mat)[j*seqNum + i] = 44;
      }
      else
      {
        cout << "ERROR: make sure the probes has been converted to characters using <as.character()>" <<endl;
        break;
      }

    }
  }

  UNPROTECT(2);
  return mat;
}



} /** end extern "C" **/

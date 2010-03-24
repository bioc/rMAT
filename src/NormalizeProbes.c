#include <R.h>
#include <Rmath.h>
#include <float.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_sort_vector.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if defined(__APPLE__) || defined(macintosh)
  #ifdef __MAC_OS_X_VERSION_MIN_REQUIRED
    #if __MAC_OS_X_VERSION_MIN_REQUIRED >= 1060
      #include <dispatch/dispatch.h>
      #define HAVE_DISPATCH 1
    #endif
  #endif
#endif

int convertSeq(char seq);
int convertNum(char seq1, char seq2);

void createDesignMatrixMAT(gsl_matrix *seqNumCount, gsl_vector *copyNumber, gsl_matrix *X, char** seq);
void createDesignMatrixMATRow(gsl_matrix *seqNumCount, gsl_vector *copyNumber, gsl_vector *XRow, int i, char** seq);

void createSeqMatrixCount(gsl_matrix *seqNumCount, char **seq);
void createDesignMatrixPairBinned(gsl_matrix *seqNumCount, 
  gsl_matrix *pairNumCount1,
  gsl_matrix *pairNumCount2,
  gsl_matrix *pairNumCount3,
  gsl_matrix *pairNumCount4,
  gsl_vector *copyNumber, 
  gsl_matrix *X);

void createDesignMatrixPairBinnedRow(gsl_matrix *seqNumCount, 
  gsl_matrix *pairNumCount1,
  gsl_matrix *pairNumCount2,
  gsl_matrix *pairNumCount3,
  gsl_matrix *pairNumCount4,
  gsl_vector *copyNumber, 
  gsl_vector *XRow, int i);

void createPairMatrixCount(gsl_matrix *pairNumCount1, 
  gsl_matrix *pairNumCount2, 
  gsl_matrix *pairNumCount3,
  gsl_matrix *pairNumCount4,
  char** seq);

void createDesignMatrixPair( gsl_matrix *seqNumCount, gsl_matrix *pairNumMatrix, gsl_vector *copyNumber, gsl_matrix *X, char **seq);
void createDesignMatrixPairRow( gsl_matrix *seqNumCount, gsl_matrix *pairNumMatrix, gsl_vector *copyNumber, gsl_vector *XRow, int i, char **seq);

void MATScore(double *C, double *I, int *nProbes, int *nArraysC, int *nArraysI, int *position, double *dMax, double *MATScores, int *seqNum);
void MATNullDistribution(int *position, int *nProbes, double *dMax, double *MATScore, double *sigma0, double *mu0, int  *seqNum);
void MATpValue(int nProbes, double *MATScore, double sigma0, double mu0, double *pValues);
double MATcutoffFDR(int *position, int nProbes, double dMax, double *MATScores, double mu0, double FDR, int *regions, int *seqNum);
int mergeMATScores(int *position, int nProbes, double dMerge, double *MATScores, double m0, double cutoff, int sign, int *regions, int *Mum);
void MAT(double *C, double *I, int *nProbes, int *nArraysC, int *nArraysI, int *position, double *dMax, double *dMerge, double *threshold, double *MATScores, double *pValues, int *method, int *regions, int *verbose, int *seqNum, int *numRegions);
void NormalizeProbes(char **seq, double *y, double *yNormalized, int *nProbes, int *nArrays, double *copyNumber, int *method, int *robust, double *adjRSquare, double *RSquare, double *BIC, double *outBeta, int *betaLength, int *all,  int *MATScaling,int *isVerbose);
void getIndices(int *regions, int *nProbes, int *numRegions, int *StartRegion, int *EndRegion);
void callEnrichedRegions(double *MATScores, int *nProbes, int *position, double *dMerge, double *dMax, double *threshold, double *pValues, int *method, int *regions, int *verbose, int *seqNum, int *numRegions);
void normArray(char **seq, double *y, double *yNormalized, int *nProbes, int *nArrays, double *copyNumber, int *method, int *robust, double *adjRSquare, double *RSquare, double *BIC, double *outBeta, int *betaLength, int *all,  int *MATScaling,int *isVerbose, 
               gsl_matrix *pairNumCount1, gsl_matrix *pairNumCount2, gsl_matrix *pairNumCount3, gsl_matrix *pairNumCount4,
               gsl_matrix *seqNumCount,
               int nProbesTotal,
               int nVariables,
               int nBins,
               int nProbesPerBin,
               gsl_vector_view yVector, gsl_vector_view copyNumberVector,
               int j);

/** Main method to normalize probes**/
      void NormalizeProbes(char **seq, double *y, double *yNormalized, int *nProbes, int *nArrays, double *copyNumber, int *method, int *robust, double *adjRSquare, double *RSquare, double *BIC, double *outBeta, int *betaLength, int *all,  int *MATScaling,int *isVerbose)
      {

      /** Declaring gsl variables **/
        gsl_vector_view yVector,copyNumberVector;
        // gsl_matrix_view pairNumMatrix,seqNumMatrix;
        gsl_matrix *seqNumCount;
        gsl_matrix *pairNumCount1=NULL, *pairNumCount2=NULL, *pairNumCount3=NULL, *pairNumCount4=NULL;

        int nProbesTotal=*nProbes;
        int j=0,nVariables;
        #ifdef HAVE_DISPATCH
          dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
        #endif

        /**if not using all probes to compute (for faster computation), then use the minimum of 300,000 or number of probes, which ever is less **/
        if((*all)==0)
        {
          *nProbes=GSL_MIN_INT(300000,*nProbes);
        }

        /** Number of bins used in the MAT normalization **/
        int nBins=100;
        /** Average number of probes per bins **/
        int nProbesPerBin=(nProbesTotal)/nBins;

        copyNumberVector=gsl_vector_view_array(copyNumber,nProbesTotal);
        seqNumCount=gsl_matrix_calloc(nProbesTotal,4);

        /*****************************************************************************************************/
        if(*isVerbose)
        {
          printf("** Create count matrices **\n");
        }
        createSeqMatrixCount(seqNumCount, seq);
        /*****************************************************************************************************/

        if(*method==1) /** This is the MAT model **/
        {
          nVariables=1+3*25+4+1;
        }
        else  /** This is the pair binned model --into 4 bins**/
        {
          nVariables=73;
          pairNumCount1=gsl_matrix_calloc(nProbesTotal,16);
          pairNumCount2=gsl_matrix_calloc(nProbesTotal,16);
          pairNumCount3=gsl_matrix_calloc(nProbesTotal,16);
          pairNumCount4=gsl_matrix_calloc(nProbesTotal,16);
          createPairMatrixCount(
            pairNumCount1, 
            pairNumCount2, 
            pairNumCount3,
            pairNumCount4,
            seq);
        }
        #ifdef HAVE_DISPATCH
          dispatch_apply(*nArrays, queue, ^(size_t jj) {normArray(seq, y, yNormalized, nProbes, nArrays, copyNumber, method, robust, adjRSquare, RSquare, BIC, outBeta, betaLength, all,  MATScaling, isVerbose, 
            pairNumCount1, pairNumCount2, pairNumCount3, pairNumCount4,
            seqNumCount,
            nProbesTotal,
            nVariables,
            nBins,
            nProbesPerBin,
            yVector, copyNumberVector,
            jj);
          });
        #else
          for(j=0;j<*nArrays;j++)
            normArray(seq, y, yNormalized, nProbes, nArrays, copyNumber, method, robust, adjRSquare, RSquare, BIC, outBeta, betaLength, all,  MATScaling, isVerbose, 
            pairNumCount1, pairNumCount2, pairNumCount3, pairNumCount4,
            seqNumCount,
            nProbesTotal,
            nVariables,
            nBins,
            nProbesPerBin,
            yVector, copyNumberVector,
            j);
        #endif

        if(*method==2)
        {
          gsl_matrix_free(pairNumCount1);
          gsl_matrix_free(pairNumCount2);
          gsl_matrix_free(pairNumCount3);
          gsl_matrix_free(pairNumCount4);    
        }

        gsl_matrix_free(seqNumCount);  

        if(*isVerbose)
        {
          printf("** End of NormalizeProbes procedure **\n");
        }
      }
      
/*\Assume sequence size = 25*/
void createDesignMatrixMAT(gsl_matrix *seqNumCount, gsl_vector *copyNumber, gsl_matrix *X, char **seq)
{

  int nProbes=X->size1;
  int i,j;

  for(i=0;i<nProbes;i++)
  {
  /** Intercept **/
    gsl_matrix_set(X,i,0,1);
  /** Nucleotide positional effects **/
    int tempInt;
    for(j=0;j<25;j++)	
    {
      tempInt = convertSeq(seq[i][j]);
      if(tempInt !=4)
      {
        gsl_matrix_set(X,i,1+3*j+tempInt-1,1);
      }
    }
    /** Nucleotide counts effects **/
    for(j=0;j<4;j++)	
      gsl_matrix_set(X,i,1+3*25+j,gsl_pow_2(gsl_matrix_get(seqNumCount,i,j)));      
    /** Copy number **/
    gsl_matrix_set(X,i,1+3*25+4,log(gsl_vector_get(copyNumber,i)));   

  }
}


void createDesignMatrixMATRow(gsl_matrix *seqNumCount, gsl_vector *copyNumber, gsl_vector *XRow, int i, char** seq)
{
  // int nVariables=XRow->size;
  int j;
  gsl_vector_set_zero(XRow);
  /** Intercept **/
  gsl_vector_set(XRow,0,1);
  /** Nucleotide positional effects **/
  int tempInt;
  for(j=0;j<25;j++)	
  {
    tempInt = convertSeq(seq[i][j]);
    if(tempInt!=4)
      gsl_vector_set(XRow,1+3*j+tempInt-1,1);
  }
  /** Nucleotide counts effects **/
  for(j=0;j<4;j++)	
    gsl_vector_set(XRow,1+3*25+j,gsl_pow_2(gsl_matrix_get(seqNumCount,i,j)));      
  /** Copy number **/
  gsl_vector_set(XRow,1+3*25+4,log(gsl_vector_get(copyNumber,i))); 

}

int convertSeq(char seq)
{
  if(seq == 'A')	
    return 1;
  else if(seq == 'G')
    return 2;
  else if(seq == 'C')
    return 3;
  else if(seq == 'T')
    return 4;
  else{
    printf("Error.. The base is not A,G, C, or T\n");
    return(0);
  }
}

int convertNum(char seq1, char seq2)
{
  int num= 0;
  if(( seq1=='A')&&( seq2=='A'))
  {
    num = 1;   /*not j+size*i*/
  }
  else if (( seq1=='A')&&(seq2=='G')){
    num = 2;
  }
  else if (( seq1=='A')&&(seq2=='C')){
    num = 3;
  }
  else if (( seq1=='A')&&(seq2=='T')){
    num = 4;
  }
  else if (( seq1=='G')&&(seq2=='A')){
    num = 5;
  }
  else if (( seq1=='G')&&(seq2=='G')){
    num = 6;
  }
  else if (( seq1=='G')&&(seq2=='C')){
    num = 7;
  }
  else if (( seq1=='G')&&(seq2=='T')){
    num = 8;
  }
  else if (( seq1=='C')&&(seq2=='A')){
    num = 9;
  }
  else if (( seq1=='C')&&(seq2=='G')){
    num = 10;
  }
  else if (( seq1=='C')&&(seq2=='C')){
    num = 11;
  }
  else if (( seq1=='C')&&(seq2=='T')){
    num = 12;
  }
  else if (( seq1=='T')&&(seq2=='A')){
    num = 13;
  }
  else if (( seq1=='T')&&(seq2=='G')){
    num = 14;
  }
  else if (( seq1=='T')&&(seq2=='C')){
    num = 15;
  }
  else if (( seq1=='T')&&(seq2=='T')){
    num = 16;
  }
  else{
    printf("error\n");
  }
  return num;   
}

void createDesignMatrixPair( gsl_matrix *seqNumCount, gsl_matrix *pairNumMatrix, gsl_vector *copyNumber, gsl_matrix *X, char** seq)
{

  int nProbes=X->size1;
  // int nVariables=X->size2;
  int i,j;
  int tempInt;
  for(i=0;i<nProbes;i++)
  {
  /** Intercept **/
    gsl_matrix_set(X,i,0,1);
  /** Nucleotide positional pair effects **/
    tempInt = convertNum(seq[i][0], seq[i][1]);
    if(tempInt !=16)
    {
      gsl_matrix_set(X,i,1+tempInt-1,1);
    }
    /*      if(gsl_matrix_get(pairNumMatrix,i,0)!=16)*/
    /*	gsl_matrix_set(X,i,1+gsl_matrix_get(pairNumMatrix,i,0)-1,1);*/
    for(j=1;j<24;j++)
    {
      tempInt = convertNum(seq[i][j], seq[i][j+1]);	
      if(tempInt <=12)
        gsl_matrix_set(X,i,1+15+12*(j-1)+tempInt-1,1);
    /*	if(gsl_matrix_get(pairNumMatrix,i,j)<=12)*/
    /*	  gsl_matrix_set(X,i,1+15+12*(j-1)+gsl_matrix_get(pairNumMatrix,i,j)-1,1);*/
    /** Nucleotide counts effects **/
    }
    for(j=0;j<4;j++)	
    {
      gsl_matrix_set(X,i,1+15+12*23+j,gsl_pow_2(gsl_matrix_get(seqNumCount,i,j)));
      gsl_matrix_set(X,i,1+15+12*23+4+j,gsl_pow_3(gsl_matrix_get(seqNumCount,i,j)));
    }
    /** Copy number **/
    gsl_matrix_set(X,i,1+15+12*23+8,gsl_vector_get(copyNumber,i));   
  }
}

void createDesignMatrixPairRow(gsl_matrix *seqNumCount, gsl_matrix *pairNumMatrix, gsl_vector *copyNumber, gsl_vector *XRow, int i, char **seq)
{
  // int nVariables=XRow->size;
  int j;
  int tempInt;

/** Intercept **/
  gsl_vector_set(XRow,0,1);
/** Nucleotide positional pair effects **/
  tempInt = convertNum(seq[i][0], seq[i][1]);
  if(tempInt !=16)
  {
    gsl_vector_set(XRow,1+tempInt-1,1);
  }
  for(j=1;j<24;j++)	
  {
    tempInt = convertNum(seq[i][j], seq[i][j+1]);	
    if(tempInt <=12)
      gsl_vector_set(XRow,1+15+12*(j-1)+tempInt-1,1);
  }
/** Nucleotide counts effects **/
  for(j=0;j<4;j++)	
  {
    gsl_vector_set(XRow,1+15+12*23+j,gsl_pow_2(gsl_matrix_get(seqNumCount,i,j)));
    gsl_vector_set(XRow,1+15+12*23+4+j,gsl_pow_3(gsl_matrix_get(seqNumCount,i,j)));
  }
/** Copy number **/
  gsl_vector_set(XRow,1+15+12*23+8,gsl_vector_get(copyNumber,i));   
}

void createDesignMatrixPairBinned(gsl_matrix *seqNumCount, 
  gsl_matrix *pairNumCount1,
  gsl_matrix *pairNumCount2,
  gsl_matrix *pairNumCount3,
  gsl_matrix *pairNumCount4,
  gsl_vector *copyNumber, 
  gsl_matrix *X)
{

  int nProbes=X->size1;
  // int nVariables=X->size2;
  int i,j;

  for(i=0;i<nProbes;i++)
  {
  /** Intercept **/
    gsl_matrix_set(X,i,0,1);
  /** Nucleotide positional pair effects **/
    for(j=0;j<15;j++)
    {
      gsl_matrix_set(X,i,1+j,gsl_matrix_get(pairNumCount1,i,j));
      gsl_matrix_set(X,i,1+15+j,gsl_matrix_get(pairNumCount2,i,j));
      gsl_matrix_set(X,i,1+30+j,gsl_matrix_get(pairNumCount3,i,j));
      gsl_matrix_set(X,i,1+45+j,gsl_matrix_get(pairNumCount4,i,j));	  	  
    }
  /** Nucleotide counts effects **/
    for(j=0;j<3;j++)	
      gsl_matrix_set(X,i,1+60+j,gsl_matrix_get(seqNumCount,i,j));
    for(j=0;j<4;j++)	
    {
      gsl_matrix_set(X,i,1+63+j,gsl_pow_2(gsl_matrix_get(seqNumCount,i,j)));
      gsl_matrix_set(X,i,1+67+j,gsl_pow_3(gsl_matrix_get(seqNumCount,i,j)));
    }     
  /** Copy number **/
    gsl_matrix_set(X,i,1+71,gsl_vector_get(copyNumber,i));   
  }
}

void createDesignMatrixPairBinnedRow(gsl_matrix *seqNumCount, 
  gsl_matrix *pairNumCount1,
  gsl_matrix *pairNumCount2,
  gsl_matrix *pairNumCount3,
  gsl_matrix *pairNumCount4,
  gsl_vector *copyNumber, 
  gsl_vector *XRow, int i)
{

  // int nVariables=XRow->size;
  int j;

/** Intercept **/
  gsl_vector_set(XRow,0,1);
/** Nucleotide positional pair effects **/
  for(j=0;j<15;j++)
  {
    gsl_vector_set(XRow,1+j,gsl_matrix_get(pairNumCount1,i,j));
    gsl_vector_set(XRow,1+15+j,gsl_matrix_get(pairNumCount2,i,j));
    gsl_vector_set(XRow,1+30+j,gsl_matrix_get(pairNumCount3,i,j));
    gsl_vector_set(XRow,1+45+j,gsl_matrix_get(pairNumCount4,i,j));	  	  
  }
/** Nucleotide counts effects **/
  for(j=0;j<3;j++)	
    gsl_vector_set(XRow,1+60+j,gsl_matrix_get(seqNumCount,i,j));
  for(j=0;j<4;j++)	
  {
    gsl_vector_set(XRow,1+63+j,gsl_pow_2(gsl_matrix_get(seqNumCount,i,j)));
    gsl_vector_set(XRow,1+67+j,gsl_pow_3(gsl_matrix_get(seqNumCount,i,j)));
  }     
/** Copy number **/
  gsl_vector_set(XRow,1+71,gsl_vector_get(copyNumber,i));   
}

void createSeqMatrixCount(gsl_matrix *seqNumCount, char **seq)
{

  int nProbes=seqNumCount->size1;
  int nNucleotides=seqNumCount->size2;
  int i=0,j=0;
  int tempInt;

  for(i=0;i<nProbes;i++)
  {
    for(j=0;j<nNucleotides;j++)
    {
      tempInt = convertSeq(seq[i][j]);
      /**printf("tempInt is %d", tempInt);**/
      gsl_matrix_set(seqNumCount,i, tempInt-1,gsl_matrix_get(seqNumCount,i,tempInt-1)+1);
    }
  }


}

/**for 25mer - can be modified for a general case**/
void createPairMatrixCount(gsl_matrix *pairNumCount1, 
  gsl_matrix *pairNumCount2, 
  gsl_matrix *pairNumCount3,
  gsl_matrix *pairNumCount4,
  char **seq)
{
  /** int nProbes=pairNum->size1;**/
  /**  int nNucleotides=pairNum->size2; **/
  int nProbes=pairNumCount1->size1;
  /**  int nNucleotides=(pairNumCount2->size2)*4; **/
  // int nNucleotides=(pairNumCount2->size2);
  int i=0,j=0;
  int tempInt;
  for(i=0;i<nProbes;i++)
  {
    for(j=0;j< 6;j++)
    {
      tempInt = convertNum(seq[i][j], seq[i][j+1]);
      gsl_matrix_set(pairNumCount1,i,tempInt-1,gsl_matrix_get(pairNumCount1,i,tempInt-1)+1);
      tempInt = convertNum(seq[i][j+6], seq[i][j+6+1]);
      gsl_matrix_set(pairNumCount2,i,tempInt-1,gsl_matrix_get(pairNumCount2,i,tempInt-1)+1);
      tempInt = convertNum(seq[i][j+12], seq[i][j+12+1]);
      gsl_matrix_set(pairNumCount3,i,tempInt-1,gsl_matrix_get(pairNumCount3,i,tempInt-1)+1);
      tempInt = convertNum(seq[i][j+18], seq[i][j+18+1]);
      gsl_matrix_set(pairNumCount4,i,tempInt-1,gsl_matrix_get(pairNumCount4,i,tempInt-1)+1);	
    }
  }
}


void callEnrichedRegions(double *MATScores, int *nProbes, int *position, double *dMerge, double *dMax, double *threshold, double *pValues, int *method, int *regions, int *verbose, int *seqNum, int *numRegions)
{

  double sigma0=0, mu0=0, cutoff=0;
  // double totalSeqsPosition=0;

  /** Compute the associated pValues **/
  MATNullDistribution(position, nProbes, dMax, MATScores, &sigma0, &mu0, seqNum);
  
  /** Compute the FDR cutoff **/
  if(*method==1) /** Based on MAT scores **/
  {
    if(*verbose)
    {
      printf("** Merging regions **\n");  
    }        
    numRegions[0] = mergeMATScores(position, *nProbes, *dMerge, MATScores, mu0, *threshold, 1, regions, seqNum);
  }
  else if(*method==2)
  {
    if(*verbose)
    {
      printf("** Calculating P-values **\n");  
    }    
    MATpValue(*nProbes, MATScores, sigma0, mu0, pValues);
    if(*verbose)
    {
      printf("** Merging regions **\n");  
    }        
    numRegions[0] = mergeMATScores(position, *nProbes, *dMerge, pValues, 0, *threshold, -1, regions, seqNum);
  }
  else if(*method==3)
  {
    if(*verbose)
    {
      printf("** Calculating FDR **\n");  
    }    
    
    cutoff=MATcutoffFDR(position, *nProbes, *dMerge, MATScores, mu0, *threshold, regions, seqNum);
    if(*verbose)
    {
      printf("** Merging regions **\n");  
    }    
    /** Form the regions using the threshold found **/
    numRegions[0] = mergeMATScores(position, *nProbes, *dMerge, MATScores, mu0, cutoff, 1, regions, seqNum);
  }
}

void getIndices(int *regions, int *nProbes, int *numRegions, int *StartRegion, int *EndRegion)
{
  int i=0;
  int count=0;
  
  for(i=1;i<=*numRegions;i++)
  {
    while(((regions[count]<i) | (regions[count]==0)) & (count<*nProbes))
    {
      count++;
    }
    StartRegion[i-1]=count+1;
    while((regions[count]==i) & (count<*nProbes))
    {
      count++;
    }
    EndRegion[i-1]=count;
  }
}


/*\ Assume data already sorted by sequence number, then by position */
void MAT(double *C, double *I, int *nProbes, int *nArraysC, int *nArraysI, int *position, double *dMax, double *dMerge, double *threshold, double *MATScores, double *pValues, int *method, int *regions, int *verbose, int *seqNum, int *numRegions)
{

  double sigma0=0, mu0=0, cutoff=0;
  // double totalSeqsPosition=0;

  MATScore(C, I, nProbes, nArraysC, nArraysI, position, dMax, MATScores, seqNum);
  if(*verbose)
  {
    printf("** Estimate Null distribution **\n");
  }
  /** Estimate the Null distributions **/
  MATNullDistribution(position, nProbes, dMax, MATScores, &sigma0, &mu0, seqNum);

  if(*verbose)
  {
    printf("** Calculate P-values **\n");  
  }
  /** Compute the associated pValues **/
  MATpValue(*nProbes, MATScores, sigma0, mu0, pValues);
  
  /** Compute the FDR cutoff **/
  if(*method==1) /** Based on MAT scores **/
  {
    numRegions[0] = mergeMATScores(position, *nProbes, *dMerge, MATScores, mu0, *threshold, 1, regions, seqNum);
  }
  else if(*method==2)
  {
    numRegions[0] = mergeMATScores(position, *nProbes, *dMerge, pValues, 0, *threshold, -1, regions, seqNum);
  }
  else if(*method==3)
  {
    cutoff=MATcutoffFDR(position, *nProbes, *dMerge, MATScores, mu0, *threshold, regions, seqNum);
    // printf("** The FDR cutoff is %lf **\n",cutoff);

      /** Form the regions using the threshold found **/
    numRegions[0] = mergeMATScores(position, *nProbes, *dMerge, MATScores, mu0, cutoff, 1, regions, seqNum);
  }
/*printf("num regions is %d", numRegions[0]);*/
}

void MATScore(double *C, double *I, int *nProbes, int *nArraysC, int *nArraysI, int *position, double *dMax, double *MATScores, int *seqNum)
{
  int pMin=0,pMax=0;
  double cMin=0,cMax=0;
  double MC=0,MI=0;
  double *dataInRegion;
  gsl_vector *vectorInRegion;
  gsl_vector_view vectorTmp;
  int nProbesRegion=0, nProbesNotTrimmed=0;
  int i,j,k;

  
  pMin=0;
  pMax=0;

  for(i=0;i<*nProbes;i++)
  {
    /* Current implementation is such that it is within 600 bp from pMin to i and another 600bp from i to pMax*/
    if((seqNum[pMin]!=seqNum[i]) && (seqNum[pMax]!=seqNum[i]))
    {
      pMin = i;
      pMax = i;
    }
    /* with this implementation, pMax and pMin must stay within the same sequence as i. The only time pMax and pMin are not the same as i is when moves out to the next sequence, and by then pMax and pMin should both be set as i to initialize the beginning of the new sequence. The else if clause checks the improbable state when either one of seqNum[pMin] or seqNum[pMax] is not the same as i*/
    else if((seqNum[pMin]!=seqNum[i])||(seqNum[pMax]!=seqNum[i]))
    {
      error("Check that your intensities are ordered by chromosome then by position \n");
    }
    /* It was sorted.. so no need abs(); pMin keeps moving WRT to i */
    while((pMin < i) &&  ((position[i] - position[pMin])>*dMax/2.) && (seqNum[pMin]==seqNum[i]))
    {
      pMin++;
    }
    /* should be position[pMin+1] so it is always within */
    while((pMax<*nProbes) && ((position[pMax+1]-position[i])<=*dMax/2.) && (seqNum[pMax+1]==seqNum[i]) && (pMax+1 < *nProbes))
    {
      pMax++;
    }
    /* If there is only one probe we skip it */
    if(pMax-pMin>0)
    {
      MC=0,MI=0;
      nProbesNotTrimmed=0;

      if(*nArraysC>0)
      {
        /*If information about control is available, we should integrate it into our study */
        nProbesRegion=(pMax-pMin)**nArraysC;
        dataInRegion=C+pMin**nArraysC;
        /** Vector view of the window data **/
        vectorTmp=gsl_vector_view_array(dataInRegion, nProbesRegion);
        /** Allocate a vector to store the probe measurements within the window **/
        vectorInRegion=gsl_vector_alloc(nProbesRegion);
        /** Copy the probe measurements within the window **/
        gsl_vector_memcpy(vectorInRegion,&vectorTmp.vector);
        /** Sort the probe measurements within the window **/
        gsl_sort_vector(vectorInRegion);
        /** Compute the .1 and .9 quantiles **/
        cMin=gsl_stats_quantile_from_sorted_data(vectorInRegion->data,1, nProbesRegion, .1);
        cMax=gsl_stats_quantile_from_sorted_data(vectorInRegion->data,1, nProbesRegion, .9);
        /** Desallocate the memory **/
        gsl_vector_free(vectorInRegion);
        /** Computing the Control's Trimmed Mean **/
        for(k=pMin;k<pMax;k++)
        {
          for(j=0;j<*nArraysC;j++)
          {
            if((C[k**nArraysC+j]>=cMin) & (C[k**nArraysC+j]<=cMax))
            {
              MC+=C[k**nArraysC+j];
              nProbesNotTrimmed++;
            }
          }
        }
        if(nProbesNotTrimmed>0)
          MC=MC/nProbesNotTrimmed;
      }

      /* For the Immunoprecipitated Array*/
      nProbesNotTrimmed=0;
      nProbesRegion=(pMax-pMin)**nArraysI;
      dataInRegion=I+pMin**nArraysI;

      /** Vector view of the window data **/
      vectorTmp=gsl_vector_view_array(dataInRegion, nProbesRegion);
      /** Allocate a vector to store the probe measurements within the window **/
      vectorInRegion=gsl_vector_alloc(nProbesRegion);
      /** Copy the probe measurements within the window **/
      gsl_vector_memcpy(vectorInRegion,&vectorTmp.vector);
      /** Sort the probe measurements within the window **/
      gsl_sort_vector(vectorInRegion);
      /** Compute the .1 and .9 quantiles **/
      cMin=gsl_stats_quantile_from_sorted_data(vectorInRegion->data,1, nProbesRegion, .1);
      cMax=gsl_stats_quantile_from_sorted_data(vectorInRegion->data,1, nProbesRegion, .9);

      /** Desallocate the memory **/
      gsl_vector_free(vectorInRegion);
      /** Computing the IP treated sample's trimmed mean **/
      for(k=pMin;k<pMax;k++)
      {
        for(j=0;j<*nArraysI;j++)
        {
          if((I[k**nArraysI+j]>=cMin) & (I[k**nArraysI+j]<=cMax))
          {
            MI+=I[k**nArraysI+j];
            nProbesNotTrimmed++;
          }
        }
      }

      /** Computing the MAT Score **/
      if(nProbesNotTrimmed>0)
      {
        MI=MI/nProbesNotTrimmed;
        MATScores[i]=sqrt((pMax-pMin)*.8)*(sqrt(*nArraysI)*MI-sqrt(*nArraysC)*MC);
      }
      else
      {
        MATScores[i]=0.0;
      }
    }
    else
    {
      MATScores[i]=0.0;
    }
  }
}

void MATNullDistribution(int *position, int *nProbes, double *dMax, double *MATScores, double *sigma0, double *mu0, int *seqNum)
{
  /** Number of regions used to estimate the null distribution **/
  /**Assume position vector sorted **/
  /**  int nRegionsMax=(int) (((position[nProbes-1]-position[0])/dMax)+5),nRegions=0; **/
  int nRegionsMax = 0, nRegions=0;
  double totalSeqsPosition=0;
  int i,p0,p1;
  int curSeqStartPos=-1;
  int curSeqNum=-1;
  
  for(i=0;i<*nProbes;i++)
  {
    if(curSeqNum!=seqNum[i])
    {
      curSeqNum=seqNum[i];
      curSeqStartPos=position[i];
    }
    if((i+1==*nProbes) | (curSeqNum!=seqNum[i+1]))
    {
      totalSeqsPosition+=(double)(position[i]-curSeqStartPos);
    }
    nRegionsMax = (int)(totalSeqsPosition/(*dMax))+5;
  }

  gsl_vector *NoOverLapMATScores=gsl_vector_calloc(nRegionsMax);
  gsl_vector *NoOverLapMATScoresMinus;

  p0=0;p1=0;
  while(p1<*nProbes)
  {
      /** If a score is zero, keep going **/
    while((p1<*nProbes) && (MATScores[p1]==0))
    {
      p1++;
    }
    gsl_vector_set(NoOverLapMATScores,nRegions,MATScores[p1]);
    nRegions++;
    
    while((p1<*nProbes) && (position[p1]-position[p0])<=(*dMax) && (seqNum[p1]==seqNum[p0]))
    {
      p1++;
    }
    p0=p1;	    
  }

  gsl_sort(NoOverLapMATScores->data, 1, nRegions);
  *mu0=gsl_stats_median_from_sorted_data(NoOverLapMATScores->data, 1, nRegions);
  /*\ Half the Data? */
  /** Multiply by sqrt(2) since we estimated the variance from half of the data **/
  
  NoOverLapMATScoresMinus=gsl_vector_calloc(nRegions);
  
  for(i=0;i<(int)nRegions/2;i++)
  {
    gsl_vector_set(NoOverLapMATScoresMinus,i,gsl_vector_get(NoOverLapMATScores,i));
  }
  for(i=(int)nRegions/2;i<nRegions;i++)
  {
    gsl_vector_set(NoOverLapMATScoresMinus,i,2**mu0-gsl_vector_get(NoOverLapMATScores,i-(int)nRegions/2));
  }
  // gsl_sort(NoOverLapMATScoresMinus->data, 1, nRegions);
  *sigma0=gsl_stats_sd(NoOverLapMATScoresMinus->data, 1, nRegions);
  // *sigma0=sqrt(2)*gsl_stats_sd(NoOverLapMATScores->data, 1, nRegions/2);
  
  gsl_vector_free(NoOverLapMATScores);
  gsl_vector_free(NoOverLapMATScoresMinus);
}

void MATpValue(int nProbes, double *MATScores, double sigma0, double mu0, double *pValues)
{
  /** Number of regions used to estimate the null distribution **/
  int i;
  for(i=0;i<nProbes;i++)
  {
    pValues[i]=1-gsl_cdf_gaussian_P(MATScores[i]-mu0, sigma0);
  }
}

double MATcutoffFDR(int *position, int nProbes, double dMerge, double *MATScores, double mu0, double FDR, int *regions, int *seqNum)
{

  /** Number of regions used to estimate the null distribution **/
  double step=0.05, estimatedFDR=1;
  int nPositive=0,nNegative=0;

  // int INITIALSTEP=1;
  // int notInitial=0;
  double proposeCutoff=0.1;

  while((proposeCutoff<50) && (estimatedFDR > FDR))
  {
    nPositive=mergeMATScores(position, nProbes, dMerge, MATScores, mu0, proposeCutoff, 1, regions, seqNum);
    nNegative=mergeMATScores(position, nProbes, dMerge, MATScores, mu0, -proposeCutoff, -1, regions, seqNum);
    if(nPositive>0)
    {
      estimatedFDR=fmin2(nNegative/((double)nPositive),1.0);
    }
    else
    {
      estimatedFDR=0;
    }
    /* Increase the cutoff */
    proposeCutoff=proposeCutoff+step;
 
  }
  return(proposeCutoff);  
}



int mergeMATScores(int *position, int nProbes, double dMerge, double *MATScores, double m0, double cutoff, int sign, int *regions, int *seqNum)
{
  int i=0,p=0;
  int w=0,ww=0,max=0;
  int nRegions=0;

  while(p<nProbes)
  {

    /*First detection of a MAT Region */
    if(((MATScores[p]-m0>cutoff) && (sign==1)) | ((MATScores[p]-m0<cutoff) && (sign==-1)))
    {
      nRegions++;
      regions[p]=nRegions;
      w=p+1;
      max=p;
      ww=p;
      /*Scan through the whole MAT Region */
      /** Here we make sure we are on the same chromosome **/
      while((w<nProbes) && ((position[w]-position[ww])<=dMerge) && (seqNum[w]==seqNum[ww]))
      {
        if(((MATScores[w]-m0>cutoff) && (sign==1)) || ((MATScores[w]-m0<cutoff) && (sign==-1)))
        {
          max=w;
          ww=max;
        }
        w++;
      }
      for(i=p;i<=max;i++)
      {
        regions[i]=nRegions;
      }
      p=w;
    }
    else
    {
      regions[p]=0;
      p++;
    }
  }
  return(nRegions);
}

void normArray(char **seq, double *y, double *yNormalized, int *nProbes, int *nArrays, double *copyNumber, int *method, 
               int *robust, double *adjRSquare, double *RSquare, double *BIC, double *outBeta, int *betaLength, int *all,  int *MATScaling,int *isVerbose, 
               gsl_matrix *pairNumCount1, gsl_matrix *pairNumCount2, gsl_matrix *pairNumCount3, gsl_matrix *pairNumCount4,
               gsl_matrix *seqNumCount,
               int nProbesTotal,
               int nVariables,
               int nBins,
               int nProbesPerBin,
               gsl_vector_view yVector, gsl_vector_view copyNumberVector,
               int j)
{
  double RSS=0,TSS=0,meanY=0;
  int i=0,k=0,iterMax=10;
  gsl_vector_view xRow;
  gsl_vector *beta, *betaTmp, *weight, *fittedSorted;
  gsl_vector *sdBins=NULL, *yBin=NULL, *XR;
  gsl_permutation *indexSorted;
  gsl_matrix *X, *H;

  weight=gsl_vector_calloc(*nProbes);
  fittedSorted=gsl_vector_calloc(nProbesTotal);

  if((*MATScaling)==1)
  {
    sdBins=gsl_vector_calloc(nBins);
    yBin=gsl_vector_calloc(nProbesPerBin+nProbesTotal-nProbesPerBin*(nBins-1)+1);
  }

  copyNumberVector=gsl_vector_view_array(copyNumber,nProbesTotal);
  XR=gsl_vector_calloc(nVariables);

  if(*method==1) /** This is the MAT model **/
  {
    nVariables=1+3*25+4+1;
    X=gsl_matrix_calloc(*nProbes,nVariables);
    createDesignMatrixMAT(seqNumCount, &copyNumberVector.vector, X, seq);
  }
  else  /** This is the pair binned model --into 4 bins**/
  {
    nVariables=73;
    X=gsl_matrix_calloc(*nProbes, nVariables);
    createDesignMatrixPairBinned(seqNumCount, 
      pairNumCount1,
      pairNumCount2,
      pairNumCount3,
      pairNumCount4,
      &copyNumberVector.vector, X);
  }

  /*****************************************************************************************************/

  beta=gsl_vector_alloc(nVariables);
  betaTmp=gsl_vector_alloc(nVariables);
  H=gsl_matrix_calloc(nVariables,nVariables);
  indexSorted=gsl_permutation_calloc(nProbesTotal);

  gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X, 0, H);
  /** Cholesky decomposition of H **/
  gsl_linalg_cholesky_decomp(H);

  if(*isVerbose)
  {
    printf("** Start normalization of array %d **\n",j);
  }

  /** Select the array **/
  yVector=gsl_vector_view_array(y+j*nProbesTotal, *nProbes);
  /** Compute X'y **/
  gsl_blas_dgemv(CblasTrans, 1.0, X, &yVector.vector, 0.0, betaTmp);

  /** Compute (X'X)^{-1}X'y **/
  gsl_linalg_cholesky_solve(H, betaTmp, beta);

  meanY=gsl_stats_mean(y+j*nProbesTotal,1,*nProbes);

  for(i=0;i<nProbesTotal;i++)
  {
    if(i<(*nProbes))
    {
      xRow=gsl_matrix_row(X,i);
      /** Compute the fitted data X(X'X)^{-1}X'y **/
      gsl_blas_ddot(&xRow.vector, beta, yNormalized+j*nProbesTotal+i);
      RSS+=gsl_pow_2(y[j*nProbesTotal+i]-yNormalized[j*nProbesTotal+i]);
      TSS+=gsl_pow_2(y[j*nProbesTotal+i]-meanY);
    }
    else
    {
      /** Need to get the corresponding row of the design matrix **/
      if(*method==1)
      {
        createDesignMatrixMATRow(seqNumCount, &copyNumberVector.vector, XR, i, seq);
      }
      else
      {
        createDesignMatrixPairBinnedRow(seqNumCount, 
          pairNumCount1,
          pairNumCount2,
          pairNumCount3,
          pairNumCount4,
          &copyNumberVector.vector, XR,i);
      }
      gsl_blas_ddot(XR, beta, yNormalized+j*nProbesTotal+i); 
    }
  }

  
  RSquare[j]=1.0-RSS/TSS;
  adjRSquare[j]=1.0-(nProbesTotal-1.0)/(nProbesTotal-nVariables-1.0)*RSS/TSS;
  BIC[j]=-*nProbes*log(RSS/(*nProbes))-nVariables*log(*nProbes);

  // gsl_vector_fprintf(stdout, beta, "%.5g");

  /* If we specify to use the robust estimate */
  if(*robust==1)
  {
    if(*isVerbose){
      printf("** Start robust normalization of array %d **\n",j);      
    }
    
    for(k=0;k<=iterMax;k++)
    {
      for(i=0;i<*nProbes;i++)
      {
        /** Compute the weights (with 4 degrees of freedom) **/
        gsl_vector_set(weight,i,(4.0+1.0)/(4.0+*nProbes/RSS*gsl_pow_2(y[j*nProbesTotal+i]-yNormalized[j*nProbesTotal+i])));
        /** Multiply each row of X by the weight **/
        xRow=gsl_matrix_row(X,i);
        gsl_vector_scale(&xRow.vector, sqrt(gsl_vector_get(weight,i)));
        /** Multiply each observation by its weight **/
        y[j*nProbesTotal+i]=sqrt(gsl_vector_get(weight,i))*y[j*nProbesTotal+i];
      }

      gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, X, 0, H);
      /** Compute X'y **/
      gsl_blas_dgemv(CblasTrans, 1.0, X, &yVector.vector, 0.0, betaTmp);
      /** Cholesky decomposition of H **/
      gsl_linalg_cholesky_decomp(H);
      /** Compute (X'X)^{-1}X'y **/
      gsl_linalg_cholesky_solve(H, betaTmp, beta);

      RSS=0;TSS=0;
      for(i=0;i<nProbesTotal;i++)
      {
        if(i<(*nProbes))
        {
          xRow=gsl_matrix_row(X,i);
          /** Compute the fitted data X(X'X)^{-1}X'y **/
          gsl_blas_ddot(&xRow.vector, beta, yNormalized+j*nProbesTotal+i);
          /** reweight the data **/
          y[j*nProbesTotal+i]=y[j*nProbesTotal+i]/sqrt(gsl_vector_get(weight,i));
          yNormalized[j*nProbesTotal+i]=yNormalized[j*nProbesTotal+i]/sqrt(gsl_vector_get(weight,i));
          RSS+=gsl_pow_2(y[j*nProbesTotal+i]-yNormalized[j*nProbesTotal+i]);
          TSS+=gsl_pow_2(y[j*nProbesTotal+i]-meanY);          
          /** reweight X **/
          gsl_vector_scale(&xRow.vector, 1./sqrt(gsl_vector_get(weight,i)));
        }
        else if(k==iterMax & i>=(*nProbes))
        {
          /** Need to get the corresponding row of the design matrix **/
          if(*method==1)
          {
            createDesignMatrixMATRow(seqNumCount, &copyNumberVector.vector, XR, i, seq);
          }
          else
          {
            createDesignMatrixPairBinnedRow(seqNumCount,
              pairNumCount1,
              pairNumCount2,
              pairNumCount3,
              pairNumCount4,
              &copyNumberVector.vector, XR,i);
          }
          gsl_blas_ddot(XR, beta, yNormalized+j*nProbesTotal+i);
        }
      }
    }
    /** Compute the residual and total SS **/
    RSquare[j]=1.0-RSS/TSS;
    adjRSquare[j]=1.0-(nProbesTotal-1.0)/(nProbesTotal-nVariables-1.0)*RSS/TSS;
    BIC[j]=-*nProbes*log(RSS/(*nProbes))-nVariables*log(*nProbes);
  }

    // gsl_vector_fprintf(stdout, beta, "%.5g");

  /** Need to scale the values for MAT **/
  if((*MATScaling)==1)
  {
    if(*isVerbose)
    {
      printf("** Start scaling of array %d **\n",j);
    }
    gsl_permutation_init(indexSorted);
    for(i=0;i<nProbesTotal;i++)
      gsl_vector_set(fittedSorted, i, yNormalized[j*nProbesTotal+i]);
        /** Sort the fitted values **/
    gsl_sort_vector_index(indexSorted, fittedSorted);
    /** Compute the standard deviation per bin and scale the values **/
    for(i=0;i<(nBins-1);i++)
    {
      for(k=0;k<nProbesPerBin;k++)
        gsl_vector_set(yBin,k,y[j*nProbesTotal+indexSorted->data[nProbesPerBin*i+k]]);
      gsl_vector_set(sdBins,i,gsl_stats_sd(yBin->data,1,nProbesPerBin));

      for(k=0;k<nProbesPerBin;k++)
        yNormalized[j*nProbesTotal+indexSorted->data[nProbesPerBin*i+k]]=(y[j*nProbesTotal+indexSorted->data[nProbesPerBin*i+k]]-yNormalized[j*nProbesTotal+indexSorted->data[nProbesPerBin*i+k]])/gsl_vector_get(sdBins,i);
    }
        /** For the last bin, we actually used a bit more probes **/
    for(k=0;k<(nProbesTotal-nProbesPerBin*(nBins-1));k++)
      gsl_vector_set(yBin,k,y[j*nProbesTotal+indexSorted->data[nProbesPerBin*(nBins-1)+k]]);
    gsl_vector_set(sdBins,nBins-1,gsl_stats_sd(yBin->data,1,nProbesTotal-nProbesPerBin*(nBins-1)));
    for(k=0;k<(nProbesTotal-nProbesPerBin*(nBins-1));k++)
      yNormalized[j*nProbesTotal+indexSorted->data[nProbesPerBin*(nBins-1)+k]]=(y[j*nProbesTotal+indexSorted->data[nProbesPerBin*(nBins-1)+k]]-yNormalized[j*nProbesTotal+indexSorted->data[nProbesPerBin*(nBins-1)+k]])/gsl_vector_get(sdBins,nBins-1);
  } 
  else /** if no MAT Scaling for MAT or not using MAT **/
  {
    for(i=0;i<nProbesTotal;i++)
    {
      yNormalized[j*nProbesTotal+i]=y[j*nProbesTotal+i]-yNormalized[j*nProbesTotal+i];
    }
  }
  /*****************************************************************************************************/ 

  *betaLength = nVariables;
  for(k=0; k<nVariables; k++)
  {
    outBeta[k+j*nVariables] = beta->data[k];
  }
  /*****************************************************************************************************/
  
  gsl_vector_free(beta);
  gsl_vector_free(betaTmp);
  gsl_matrix_free(X);
  gsl_matrix_free(H);
  gsl_vector_free(weight);
  gsl_vector_free(fittedSorted);
  if((*MATScaling==1))
  {
    gsl_vector_free(sdBins);
    gsl_vector_free(yBin);
  }
  gsl_permutation_free(indexSorted);
  gsl_vector_free(XR);
}

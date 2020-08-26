#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "wave.h"
#include "LBG.h"
#include "PCA.h"

void main(void)
{
	double **fv, **cb, **m;
	int i,j,k,n=0;

	
	/*** read wav file and cut the speech signal into several segments with a certain length ***/
	readWavFile("wav/scientific america.wav", 20);//the first parameter is the address of input wav file and the second parameter is the length of each segment 

	
	/*** feature extraction ***/
	fv=getFeatureVectors(100,0.9);//the first parameter is the threshold to remove useless components and the second parameter is correlation threshold
	
	
	/*** optimized LBG training & codebook adjustment ***/
	cb=training(fv,0.4,0.5,10);//the second parameter is variance threshold, the third is distortion threshold and the last is cluster size threshold

	
	/*** combine fv and cb to form a big matrix ***/
	m=(double**)malloc((no_segments+N)*sizeof(double*));
	for(i=0;i<(no_segments+N);i++)
	{
		m[i]=(double*)malloc(K*sizeof(double));
	}	
	
	for(i=0;i<no_segments;i++)
	{
		for(k=0;k<K;k++)
		{
			m[i][k]=fv[i][k];
		}
	}
	
	for(i=no_segments;i<(no_segments+N);i++)
	{
		for(j=0;j<K;j++)
		{
			m[i][j]=cb[i-no_segments][j];
		}
	}
	
	
	/*** PCA 2d plot ***/
	printf("start PCA 2d\n");
	getMatrix(m, no_segments+N, K);
	zeroMean();
	covMatrix();
	findEigValue();
	get2dCoordinate();
	release();
	printf("2d coordinates generated\n");

		
	/*** Re-arranged test wav file***/
	writeTestWavFile();
	
	
	releaseMemory();
}
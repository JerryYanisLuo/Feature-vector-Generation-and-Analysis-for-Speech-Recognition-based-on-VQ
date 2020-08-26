#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "wave.h"
#include "LBG.h"

Wav				infile;				//the header of the input wav file
FILE			*wav;				//input wav file
short 			**wav_data;			//wav sample data
int				no_sample_ps;		//number of samples per segment
int				no_segments;		//number of segments
unsigned short	size_sample;		//in bytes
double 			**autocorrelation;	//the result of acf
double			**fv;				//featrure vector
int				K;					//dimension of feature vector

/*
	Function description: 
	read wav file, print the info of wav header and record all sample data in **wav_data
	
	Parameters: 
	@char wavfile[] : the name of wav file
	@int delta_t 	: the period of each segment, usually between 10 - 20 ms
	
	return value	: none
*/
void readWavFile(char wavfile[], int delta_t)
{
	double			duration;		//total duration of wav file in seconds
	unsigned char	*riff;
	unsigned long	riffsize;		//size of overall file in bytes
	unsigned char 	*wave;			//"wavefmt"
	unsigned long	hex;			//length of format data listed above
	unsigned short	pcm;			//Type of format (1 is PCM)
	unsigned short	nochannel;		//Number of Channels
	unsigned long	fs;				//Sample Rate
	unsigned long 	fsch;			//(Sample Rate * BitsPerSample * Channels) / 8
    unsigned short 	ch;				//(BitsPerSample * Channels) / 8
    unsigned short 	bps;			//BitsPerSample
    unsigned char 	*data;			//"data"	
	unsigned long	no_sample;		//size of data chunk
	int				i,j;			//used in for loop
	
	//get wav file
    wav=fopen(wavfile,"rb");
	
	//read wav head
    fread(&infile, sizeof(infile), 1, wav);
	
	//prepare for some property data
	riff = infile.riff;
	riffsize=infile.riffsize;
	wave=infile.wave;
	hex=infile.hex;
	pcm=infile.pcm;
	nochannel=infile.nochannel;
    fs=infile.fs;
	fsch=infile.fsch;
	ch=infile.ch;
	bps=infile.bps;
	data=infile.data;
	no_sample=infile.Subchunk2Size*8/bps;
    size_sample=infile.bps/8;
 
	//print wav head info
	printf("riff     : %c%c%c%c\n",  	riff[0],riff[1],riff[2],riff[3]);
    printf("riffsize : %d\n",   		riffsize);
    printf("wave     : %c%c%c%c\n",   	wave[0],wave[1],wave[2],wave[3]);
	printf("fmt      : %c%c%c%c\n",   	wave[4],wave[5],wave[6],wave[7]);
    printf("hex      : %d\n",   		hex);
    printf("pcm      : %d\n",  			pcm);
    printf("nochannel: %d\n",   		nochannel);
    printf("fs       : %d\n",   		fs);
    printf("fsch     : %d\n",   		fsch);
    printf("ch       : %d\n",   		ch);
    printf("bps      : %d\n",   		bps);
    printf("data     : %c%c%c%c%\n",   	data[0],data[1],data[2],data[3]);
    printf("nosample : %d\n", 			no_sample);
	printf("--------------------------------------------------\n");

	//calculate the total duration of the wav file
	duration=1.0*riffsize/fsch;

    //calculate number of samples in t ms
	no_sample_ps=fs*nochannel*delta_t/1000;
	
	//adjust number of samples per segment
	if(nochannel==2)
	{
		if(no_sample_ps%4!=0)
		{
			no_sample_ps=no_sample_ps+no_sample_ps%4;
		}
	}

	//calculate number of segments
	no_segments=no_sample/no_sample_ps;

	printf("duration                     : %f seconds\n",  duration);
	printf("segment period               : %d ms\n",  delta_t);
	printf("number of samples per segment: %d\n",  no_sample_ps);
	printf("number of segments           : %d\n",  no_segments);
	printf("--------------------------------------------------\n");
	
	//define data array
	wav_data=(short**)malloc(no_segments*sizeof(short*));
	for(i=0;i<no_segments;i++)
	{
		wav_data[i]=(short*)malloc(no_sample_ps*sizeof(short));
	}
	
	//read wav data
	printf("Start reading wav data\n");	
	for(i=0;i<no_segments;i++)
	{
		fread(wav_data[i], size_sample, no_sample_ps, wav);
	}
	printf("Finish reading wav data\n");
	
	FILE *wavData;
	wavData = fopen("wav data/wav_data.txt","w");
	
	//check data
	for(i=0;i<no_segments;i++)
	{
		for(j=0;j<no_sample_ps;j++)
		{
			fprintf(wavData,"%d ",wav_data[i][j]);
			if(wav_data[i][j]>32767)
			{
				printf("Error data %d\n",wav_data[i][j]);
				wav_data[i][j]=32767;
				printf("Data corrected\n");
			}
			else if(wav_data[i][j]<-32768)
			{
				printf("Error data %d\n",wav_data[i][j]);
				wav_data[i][j]=-32768;
				printf("Data corrected\n");
			}
		}
	}
	fclose(wavData);
}



/*
	Function description: 
	calculate the variance of the given dimension of the input signal
	
	Parameters: 
	@double *signal	: the input signal
	@int size		: the number of signals
	
	return value:
	@double v		: the calculated variance
*/
double getVariance(double *signal, int size)
{
	double temp=0, mean, v = 0;
	int i;
	
	//calculate mean value
	for(i=0;i<size;i++)
	{
		temp+=signal[i];
	}
	mean=temp/size;

	//calculate variance
	for(i=0;i<size;i++)
	{
		v += (signal[i]-mean)*(signal[i]-mean);
	}
	v=v/size;
	
	return v;
}



/*
	Function description: 
	calculate the normalized correlation of signel x and y. x and y have the same length
	
	Parameters: 
	@double *x	: signel x
	@double *y	: signel y
	@int size	: the size of each signal
	
	return value:
	@double corr: the result correlation value
*/
double normCorr(double *x, double *y, int size)
{
	double numerator;
	double denominator;
	double corr;
	double x2=0, y2=0, xy=0;
	int i,j,k;
	
	for(i=0;i<size;i++)
	{
		xy+=x[i]*y[i];
		x2+=x[i]*x[i];
		y2+=y[i]*y[i];
	}
	
	numerator = xy;
	denominator = sqrt(x2*y2);
	corr = numerator/denominator;
	
	return corr;
}



/*
	Function description: 
	remove one given vector from **autocorrelation
	
	Parameters: 
	@int r		: the index of the component needed to be removed
	
	return value: none
*/
void removeComponent(double **a, int r)
{
	int i,j;
	for(i=0;i<no_segments;i++)
	{
		for(j=r;j<K-1;j++)
		{
			a[i][j]=a[i][j+1];
		}
	}
	K--;
}



/*
	Function description: 
	signal normalization with rout mean squre 
	
	Parameters: 
	@double *signal	: the input signal
	@int k			: the length of the signal
	
	return value: none
*/
void normalizationSignal(double *signal,int k)
{
	int i;
	double rms=0;
	for(i=0;i<k;i++)
	{
		rms+=signal[i]*signal[i];
	}
	
	rms=sqrt(rms/k);

	for(i=0;i<k;i++)
	{
		signal[i]=signal[i]/rms;
	}
}


/*
	Function description: 
	transform the input signal to power spectrum based on discret fourier transform
	
	Parameters: 
	@double *signal	: the input signal
	@int num		: the length of the signal
	
	return value: none
*/
void dft(double *signal, int num)
{
	int i,j;
	double f[num];
	double sin_x=0,cos_x=0;
	for(i=0;i<num;i++)
	{
		sin_x=0;
		cos_x=0;
		for(j=0;j<num;j++)
		{
			sin_x+=signal[j]*cos(2*3.14/num*j*i);
			cos_x+=signal[j]*sin(2*3.14/num*j*i);
		}
		f[i]=(sin_x*sin_x+cos_x*cos_x);
	}
	
	for(i=0;i<num;i++)
	{
		signal[i]=f[i];
	}
}



/*
	Function description: 
	calculate normalized autocorrelation function for each segment of sample data, and record the result in **autocorrelation
	
	Parameters:
	@double *data	: the input sample data
	@int num_seg	: the total number of segments
	@dimension		: the number of samples within each segment
	
	return value	: none
*/
void acf(double **data, int num_seg, int dimension)
{
	int i,n; //counter in loop
	int t; //time shift
	double max;

	autocorrelation=(double**)malloc(num_seg*sizeof(double*));
	for(i=0;i<num_seg;i++)
	{
		autocorrelation[i]=(double*)malloc(dimension*sizeof(double));
	}
	
	for(i=0;i<num_seg;i++)
	{
		for(t=0;t<dimension;t++)
		{
			autocorrelation[i][t]=0;
			for(n=0;n<dimension;n++)
			{
				if(n>=t)
				{
					autocorrelation[i][t]+=1.0*data[i][n]*data[i][n-t];
				}
			}
			autocorrelation[i][t] = autocorrelation[i][t]/(1.0*dimension);
		}

		max = autocorrelation[i][0];
		for(t=0;t<dimension;t++)
		{
			if(max!=0)
			{
				autocorrelation[i][t]=autocorrelation[i][t]/max;
			}
		}
	}
	printf("acf finished\n");
}



/*
	Function description: 
	add Hamming window function to the input signal
	
	Parameters:
	@double *data	: the input signal
	@int n			: the length of the signal
	
	return value	: none
*/
void hammingWin(double *data, int n)
{
	int i;
	for(i=0;i<n;i++)
	{
		data[i]=(0.54-0.46*cos(2*3.1415926*(i-(n-1)/2)/(n-1)))*data[i];
	}
}



/*
	Function description: 
	add Hanning window function to the input signal, and record the result in *w
	
	Parameters:
	@short *data	: the input signal
	@int n			: the length of the signal
	@double *w		: the output signal
	
	return value	: none
*/
void hanningWin(short *data, int n,double *w)
{
	int i;
	for(i=0;i<n;i++)
	{
		w[i]=(0.5-0.5*cos(2*3.1415926*i/(n-1)))*data[i];
	}
}


/*
	Function description: 
	add Blackman window function to the input signal, and record the result in *w
	
	Parameters:
	@short *data	: the input signal
	@int n			: the length of the signal
	@double *w		: the output signal
	
	return value	: none
*/
void blackmanWin(short *data, int n, double *w)
{
	int i;
	for(i=0;i<n;i++)
	{
		w[i]=(0.42+0.5*cos(2*3.1415926*(i-(n-1)/2)/(n-1))+0.08*cos(4*3.1415926*(i-(n-1)/2)/(n-1)))*data[i];
	}
}



/*
	Function description: 
	add Rectangular window function to the input signal, and record the result in *w
	
	Parameters:
	@short *data	: the input signal
	@int n			: the length of the signal
	@double *w		: the output signal
	
	return value	: none
*/
void rectWin(short *data, int n,double *w)
{
	int i;
	for(i=0;i<n;i++)
	{
		w[i]=data[i];
	}
}



/*
	Function description: 
	derive feature vectors from the result of autocorrelation, 
	decide the dimension and components of fv according to the 
	variance and correlation between fv
	
	Parameters:
	@int n_threshold				: keep the first n components after removing useless components
	@double correlation_threshold	: correlation threshold for removing redundant components
	
	return value: 
	@double **fv					: feature vectors
*/
double** getFeatureVectors(int n_threshold, double correlation_threshold)
{
	int i,j,k,n,m;
	double **w;
	w=(double**)malloc(no_segments*sizeof(double*));
	for(i=0;i<no_segments;i++)
	{
		w[i]=(double*)malloc(no_sample_ps*sizeof(double));
	}

	

	/************************** windowing **************************/
	for(i=0;i<no_segments;i++)
	{
		rectWin(wav_data[i],no_sample_ps,w[i]);
		//hanningWin(wav_data[i],no_sample_ps,w[i]);
		//blackmanWin(wav_data[i],no_sample_ps,w[i]);
	}
	printf("adding window\n");
	

	
	/***************************** acf *****************************/
	acf(w, no_segments, no_sample_ps);

	
	
	/********************* calculate variance **********************/
	K=no_sample_ps;
	double variance[K];
	double signal[no_segments];
	double temp;
	int r[K];
	for(j=0;j<K;j++)
	{	
		for(n=0;n<no_segments;n++)
		{
			signal[n]=autocorrelation[n][j];
		}
		variance[j]=getVariance(signal,no_segments);
	}

	
	
	/*********** descend the result of acf via variance ************/
	for(j=0;j<K;j++)		
	{
		for(i=j+1;i<K;i++)
		{
			if(variance[j]<variance[i])
			{
				temp=variance[j];
				variance[j]=variance[i];
				variance[i]=temp;
				
				for(n=0;n<no_segments;n++)
				{
					temp=autocorrelation[n][j];
					autocorrelation[n][j]=autocorrelation[n][i];
					autocorrelation[n][i]=temp;
				}
			}
		}
	}
	printf("components ordering\n");
	
	
	/*********** keep the first n components ************/
	K=n_threshold;
	
	
	/**** calculate normalized correlation of each **autocorrelation to adjust the dimension K of fv *****/
	/*** in this case, the generated feature vectors will differ from each other to the largest extent ***/
	double signalx[no_segments], signaly[no_segments];
	double corr;
	for(i=0;i<K;i++)
	{
		r[i]=0;
	}
	
	for(i=0;i<K;i++)
	{
		for(j=i+1;j<K;j++)
		{
			for(n=0;n<no_segments;n++)
			{
				signalx[n]=autocorrelation[n][i];
				signaly[n]=autocorrelation[n][j];
			}
			corr=normCorr(signalx, signaly, no_segments);
			
			if(corr>correlation_threshold)
			{
				r[j]=1;
			}
		}

		for(m=K-1;m>=0;m--)
		{
			if(r[m]==1)
			{
				removeComponent(autocorrelation,m);
			}
			r[m]=0;
		}	
	}
    printf("K:%d\n",K);

	
	
	/************* get feature vectors and record all feature vectors in txt file *************/
	FILE *file_fv;
	file_fv=fopen("K-fv/Kfv.txt","w");
	fv=(double**)malloc(no_segments*sizeof(double*));
	for(i=0;i<no_segments;i++)
	{
		fv[i]=(double*)malloc(K*sizeof(double));
	}
	
	for(i=0;i<no_segments;i++)
	{
		for(j=0;j<K;j++)
		{
			fv[i][j]=autocorrelation[i][j];
			fprintf(file_fv,"%f ",fv[i][j]);
		}
		if(i!=(no_segments-1))
		{
			fprintf(file_fv,"\n");
		}
	}
	fclose(file_fv);
	printf("Feature vectors generated\n");
	return fv;
}



/*
	Function description: 
	analyze similarity of each cluster by measuring the segment variance
	
	Parameters:
	@double *noSet: the set of the id of each feature vector
	@double *fvSet: mirror the id of feature vector and the id of cluster, 
					with the index representing the id of feature vector, 
					and the value representing the id of cluster
	
	return value: none
*/
void segmentVariance(int *noSet, int *fvSet)
{
	/****************************** compute segment variance  ***********************************/
	/*** re-arrange sample date based on clusters ***/
	int i,j,k;
	double ***clu_data;
	clu_data=(double***)malloc(50*sizeof(double**));
	for(i=0;i<50;i++)
	{
		clu_data[i]=(double**)malloc(no_segments*sizeof(double*));
		for(j=0;j<no_segments;j++)
		{
			clu_data[i][j]=(double*)malloc(no_sample_ps*sizeof(double));
		}
	}

	int size[no_segments];
	int counter=0;
	int class=0;
	for(i=0;i<no_segments;i++)
	{
		for(j=0;j<no_sample_ps;j++)
		{
			clu_data[class][counter][j]=wav_data[noSet[i]][j];

		}
		counter++;
		if(fvSet[i+1]!=fvSet[i])
		{
			size[class]=counter;
			counter=0;
			class++;
		}
	}

	
	/*** add Hamming window ***/
	for(i=0;i<class;i++)
	{
		for(j=0;j<size[i];j++)
		{
		hammingWin(clu_data[i][j],no_sample_ps);
		}
	}
	printf("hamming window\n");
	
	
	/*** DFT get power spectrum ***/
	for(i=0;i<class;i++)
	{
		for(j=0;j<size[i];j++)
		{
			dft(clu_data[i][j],no_sample_ps);
		}
	}
	printf("dft\n");
	
	
	/*** signal normalizarion ***/
	for(i=0;i<class;i++)
	{
		for(j=0;j<size[i];j++)
		{
			normalizationSignal(clu_data[i][j],no_sample_ps);
		}
	}
	printf("nomalizarion\n");

	
	/*** compute mean ***/
	double mean[class][no_sample_ps];
	for(i=0;i<class;i++)
	{
		for(k=0;k<no_sample_ps;k++)
		{
			mean[i][k]=0;
		}
	}
	
	for(i=0;i<class;i++)
	{
		for(k=0;k<no_sample_ps;k++)
		{
			for(j=0;j<size[i];j++)
			{
				mean[i][k]+=clu_data[i][j][k];
			}
			mean[i][k]=mean[i][k]/size[i];
		}
	}
	printf("mean\n");
	
	
	/*** compute segment variance ***/
	double avg_v=0;
	double v=0;
	double d=0;
	//class=class-1;
	for(i=0;i<class;i++)
	{
		v=0;
		for(j=0;j<size[i];j++)
		{
			d=0;
			for(k=0;k<no_sample_ps;k++)
			{
				d+=(clu_data[i][j][k]-mean[i][k])*(clu_data[i][j][k]-mean[i][k]);
			}
			v+=d;
		}
		v=v/size[i];
		avg_v+=v;
	}
	avg_v=avg_v/class;
	printf("avg segment spectrum variance:%f\n",avg_v);
}



/*
	Function description: 
	measure the average simiarity of the segments cluster by computing segment variance
	and generate the re-arranged test wav file
	
	Parameters		: none
	return value	: none
*/
void writeTestWavFile()
{
	int no[no_segments],temp;
	int	i,j,k;
	for(i=0;i<no_segments;i++)
	{
		no[i]=i;
	}

	
	/******* ordering the fv set  *******/
	for(i=0;i<no_segments;i++)
	{
		for(j=i+1;j<no_segments;j++)
		{
			if(fv_set[i]>fv_set[j])
			{
				temp=fv_set[i];
				fv_set[i]=fv_set[j];
				fv_set[j]=temp;
				
				temp=no[i];
				no[i]=no[j];
				no[j]=temp;
			}
		}
	}

	
	/*** compute the segment variance ***/
	//segmentVariance(no, fv_set);

	
	/*** generate re-arranged test wav file ***/
	FILE *file;	
	file=fopen("TestWavFile.wav","wb");
	int space = 10;
	int inter = 0;
	short empty[100000];
	int spaceSize = (N*space+inter*no_segments)*no_sample_ps*size_sample;
	infile.riffsize=infile.riffsize+spaceSize;
	infile.Subchunk2Size=infile.Subchunk2Size+spaceSize;
	fwrite(&infile, sizeof(infile), 1, file);
	for(i=0;i<no_segments;i++)
	{
		fwrite(wav_data[no[i]], size_sample, no_sample_ps, file);
		for(j=0;j<inter;j++)
		{
			fwrite(empty, size_sample, no_sample_ps, file);
		}
		if(fv_set[i+1]!=fv_set[i])
		{
			for(j=0;j<space;j++)
			{
				fwrite(empty, size_sample, no_sample_ps, file);
			}
		}
	}
	fclose(file);
	printf("test wav file generated\n");
}



/*
	Function description: 
	free all the arrays used in this file
	
	Parameters		: none
	return value	: none
*/
void releaseMemory()
{
	int i;
	
	//free array memory location
	for(i=0;i<no_segments;i++)
	{
		free(wav_data[i]);
	}
	free(wav_data);

	for(i=0;i<no_segments;i++)
	{
		free(fv[i]);
	}
	free(fv);
	
	//close wav file
	fclose(wav);
}

#ifndef _WAVE_H_
#define _WAVE_H_
/*
	wav header
*/
typedef struct WAVheader {
         unsigned char riff[4];
         unsigned long riffsize; 		//size of overall file in bytes
         unsigned char wave[8];			//"wavefmt"
         unsigned long hex;				//length of format data listed above
         unsigned short pcm;			//Type of format (1 is PCM)
         unsigned short nochannel;		//Number of Channels
         unsigned long fs;				//Sample Rate
         unsigned long fsch;			//(Sample Rate * BitsPerSample * NumChannels) / 8
         unsigned short ch;				//(BitsPerSample * NumChannels) / 8  The number of bytes for one sample including all channels
         unsigned short bps;			//BitsPerSample
         unsigned char data[4];			//"data"
         unsigned long Subchunk2Size;	//NumSamples * NumChannels * BitsPerSample/8
} Wav;



/*
	feature vector dimension
*/
extern int K;



/*
	number of feature vectors, or the length of training sequence
*/
extern int no_segments;



/*
	Function description: 
	read wav file, print the info of wav header and record all sample data in **wav_data
	
	Parameters: 
	@char wavfile[] : the name of wav file
	@int delta_t 	: the period of each segment, usually between 10 - 20 ms
	
	return value	: none
*/
extern void readWavFile(char wavfile[], int delta_t);



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
extern double** getFeatureVectors(int n_threshold, double correlation_threshold);



/*
	Function description: 
	generate re-arranged test wav file, gather the segments belonging to the same cluster
	
	Parameters		: none
	return value	: none
*/
extern void writeTestWavFile();



/*
	Function description: 
	free all the array
	
	Parameters		: none
	return value	: none
*/
extern void releaseMemory();
#endif
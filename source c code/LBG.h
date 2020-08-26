#ifndef _LBG_H_
#define _LBG_H_
/*
	According to LBG, D(-1) is defined as infinitive to be further compared with D(0) to calculate change rate.
	Since all the fvs have been processed through normalization, thus it is reasonable to just assign number 100
	to represent infinitive, which is totally enough in this project.
*/
#define INF 100



/*
	The set of the id of clusters corresponding to the id of feature vectors.
	The index represent the id of fv, and the value represent represent the id of cluster it belongs to.
	This pointer is later used in producing re-arranged test wav file to figure out the relation between clusters and feature vectors
*/
extern int *fv_set;



/*
	the final number of the representation vectors included in the codebook
*/
extern int N;



/*
	Function description: 
	the complete process of the optimized LBG training
	
	Parameters: 
	@double **fv					: input training sequence of all feature vectors
	@double variance_threshold		: variance threshold for determine the stop point of training
	@double distortion_threshold	: distortion threshold for combining similar clusters
	@double size_threshold			: size_threshold for removing small clusters
	
	return value:
	@double **A		: output codebook after codebook adjustment
*/
extern double** training(double** fv, double variance_threshold, double distortion_threshold, double size_threshold);
#endif
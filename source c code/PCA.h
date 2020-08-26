#ifndef _PCA_H_
#define _PCA_H_
/*
	Function description: 
	get the matrix of feature vectors and codebook
	
	Parameters: 
	@double **v			: input matrix of feature vectors and codebook
	@int num			: the number of all vectors included in the matrix 
	@int dimension		: vector dimension
	
	return value		:none
*/
extern void getMatrix(double **v, int num, int dimension);



/*
	Function description: 
	process zero mean for the matrix **matrix
	
	Parameters		:none
	return value	:none
*/
extern void zeroMean();



/*
	Function description: 
	compute the covariance matrix
	
	Parameters		:none
	return value	:none
*/
extern void covMatrix();



/*
	Function description: 
	compute the first and second eigenvalues
	
	Parameters		:none
	return value	:none
*/
extern void findEigValue();



/*
	Function description: 
	create txt files for recording the coordinates of two-dimensional feature vectors and codebook
	
	Parameters:
	@char* filename		:the name the file to be created
	
	return value:
	@FILE *file 		:the created file
*/
FILE* createFile(char* filename);



/*
	Function description: 
	compute coordinates for the two-dimensional feature vectors and codebook, and record them in txt files
	
	Parameters		:none
	return value	:none
*/
extern void get2dCoordinate();



/*
	Function description: 
	free arrays
	
	Parameters		:none
	return value	:none
*/
extern void release();
#endif
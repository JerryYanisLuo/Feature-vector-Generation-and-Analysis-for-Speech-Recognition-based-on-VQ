#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include "LBG.h"
#include "wave.h"


int					*fv_set;			//the index represent the id of fv, and the value represent the id of cluster it belongs to
double				***S;				//minimum distortion partition
int					n;					//number of training sequence
int					N;					//final training level of LBG
int					M;					//the current training level

double				**split_set;		//splitting set of representation vectors
int					split_num;			//the size of the splitting set
double				**no_split_set;		//non-splitting set of representation vectors
int					no_split_num;		//the size of the non-splitting set


/*
	Function description: 
	Squared-error distortion between two vectors
	
	Parameters: 
	@double *x		: first vector
	@double *y 		: second vector
	@int k 			: dimension of vectors
	
	return value:
	@double d		: resulting distortion
*/
double distortion(double *x, double *y, int k)
{
	double d=0;
	int i;
	for(i=0;i<k;i++)
	{
		d+=(x[i]-y[i])*(x[i]-y[i]);
	}
	return d;
}



/*
	Function description: 
	split the cb to double the number of representation vectors
	
	Parameters: 
	@double **v		: input cb
	@int m 			: number of representation vectors
	@int k 			: dimension of representation vectors
	@int e 			: component of pertubation vector
	
	return value:
	@double **a		: output cb
*/
double** splitting(double **v, int m, int k, double e)
{
	int i,j;
	double** a;
	
	a=(double**)malloc(2*m*sizeof(double*));
	for(i=0;i<2*m;i++)
	{
		a[i]=(double*)malloc(k*sizeof(double));
	}
	
	for(i=0;i<m;i++)
	{
		for(j=0;j<k;j++)
		{
			a[2*i+1][j]=v[i][j]-e;
			a[2*i][j]=v[i][j]+e;
		}
	}
	return a;
}



/*
	Function description: 
	Find optimum partition
	
	Parameters: 
	@double **fv	: feature vectors
	@int n 			: total number of feature vectors
	@double **a		: codebook
	@int m 			: current number of representation vectors
	@int block_size	: codebook dimension
	
	return value:
	@double **a		: output cb
*/
double*** optimumPartition(double** fv, int n, double** a, int m, int block_size )
{
	double d[m], temp;
	int no,i,j,k,l;
	double ***s;
	int count[m];
	fv_set=(int*)malloc(n*sizeof(int));
	s=(double***)malloc(m*sizeof(double**));
	for(i=0;i<m;i++)
	{
		s[i]=(double**)malloc((n+1)*sizeof(double*));
		count[i]=0;
		for(j=0;j<n+1;j++)
		{
			s[i][j]=(double*)malloc(block_size*sizeof(double));
		}
	}
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
		{
			d[j]=distortion(fv[i],a[j],block_size); //square-error distortion
		}
		
		for(j=0;j<m;j++)
		{
			if(j==0)
			{
				temp=d[j];
				no=0;
			}
			else
			{
				if(d[j]<temp)
				{
					no=j;
					temp=d[j];
				}
			}
		}
		fv_set[i]=no;
		
		for(k=0;k<block_size;k++)
		{
			s[no][count[no]][k]=fv[i][k];
		}	
		count[no]++;	
	}
	
	for(i=0;i<m;i++)
	{
		s[i][n][0]=count[i]; //int s[i][n][0], recording the total number of fv in set Si
	}
	
	return s;
}



/*
	Function description: 
	reset the vector, assign 0 to all components
	
	Parameters: 
	@double *a		: input vector
	@int n 			: vector dimension
	
	return value	: none
*/
void clean(int *a, int n)
{
	int i;
	for(i=0;i<n;i++)
	{
		a[i]=0;
	}
}



/*
	Function description: 
	remove a certain vector from a matrix
	
	Parameters: 
	@int index			: the index of the removing vector
	@double **vector	: the input matrix
	@int length			: the number of the vectors inside the matrix
	
	return value		: none
*/
void removeVector(int index, double **vector, int length)
{
	int i,j;
	for(i=index;i<length-1;i++)
	{
		for(j=0;j<K;j++)
		{
			vector[i][j]=vector[i+1][j];
		}
	}
}



/*
	Function description: 
	adjust the pointer fv_set, set certain fvs to an additional cluster which has the largest id   
	
	Parameters: 
	@int i				: the index of the feature vector which should be pushed to the additional cluster
	
	return value		: none
*/
void adjustFvSet(int i)
{
	int j;
	for(j=0;j<no_segments;j++)
	{
		if(fv_set[j]==i)
		{
			fv_set[j]=9999;
		}
	}
}



/*
	Function description: 
	removing small cluster from codebook adjustment, those cluster thought to be the small clusters will be pushed to an additional cluster, 
	and the corresponding representation vector will be removed from the codebook 
	
	Parameters: 
	@double **cb		: codebook
	@int *size			: the set of sizes of all clusters 
	@int threshold		: the size distortion
	
	return value		: none
*/
void removeSmallCluster(double **cb, int *size, int threshold)
{
	int r[M];
	int i,j,m;
	clean(r,M);
	for(i=0;i<M;i++)
	{
		if(size[i]<=threshold)
		{
			r[i]=1;
			adjustFvSet(i);
		}
	}
				
	for(m=M-1;m>=0;m--)
	{
		if(r[m]==1)
		{
			removeVector(m,cb,M);
			M--;
		}
	}	
	printf("reduce small cluster M:%d\n",M);
}



/*
	Function description: 
	combine the splitting set and non-splitting set to get the current codebook 
	
	Parameters: 
	@double **cb			: codebook
	@double **splitSet		: splitting set
	@double **nonSplitSet	: non-splitting set
	@int n					: number of splitted representation vectors
	
	return value		: none
*/
void combineGetCb(double **cb, double **splitSet, double **nonSplitSet, int n)
{
	int i,j;
	for(i=0;i<n*2;i++)
	{
		for(j=0;j<K;j++)
		{
			cb[i][j]=splitSet[i][j];
		}
	}
		
	for(i=n*2;i<M;i++)
	{
		for(j=0;j<K;j++)
		{
			cb[i][j]=nonSplitSet[i-n*2][j];
		}
	}
}



/*
	Function description: 
	get the size set of all clusters
	
	Parameters		: none
	
	return value:
	@int *size		: the array containing the size of each cluster
*/
int* getSize()
{
	int *size;
	int i;
	size=(int*)malloc(M*sizeof(int));
	for(i=0;i<M;i++)
	{
		size[i]=(int)S[i][n][0];
	}
	return size;
}



/*
	Function description: 
	compute the resulting distortion
	
	Parameters: 
	@double ***cluster	: set of clusters
	@int *size			: set of cluster sizes
	@double **cb		: codebook
	
	return value:
	@double d			:resulting distortion
*/
double getResultDistortion(double ***cluster, int *size, double **cb)
{
	int i,j;
	double d=0;
	for(i=0;i<M;i++)
	{
		for(j=0;j<size[i];j++)
		{
			d+=distortion(cluster[i][j], cb[i], K);
		}	
	}
	d=d/n;
	return d;
}



/*
	Function description: 
	update the codebook by updating the cluster centroids
	
	Parameters: 
	@double ***cluster	: set of clusters
	@int *size			: set of cluster sizes
	@double **cb		: codebook
	
	return value		:none
*/
void updateClusterCentroid(double ***cluster, int *size, double **cb)
{
	int i,j,k;
	for(i=0;i<M;i++)
	{
		for(k=0;k<K;k++)
		{
			cb[i][k]=0;
			for(j=0;j<size[i];j++)
			{
				cb[i][k]+=cluster[i][j][k];
			}
			
			if(size[i]==0)
			{
				srand((unsigned int)time(NULL));
				cb[i][k]=(rand()%99999)/100000.0;
			}
			else
			{
				cb[i][k]=cb[i][k]/size[i];	
			}
		}	
	}
}



/*
	Function description: 
	combine similar clusters from codebook adjustment
	
	Parameters: 
	@double **cb		: codebook
	@double threshold	: distortion threshold
	
	return value:
	@int combined		: used as boolean to imply whether there exists a combination or not
*/
int combineSimilarCb(double **cb, double threshold)
{
	double dist;
	int combined=0;
	int i,j,k;
	for(i=0;i<M;i++)
	{
		for(j=i+1;j<M;j++)
		{
			dist=distortion(cb[i], cb[j], K);
			if(dist<threshold)
			{					
				for(k=0;k<K;k++)
				{
					cb[i][k]=(cb[i][k]+cb[j][k])/2;
				}
				removeVector(j,cb,M);
				M--;
				j--;
				combined++;
			}
		}
	}
	return combined;
}



/*
	Function description: 
	compute the cluster variance
	
	Parameters: 
	@double **cluster	: a certain cluster
	@int size			: size of the cluster
	@double **cb		: codebook
	
	return value:
	@double v			:resulting cluster variance
*/
double getClusterVariance(double **cluster, double *cb, int size)
{
	double d,v=0;
	int i,j,k;
	for(j=0;j<size;j++)
	{
		d=0;
		for(k=0;k<K;k++)
		{
			d+=(cluster[j][k]-cb[k])*(cluster[j][k]-cb[k]);
		}
		v+=d;
	}
	v=v/size;
	return v;
}



/*
	Function description: 
	record those clusters demand further splitting by analyze the cluster variance, 
	if the variance of the cluster is below a certain threshold, then the cluster 
	is pushed into splitting set. Otherwise, it is pushed into non-splitting cluster
	The status of each cluster is recorded in the pointer *toSplit after the cluster variance analysis
	The variance of each cluster is recorded in txt file
	
	Parameters: 
	@int *toSplit			: set of cluster status, record whether each cluster is required to be further split of not
	@double threshold		: variance threshold
	@double ***cluster		: set of all clusters
	@double **cb			: codebook
	@int *size				: set of sizes of clusters
	
	return value:
	@double d			:resulting distortion
*/
int clusterVarianceAnalysis(int *toSplit, double threshold, double ***cluster, double **cb, int *size)
{
	double v;
	int flag=1;
	int i;
	clean(toSplit,M);
	for(i=0;i<M;i++)
	{
		v=getClusterVariance(cluster[i], cb[i], size[i]);
		if(v>threshold)
		{
			toSplit[i]=1;
			flag=0;
		}
	}
	return flag;
}



/*
	Function description: 
	initialize splitting sequence and non-splitting sequence based on the cluster status pointer and codebook
	
	Parameters: 
	@int *toSplit		: set of cluster status, record whether each cluster is required to be further split of not
	@double **cb		: codebook
	
	return value		:none
*/
void get_split_nonSplit_sequence(double **cb, int *toSplit)
{
	int i,k;
	split_set=(double**)malloc(M*sizeof(double*));
	no_split_set=(double**)malloc(M*sizeof(double*));
	for(i=0;i<M;i++)
	{
		split_set[i]=(double*)malloc(K*sizeof(double));
		no_split_set[i]=(double*)malloc(K*sizeof(double));
	}
		
	split_num=0;
	no_split_num=0;
	for(i=0;i<M;i++)
	{
		if(toSplit[i])
		{
			for(k=0;k<K;k++)
			{
				split_set[split_num][k]=cb[i][k];
			}
			split_num++;
		}
		else
		{
			for(k=0;k<K;k++)
			{
				no_split_set[no_split_num][k]=cb[i][k];
			}
			no_split_num++;
		}
	}
}



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
double** training(double** fv, double variance_threshold, double distortion_threshold, double size_threshold)
{
	/********************** initialization *************************/
	float			e=0.005;			//the value where distortion should be finally lower than
	double			e_vector=0.005;	
	int				*s_size;
	double			**A;				//codebook sequence
	double			D = 0;				//resulting distortion
	double			d_pre;
	int				i,j,k,m;			//used in for loop
	double			isEnd=0;
	int				record=0;
	
	n = no_segments;
	
	A=(double**)malloc(1024*sizeof(double*));
	for(i=0;i<1024;i++)
	{
		A[i]=(double*)malloc(K*sizeof(double));
	}
	
	split_set=(double**)malloc(1*sizeof(double*));
	for(i=0;i<1;i++)
	{
		split_set[i]=(double*)malloc(K*sizeof(double));
	}
	
	split_num=1;
	no_split_num=0;

	for(j=0;j<K;j++)
	{
		split_set[0][j]=0;
		for(i=0;i<n;i++)
		{
			split_set[0][j]+=fv[i][j];
		}
		split_set[0][j]=split_set[0][j]/n;
	}
	
	
	FILE *file_distortion;
	file_distortion=fopen("distortion/distortion.txt","w");
	
	for(M=1;;)
	{
		/********************** splitting *************************/
		split_set=splitting(split_set, split_num, K, e_vector);
		M=M+split_num;
		combineGetCb(A, split_set, no_split_set, split_num);
		d_pre=INF;
		printf("M:%d\n",M);
		while(1)
		{
			/************* compute distortion to find optimum partition **********/
			S=optimumPartition(fv, n, A, M, K);
			s_size=getSize();
			
			/******************** compute resulting distortion *******************/
			D=getResultDistortion(S, s_size, A);
			printf("distortion %f\n",D);

			if((d_pre-D)/D<=e)
			{
				fprintf(file_distortion,"%f ",D);
				break;
			}
			d_pre=D;
			D=0;
			
			/*********************** update cluster centers **********************/
			updateClusterCentroid(S, s_size, A);
			
			free(s_size);
			free(fv_set);
		}

		
		/************* Cluster variance analysis *************/
		int	toSplit[M];
		isEnd=clusterVarianceAnalysis(toSplit,variance_threshold,S,A,s_size);
		if(!isEnd)
		{
			get_split_nonSplit_sequence(A, toSplit);
		}
		else
		{
			printf("Stop training M:%d\n",M);
			
			/*************************** Codebook adjustment ****************************/
			printf("codebook adjustment\n");			
			
			/************* Combine similar clusters *************/
			int combined=0;
			while(1)
			{
				/************ combine similar codeword **********/
				combined=combineSimilarCb(A, distortion_threshold);
				if(combined)
				{
					d_pre=INF;
					while(1)
					{
						/************ Find optimum partition **********/
						S=optimumPartition(fv, n, A, M, K);
						s_size=getSize();
						/********** compute result distortion *********/
						D=getResultDistortion(S, s_size, A);
						if((d_pre-D)/D<=e)
						{
							break;
						}
						d_pre=D;
						D=0;
						/************ update cluster centers **********/
						updateClusterCentroid(S, s_size, A);
						free(s_size);
						free(fv_set);
					}
				}
				else
				{
					break;
				}
			}
			printf("combine similar clusters M:%d\n",M);
			
			
			/************* record each cluster variance *************/
			double v;
			double sum=0;
			int i;
			FILE *file_cv;
			file_cv=fopen("cluster variance/cluster_variance.txt","w");
			clean(toSplit,M);
			for(i=0;i<M;i++)
			{
				if(s_size[i]>size_threshold)
				{
					v=getClusterVariance(S[i], A[i], s_size[i]);
					fprintf(file_cv,"%f ",v);
					sum+=v;
				}
			}
			fclose(file_cv);
			
			
			/********** Remove small clusters **********/
			removeSmallCluster(A,s_size,size_threshold);
			printf("avg cluster variance: %f\n",sum/M);

			
			/************* free arrays *************/
			for(i=0;i<M;i++)
			{
				for(j=0;j<n+1;j++)
				{
					free(S[i][j]);
				}
				free(S[i]);
			}
			free(S);
			
			for(i=0;i<split_num;i++)
			{
				free(split_set[i]);
			}
			free(split_set);
			
			for(i=0;i<no_split_num;i++)
			{
				free(no_split_set[i]);
			}
			free(no_split_set);
			
			free(s_size);
			
			N=M;
			break;
		}
	}
	fclose(file_distortion);
	printf("---------------------------------------------\n");
	printf("codebook generated\n");

	
	/*** record codebook ***/
	FILE *file_cb;
	file_cb=fopen("codebook/codebook.txt","w");
	for(i=0;i<M;i++)
	{
		for(k=0;k<K;k++)
		{
			fprintf(file_cb,"%f ",A[i][k]);
		}
		fprintf(file_cb,"\n");
	}
	return A;
}


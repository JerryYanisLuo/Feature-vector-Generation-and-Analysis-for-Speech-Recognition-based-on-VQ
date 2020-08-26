#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "PCA.h"
#include "wave.h"
#include "lbg.h"


double **matrix;	//input matrix
double **c;			//covariance matrix
double **v;			//eigenvectors
int n,k;


/*
	Function description: 
	get the matrix of feature vectors and codebook
	
	Parameters: 
	@double **v			: input matrix of feature vectors and codebook
	@int num			: the number of all vectors included in the matrix 
	@int dimension		: vector dimension
	
	return value		:none
*/
void getMatrix(double **v, int num, int dimension)
{
	int i,j;	
	n=num;
	k=dimension;
	matrix=(double**)malloc(n*sizeof(double*));
	for(i=0;i<n;i++)
	{
		matrix[i]=(double*)malloc(k*sizeof(double));
	}
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<k;j++)
		{
			matrix[i][j]=v[i][j];
		}
	}
}



/*
	Function description: 
	process zero mean for the matrix **matrix
	
	Parameters		:none
	return value	:none
*/
void zeroMean()
{
	double mean;
	int i,j;
	for(j=0;j<k;j++)
	{
		mean=0;
		for(i=0;i<n;i++)
		{
			mean+=matrix[i][j];
		}
		mean=mean/n;
		
		for(i=0;i<n;i++)
		{
			matrix[i][j]=matrix[i][j]-mean;
		}
	}
}



/*
	Function description: 
	compute the covariance matrix
	
	Parameters		:none
	return value	:none
*/
void covMatrix()
{
	int i,j,l;
	c=(double**)malloc(k*sizeof(double*));
	for(i=0;i<k;i++)
	{
		c[i]=(double*)malloc(k*sizeof(double));
	}
	
	for(i=0;i<k;i++)
	{
		for(j=0;j<k;j++)
		{
			c[i][j]=0;
			for(l=0;l<n;l++)
			{
				c[i][j]+=matrix[l][i]*matrix[l][j];
			}
			c[i][j]=c[i][j]/(n-1);
		}
	}
}



/*
	Function description: 
	Find the first and second eigenvalues and derive the corresponding eigenvectors
	In this function, the power method is used to approximate the extremal eigenvalues of the matrix
	
	Parameters		:none
	return value	:none
*/
void findEigValue()
{
	double x[k], z[k], y[k], d, temp, eigen;
	int i,j,p;
	v=(double**)malloc(2*sizeof(double*));
	for(i=0;i<2;i++)
	{
		v[i]=(double*)malloc(k*sizeof(double));
	}
	
	for(p=0;p<2;p++)
	{
		//initial first random vector
		for(i=0;i<k;i++)
		{
			x[i]=0;
		}
		x[0]=1;
		
		while(1)
		{
			temp=0;
			//normalization
			for(i=0;i<k;i++)
			{
				if(x[i]>temp)
				{
					temp=x[i];
				}
			}
			
			for(i=0;i<k;i++)
			{
				z[i]=x[i]/temp;
			}
			
			//literation
			for(i=0;i<k;i++)
			{
				y[i]=0;
			}
			
			for(i=0;i<k;i++)
			{
				for(j=0;j<k;j++)
				{
					y[i]+=c[i][j]*z[j];
				}
			}
		
			//compute distortion
			d=0;
			for(i=0;i<k;i++)
			{
				d+=(y[i]-x[i])*(y[i]-x[i]);
			}
			
			d= sqrt(d);
			
			temp=0;
			if(d<0.00001)
			{
				//find eig
				for(i=0;i<k;i++)
				{
					if(y[i]>temp)
					{
						temp=y[i];
					}
				}
				eigen=temp;
				printf("eigen: %f\n",eigen);
				
				for(i=0;i<k;i++)
				{
					y[i]=y[i]/temp;
				}
				
				temp=0;
				for(i=0;i<k;i++)
				{
					temp+=y[i]*y[i];
				}
				
				temp= sqrt(temp);
				
				for(i=0;i<k;i++)
				{
					z[i]=y[i]/temp;
					v[p][i]=z[i];
				}
				break;
			}
			
			for(i=0;i<k;i++)
			{
				x[i]=y[i];
			}
		}
		
		for(i=0;i<k;i++)
		{
			for(j=0;j<k;j++)
			{
				c[i][j]=c[i][j]-eigen*z[i]*z[j];
			}
		}
	}
}



/*
	Function description: 
	create txt files for recording the coordinates of two-dimensional feature vectors and codebook
	
	Parameters:
	@char* filename		:the name the file to be created
	
	return value:
	@FILE *file 		:the created file
*/
FILE* createFile(char* filename)
{
	FILE *file=fopen(filename,"w");	
	return file;
}



/*
	Function description: 
	compute coordinates for the two-dimensional feature vectors and codebook, and record them in txt files
	
	Parameters		:none
	return value	:none
*/
void get2dCoordinate()
{
	int i,j,l;
	double coordi[n][2];
	FILE *file;
	
	for(i=0;i<n;i++)
	{
		for(l=0;l<2;l++)
		{
			coordi[i][l]=0;
			for(j=0;j<k;j++)
			{
				coordi[i][l]+=matrix[i][j]*v[l][j];
			}
		}
		
	}
	
	for(l=0;l<2;l++)
	{
		if(l==0)
		{
			file=createFile("PCA2d plot/fvx.txt");
		}
		else if(l==1)
		{
			file=createFile("PCA2d plot/fvy.txt");
		}

		for(i=0;i<n;i++)
		{
			fprintf(file,"%f ",coordi[i][l]);
			if(i==(no_segments-1))
			{
				fclose(file);
				printf("file fv%d complete\n",l);
				if(l==0)
				{
					file=createFile("PCA2d plot/cbx.txt");
				}
				else if(l==1)
				{
					file=createFile("PCA2d plot/cby.txt");
				}
			}
		}
		fclose(file);
		printf("file cb%d complete\n",l);
	}
}



/*
	Function description: 
	free arrays
	
	Parameters		:none
	return value	:none
*/
void release()
{
	int i;
	for(i=0;i<n;i++)
	{
		free(matrix[i]);
	}
	free(matrix);
	
	for(i=0;i<k;i++)
	{
		free(c[i]);
	}
	free(c);
	
	for(i=0;i<2;i++)
	{
		free(v[i]);
	}
	free(v);
}



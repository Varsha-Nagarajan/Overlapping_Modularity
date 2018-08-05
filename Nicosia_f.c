/*=========================================
Author : Varsha Nagarajan

This provides an implementation of Overlapping Modularity as proposed in the paper 
"Nicosia V, Mangioni G, Carchiolo V and Malgeri M 2009 Extending the definition of
modularity to directed graphs with overlapping communities J. Stat. Mech. P03024"
The same "two–dimensional logistic function" as mentioned in the paper has been used
and the values for p in equation 23 is taken as suggested in 
Nicosia V, Mangioni G, Carchiolo V and Malgeri M 2008 Extending modularity definition for
directed graphs with overlapping communities arXiv:0801.1647

==========================================*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

int **array;
int **array1;
int *community;
float **v_comm;
double func( float a,float b);
double funcb( int comm , int node);
int vertices,m;

int main(int argc , char ** argv)
{
	char adjlist_file[50]="", coverFile[50]="";
	int i,t,deg,node,comm=1,num_of_clusters,k,k1,j,tt,flag=0,a,b;
	
	strcat(adjlist_file,argv[1]);  // adjacency list
	strcat(coverFile,argv[2]); // final cover
	
	float j1;
	double part11=0.0,part22=0.0,alpha1=0.0,alpha2=0.0,beta=0.0,part1=0.0,part2=0.0,final=0.0,final1=0.0,part21=0.0;
	
	FILE *fp;
	fp=fopen(adjlist_file,"r");
	if(fp==NULL)
	{
		printf("%s does not exist\n", adjlist_file);
		exit(0);
	}
	fscanf(fp,"%d %d ",&vertices,&edges);
	array=(int**)calloc(vertices+1,sizeof(int*));
	while(fscanf(fp,"%d ",&node)!=EOF)
    {
		fscanf(fp,"%d ",&deg);
		array[node]=(int*)calloc(deg+1,sizeof(int));
		array[node][0]=deg;
		for(i=1;i<=deg;i++)		
        {
            fscanf(fp,"%d ",&t);
            array[node][i]=t;
        }
	}
	int countl=0,test=0;
	fclose(fp);
	int *counter=(int *)calloc(vertices+1,sizeof(int));
	fp=fopen(coverFile,"r"); // cover file
	fscanf(fp,"%d ",&num_of_clusters); // number of covers
	v_comm=(float**)calloc(vertices+1,sizeof(float*));
	// vertex X cover matrix is being constructed.
	// if v_comm[1][2] == 1, it means vertex/node label 1 is part of community number 2
    for(i=0;i<=vertices;i++) 
        v_comm[i]=(float*)calloc(num_of_clusters+1,sizeof(float));
	community=(int*)calloc(num_of_clusters+1,sizeof(int));
	while(fscanf(fp,"%d",&i)!=EOF)
	{
		if(i!=-1)
		{
			v_comm[i][comm]=1.0; // indicates that vertex label 'i' belongs to community "comm"
			community[comm]++; // counts number of nodes in a given cluster
			counter[i]++; // counts number of communities a node belongs to
		}
		else
		{
			comm++; // increments the community count. At end, this value will store the number of clusters
		}
	}
	fclose(fp);

	for(i=1;i<=vertices;i++)
	{
		for(j=1;j<=num_of_clusters;j++) 
		{
			if(i==1)			
				v_comm[j][0]=community[j]; // this holds the cluster size of each cluster.
			if(v_comm[i][j]==1.0)
				v_comm[i][j]=(float)1/(counter[i]);
		}
	}
	
	free(counter);
	free(community);
	
	int *comm1;
	double *mark=(double*)calloc(vertices+1,sizeof(double));
	double partt=0.0,mark1=1.0,mark2=0.0;
	double tot_beta=1.0,tot=0.0,tot1=0.0;
	int kij=0;
	float alpha_val,alpha2_val;int jk=0;
	int index;int temp=0;double betav=0.0;int alpa=0;
	for(k=1;k<=num_of_clusters;k++)
	{
		temp=0;	
		//printf("Community %d\n",k);
		part1=0.0,part2=0.0;	
		for(j=0;j<=vertices;j++)
		{
			mark[j]=-1;
		}
		jk=0;
		comm1=(int*)calloc((int)v_comm[k][0],sizeof(int)); // array of size equal to number of nodes in that community
		
		for(j=1;j<=vertices;j++)
		{
			if(v_comm[j][k]!=0)
			{
				comm1[jk]=j; //holds index of node that belongs to this community
				jk++;
			}
			alpha2_val=v_comm[j][k];
			if(mark[j]==-1.0)
			{		
				if(alpha2_val!=0)
					mark[j]=funcb(k,j);
				else 
					mark[j]=0;
				for(index=1;index<=vertices;index++)
				{
					if(v_comm[index][k]==alpha2_val)
						mark[index]=mark[j];
				}
			}
		}

		for(i=1;i<=vertices;i++)
		{
			for(j=1;j<=array[i][0];j++)
			{
				if(i<array[i][j])
				{
					part1+=(func(v_comm[i][k],v_comm[array[i][j]][k]));
				}
			}
		}
		for(i=0;i<(int)(v_comm[k][0]);i++)
		{
			for(j=0;j<(int)(v_comm[k][0]);j++)
			{
				part2+=((double)array[comm1[i]][0]*(double)array[comm1[j]][0]*mark[comm1[i]]*mark[comm1[j]]);
			}
		}
		
		free(comm1);
		part2=part2/4;
		part2=part2/(double)m;
		final+=((part1/(double)(m))-(part2/(double)m));
		//printf("comm --------> %d %lf\n",k,(part1/(double)(m))-(part2/(double)m));
		//printf("Community %d / %d : %lf\n",k,num,final);
	}
								
	free(mark);
	printf("Modularity is %lf :\n",final);
	
	// Freeing the resources
	for(i=0;i<=vertices;i++)
        free(array[i]);
    free(array);
	
	for(i=0;i<=vertices;i++)
        free(v_comm[i]);
    free(v_comm);

}

double func( float a, float b)
{
	double f1=(60*a)-30;
	double f2=(60*b)-30;
	f1=exp(-(f1));
	f2=exp(-(f2));
	f1+=1;
	f2+=1;
	double l=(double)1/(f1*f2);
	return(l);
}

double funcb( int comm , int node)
{
	int j;
	double beta=0.0;
	float alpha2,alpha1;
	alpha1=v_comm[node][comm];
	for(j=1;j<=vertices;j++)
	{
		alpha2=v_comm[j][comm];
		beta+=func(alpha1,alpha2);
		
	}
	beta=beta/(double)(vertices);

	return beta;
}

#include<stdio.h>

#define numKernels 512
#define hashTableWidth 12288
#define searchKernels 1024
#define numHashPerThread 24
#define numThreadsPerBlock 1


__device__ void swapTriplets(int *d_x,int *d_y,int *d_z,int *d_indexOrder, int index1,int index2)
{
	int temp1,temp2,temp3;
	int temp4;
	temp1=d_x[index1];
	temp2=d_y[index1];
	temp3=d_z[index1];
	temp4=d_indexOrder[index1];
	d_x[index1]=d_x[index2];
	d_y[index1]=d_y[index2];
	d_z[index1]=d_z[index2];
	d_indexOrder[index1]=d_indexOrder[index2];
	d_x[index2]=temp1;
	d_y[index2]=temp2;
	d_z[index2]=temp3;
	d_indexOrder[index2]=temp4;
}


__device__ int compare(int x1,int y1, int z1, int x2, int y2, int z2)
{
	if((x1==x2) && (y1==y2) && (z1==z2))
		return 0;
	if((x1>x2) ||((x1==x2) && (y1>y2)) || ((x1==x2) && (y1==y2) && (z1>z2)))
		return 1;
	return -1;
}

__global__ void createPartition(int count,int *d_x,int *d_y,int *d_z,int *d_indexOrder,int *d_partition_begin,int *d_partition_last)
{
	int partitionIndex=blockIdx.x*blockDim.x+threadIdx.x;
	int begin,last,i,j;
	int x,y,z;
	int indexOrder;
	begin=d_partition_begin[partitionIndex];
	last=d_partition_last[partitionIndex];
	if((begin>=0) && (last>=0))
	{
		if(begin<last)
		{
			swapTriplets(d_x,d_y,d_z, d_indexOrder
,last,((begin+last)/2));
			x=d_x[last];
			y=d_y[last];
			z=d_z[last];
			indexOrder=d_indexOrder[last];
			i = (begin - 1);
	
			for ( j = begin; j <= last- 1; j++)
			{
				if (compare(d_x[j],d_y[j],d_z[j],x,y,z)<=0)
				{
    					i++;
    					swapTriplets (d_x,d_y,d_z,d_indexOrder,i, j);
				}
			}
			swapTriplets (d_x,d_y,d_z,d_indexOrder,(i+1),(last));

		}
	}
	__syncthreads();
	if((begin>=0) && (last>=0))
	{
		if(begin<last)
		{
			d_partition_begin[2*partitionIndex]=begin;
			d_partition_last[2*partitionIndex]=i;
			d_partition_begin[2*partitionIndex+1]=i+1;
			d_partition_last[2*partitionIndex+1]=last;
		}
		else
		{
			d_partition_begin[2*partitionIndex]=begin;
			d_partition_last[2*partitionIndex]=last;
			d_partition_begin[2*partitionIndex+1]=-1;
			d_partition_last[2*partitionIndex+1]=-1;
		}
	}
	else
	{
		d_partition_begin[2*partitionIndex]=-1;
		d_partition_last[2*partitionIndex]=-1;
		d_partition_begin[2*partitionIndex+1]=-1;
		d_partition_last[2*partitionIndex+1]=-1;
	}
}



__global__ void countPartitions(int *d_partition_begin,int *d_partition_last,int *d_numPartitions)
{
	int i,j;
	int temp_begin,temp_last;
	j=-1;
	for(i=0;i<numKernels;i++)
	{
		if((d_partition_begin[i]!=-1)&&(d_partition_last[i]!=-1))
		{
			j++;
			temp_begin=d_partition_begin[i];
			temp_last=d_partition_last[i];
			d_partition_begin[i]=d_partition_begin[j];
			d_partition_last[i]=d_partition_last[j];
			d_partition_begin[j]=temp_begin;
			d_partition_last[j]=temp_last;

		}
	}
	d_numPartitions[0]=j+1;

}

__global__ void initializePartitionTable(int *d_partition_begin, int *d_partition_last)
{
	int index=blockIdx.x*blockDim.x+threadIdx.x;
	if(index<numKernels)
	{
		d_partition_begin[index]=-1;
		d_partition_last[index]=-1;
	}
}


__device__ int binarySearch(int *d_x, int *d_y, int *d_z,int begin, int offset,int hashvalue)
{

	int mid,i;
	if(begin<=offset)
	{
		mid=(begin+offset)/2;
		if(((d_x[mid]+d_y[mid]+d_z[mid])%numHashPerThread)==hashvalue)
		{
			for(i=mid+1;i<=offset;i++)
			{
				if(((d_x[i]+d_y[i]+d_z[i])%numHashPerThread)>hashvalue)
					break;
			}
			return i;
		}
		else if((((d_x[mid]+d_y[mid]+d_z[mid])%numHashPerThread)<hashvalue) && (((mid<offset)&&(((d_x[mid+1]+d_y[mid+1]+d_z[mid+1])%numHashPerThread)>hashvalue))||(mid==offset)))
		{
			return mid;
		}
		else if (((d_x[mid]+d_y[mid]+d_z[mid])%numHashPerThread)>hashvalue)
			return binarySearch(d_x,d_y,d_z,begin,mid-1,hashvalue);
		else
			return binarySearch(d_x,d_y,d_z,mid+1, offset,hashvalue);
	}
	return -1;



}

__device__ void shiftPos(int *d_x, int *d_y, int *d_z,int *d_indexOrder, int offset)
{
	int index=blockIdx.x*blockDim.x+threadIdx.x;
	int x,y,z;
	int indexOrder;
	x=d_x[offset-index];
	y=d_y[offset-index];
	z=d_z[offset-index];
	indexOrder=d_indexOrder[offset-index];
	__syncthreads();
	d_x[offset-index+1]=x;
	d_y[offset-index+1]=y;
	d_z[offset-index+1]=z;
	d_indexOrder[offset-index+1]=indexOrder;
}


__global__ void setIndexOrder(int *d_indexOrder, int count)
{
	int index=blockIdx.x*blockDim.x+threadIdx.x;
	int i;
	for(i=index*((count+numKernels-1)/numKernels);i<(index+1)*((count+numKernels-1)/numKernels);i++)
		if(i<count)
			d_indexOrder[i]=i;
}


__global__ void sortAndBFS(int *d_x,int *d_y, int *d_z,int *d_indexOrder, int *d_partition_begin, int *d_partition_last,int d_numPartitions, int count, int *d_labels,int *d_queue,int *d_front,int *d_rear,int *d_numGroups, int *d_neighbours)
{
	int begin, last,i,j,sortPos,x,y,z,numberOfLabels, indexOrder,index, front, rear,partitionIndex;
	
	int hashvalue,dx[6],dy[6],dz[6],tempHashValue,k;

	partitionIndex=blockIdx.x*blockDim.x+threadIdx.x;
	if(partitionIndex>=d_numPartitions)
		return;
	__shared__ int hashGlobalMemory[hashTableWidth];

	begin=d_partition_begin[partitionIndex];
	last=d_partition_last[partitionIndex];
	for(i=begin;i<=last;i++)
	{
		hashvalue=(d_x[i]+d_y[i]+d_z[i])%numHashPerThread;
		sortPos=binarySearch(d_x,d_y, d_z,begin,(i-1),hashvalue);
		x=d_x[i];
		y=d_y[i];
		z=d_z[i];
		indexOrder=d_indexOrder[i];
		for(j=i-1;j-numThreadsPerBlock>=sortPos;j=j-numThreadsPerBlock)
		{
			shiftPos(d_x,d_y,d_z,d_indexOrder,j);
		}
		for(;j>=sortPos;j--)
		{
			shiftPos(d_x,d_y,d_z,d_indexOrder,j);
		}
		d_x[sortPos]=x;
		d_y[sortPos]=y;
		d_z[sortPos]=z;
		d_indexOrder[sortPos]=indexOrder;

	}
	hashvalue=-1;
	for(i=partitionIndex*numHashPerThread;i<(partitionIndex+1)*numHashPerThread;i++)
	{
		hashGlobalMemory[i]=-1;
	}

	for(i=begin;i<=last;i++)
	{
		d_labels[i]=-1;
		tempHashValue=(d_x[i]+d_y[i]+d_z[i])%numHashPerThread;
		if(tempHashValue!=hashvalue)
		{
			hashGlobalMemory[partitionIndex*numHashPerThread+tempHashValue]=i;
			hashvalue=i;
		}

	}




	if(partitionIndex==0)
	{
		d_front[0]=-1;
		d_rear[0]=-1;
		d_numGroups[0]=0;
	}



	__syncthreads();

	for(i=0;i<count;i++)
	{
		if(partitionIndex==0)
		{
			if(d_labels[d_indexOrder[i]]==-1)
			{
				d_front[0]=0;
				d_rear[0]=0;
				d_numGroups[0]++;
				d_queue[d_rear[0]]=i;
			}
	
		}

		while(d_front[0]!=-1 || rear>front)
		{
			front=d_front[0];
			rear=d_rear[0];
			for(j=0;j<6;j++)
				d_neighbours[j]=-1;
			dx[0]=d_x[d_queue[front]]-1;
			dy[0]=d_y[d_queue[front]];
			dz[0]=d_z[d_queue[front]];
			dx[1]=d_x[d_queue[front]];
			dy[1]=d_y[d_queue[front]]-1;
			dz[1]=d_z[d_queue[front]];
			dx[2]=d_x[d_queue[front]]+1;
			dy[2]=d_y[d_queue[front]];
			dz[2]=d_z[d_queue[front]];
			dx[3]=d_x[d_queue[front]];
			dy[3]=d_y[d_queue[front]]+1;
			dz[3]=d_z[d_queue[front]];
			dx[4]=d_x[d_queue[front]];
			dy[4]=d_y[d_queue[front]];
			dz[4]=d_z[d_queue[front]]+1;
			dx[5]=d_x[d_queue[front]];
			dy[5]=d_y[d_queue[front]];
			dz[5]=d_z[d_queue[front]]-1;

			for(j=0;j<6;j++)
			{
				if((compare(dx[j],dy[j],dz[j],d_x[d_partition_begin[partitionIndex]],d_y[d_partition_begin[partitionIndex]],d_z[d_partition_begin[partitionIndex]])>=0)&&(compare(dx[j],dy[j],dz[j],d_x[d_partition_begin[partitionIndex]],d_y[d_partition_begin[partitionIndex]],d_z[d_partition_begin[partitionIndex]])<=0))
				{
					if(hashGlobalMemory[partitionIndex*numHashPerThread+((dx[j]+dy[j]+dz[j])%numHashPerThread)]!=-1)
					{
						for(k=hashGlobalMemory[partitionIndex*numHashPerThread+((dx[j]+dy[j]+dz[j])%numHashPerThread)];((d_x[k]+d_y[k]+d_z[k])%numHashPerThread)==((dx[j]+dy[j]+dz[j])%numHashPerThread);k++)
						{
							if(compare(d_x[k],d_y[k],d_z[k],dx[j],dy[j],dz[j])==0)
							{
								if(d_labels[d_indexOrder[k]]==-1)
								{
									d_neighbours[j]=k;
									d_labels[d_indexOrder[k]]=d_numGroups[0];
								}
								break;
							}
							if(k==count-1)
								break;

						}
					}
				}
			}
			__syncthreads();
			if(partitionIndex==0)
			{
				for(k=0;k<6;k++)
				{
					if(d_neighbours[k]!=-1)
					{
						d_rear[0]++;
						d_queue[d_rear[0]]=d_neighbours[k];
					}
				}
				d_front[0]++;
			}
			__syncthreads();
		}
	}

}


int main()
{
	int count,i,*d_front,*d_rear,h_front,h_rear,*h_labels,*d_labels,*h_partition_begin,*d_partition_begin, *orderedGroups;
	int *h_partition_last,*d_partition_last,*d_queue, *h_numGroups, *d_numGroups,*d_indexOrder,*d_neighbours;
	int nx,ny,nz,x,y,z,*h_x,*h_y,*h_z,*d_x,*d_y,*d_z,*h_numPartitions,*d_numPartitions,groupLabelJockey;
	FILE *fp, *ofp;

	fp=fopen("data.txt","r");
	count=0;
	fscanf(fp,"Nx=%d Ny=%d Nz=%d",&nx,&ny,&nz);
	printf("%d %d %d\n",nx,ny,nz);
	while(feof(fp)==0)
	{
		fscanf(fp,"%d %d %d\n",&x,&y,&z);
		printf("%d %d %d\n",x,y,z);
		count++;
	}
	printf("%d\n",count);
	fclose(fp);
	h_x=(int*)malloc(sizeof(int)*count);
	h_y=(int*)malloc(sizeof(int)*count);
	h_z=(int*)malloc(sizeof(int)*count);
	h_partition_begin=(int*)malloc(sizeof(int));
	h_partition_last=(int*)malloc(sizeof(int));
	h_labels=(int*)malloc(sizeof(int)*count);
	h_numPartitions=(int*)malloc(sizeof(int));
	h_numGroups=(int*)malloc(sizeof(int));

	cudaMalloc(&d_x,count*sizeof(int));
	cudaMalloc(&d_y,count*sizeof(int));
	cudaMalloc(&d_z,count*sizeof(int));
	cudaMalloc(&d_partition_begin,numKernels*sizeof(int));
	cudaMalloc(&d_partition_last,numKernels*sizeof(int));
	cudaMalloc(&d_queue,sizeof(int)*count);
	cudaMalloc(&d_labels,sizeof(int)*count);
	cudaMalloc(&d_front,sizeof(int));
	cudaMalloc(&d_rear,sizeof(int));
	cudaMalloc(&d_labels,sizeof(int)*count);
	cudaMalloc(&d_numPartitions,sizeof(int));
	cudaMalloc(&d_numGroups,sizeof(int));
	cudaMalloc(&d_indexOrder, sizeof(int)*count);
	cudaMalloc(&d_neighbours,sizeof(int)*6); /*0:left, 1: behind, 2:right,3: front, 4: top, 5:bottom*/


	fp=fopen("data.txt","r");
	fscanf(fp,"Nx=%d Ny=%d Nz=%d",&nx,&ny,&nz);
	printf("%d %d %d\n",nx,ny,nz);
	for(i=0;i<count;i++)
	fscanf(fp,"%d %d %d\n",&h_x[i],&h_y[i],&h_z[i]);
	
	fclose(fp);

	h_partition_begin[0]=0;
	h_partition_last[0]=count-1;

	cudaMemcpy(d_x,h_x,count*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d_y,h_y,count*sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d_z,h_z,count*sizeof(int),cudaMemcpyHostToDevice);
	//cudaMemcpy(d_front,h_front,sizeof(int),cudaMemcpyHostToDevice);
	//cudaMemcpy(d_rear,h_rear,sizeof(int),cudaMemcpyHostToDevice);
	//cudaMemcpy(d_rear,h_rear,sizeof(int),cudaMemcpyHostToDevice);


	setIndexOrder<<<numKernels,1>>>(d_indexOrder,count);

	initializePartitionTable<<<numKernels,1>>>(d_partition_begin,d_partition_last);
	
	cudaMemcpy(d_partition_begin,h_partition_begin,sizeof(int),cudaMemcpyHostToDevice);
	cudaMemcpy(d_partition_last,h_partition_last,sizeof(int),cudaMemcpyHostToDevice);
	for(i=1;i<=(numKernels/2);i=i*2)
	{

		createPartition<<<1,i>>>(count,d_x,d_y,d_z,d_indexOrder,d_partition_begin,d_partition_last);

	}

	countPartitions<<<1,1>>>(d_partition_begin,d_partition_last,d_numPartitions);

	cudaMemcpy(h_numPartitions,d_numPartitions,sizeof(int),cudaMemcpyDeviceToHost);

	sortAndBFS<<<1,numKernels>>>(d_x,d_y,d_z,d_indexOrder,d_partition_begin,d_partition_last,h_numPartitions[0],count,d_labels,d_queue,d_front,d_rear,d_numGroups,d_neighbours);

	cudaMemcpy(h_labels,d_labels,count*sizeof(int),cudaMemcpyDeviceToHost);
	cudaMemcpy(h_numGroups,d_numGroups,sizeof(int),cudaMemcpyDeviceToHost);
	orderedGroups=(int*)malloc(sizeof(int)*(h_numGroups[0]+1));
	for(i=1;i<=h_numGroups[0];i++)
		orderedGroups[i]=-1;
	groupLabelJockey=0;
	for(i=0;i<count;i++)
	if(orderedGroups[h_labels[i]]==-1)
	{
		groupLabelJockey++;
		orderedGroups[h_labels[i]]=groupLabelJockey;
	}
	ofp=fopen("result.txt","w");
	fprintf(ofp,"Nx=%d Ny=%d  Nz=%d Cluster=%d\nX Y Z id\n",nx,ny,nz,h_numGroups[0]);
	for(i=0;i<count;i++)
		fprintf(ofp,"%d %d %d %d\n",h_x[i],h_y[i],h_z[i],orderedGroups[h_labels[i]]);
	fclose(ofp);

	free(h_labels);
	free(h_x);
	free(h_y);	
	free(h_z);
	free(h_partition_begin);
	free(h_partition_last);
	
	cudaFree(d_front);
	cudaFree(d_rear);
	cudaFree(d_labels);
	cudaFree(d_partition_begin);
	cudaFree(d_partition_last);
	cudaFree(d_x);
	cudaFree(d_y);
	cudaFree(d_z);
	cudaFree(d_queue);
	cudaFree(d_numGroups);
	cudaFree(d_indexOrder);
	cudaFree(d_neighbours);
	
	return 0;
}

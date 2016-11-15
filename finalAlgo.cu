#include<stdio.h>
#include<math.h>
#include<limits.h>
#define numKernels 512
#define hashTableWidth 2048
#define searchKernels 1024
#define numThreadsPerBatch 32


__device__ void swapTriplets(int *d_x,int *d_y,int *d_z,long long d_indexOrder, long long index1,long long index2)
{
	int temp1,temp2,temp3;
	long long temp4;
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


int compare(int x1,int y1, int z1, int x2, int y2, int z2)
{
	if((x1==x2) && (y1==y2) && (z1==z2))
		return 0;
	if((x1>x2) ||((x1==x2) && (y1>y2)) || ((x1==x2) && (y1==y2) && (z1>z2)))
		return 1;
	return -1;
}

__global__ void createPartition(long long count,int *d_x,int *d_y,int *d_z,long long *d_indexOrder,int *d_partition_begin,int *d_partition_last)
{
	int partitionIndex=blockIdx.x*blockDim.x+threadIdx.x;
	long long begin,last,i,j;
	int x,y,z;
	long long indexOrder
	begin=d_partition_begin[partitionIndex];
	last=d_partition_last[partitionIndex];
	if((begin>=0) && (last>=0))
	{
		if(begin<last)
		{
			swapTriplets(dx,dy,dz,last,((begin+last)/2));
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
	_syncthreads();
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



__global__ void countPartitions(long long *d_partition_begin,long long *d_partition_last,int *d_numPartitions)
{
	int i,j;
	long long temp_begin,temp_last;
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

__global__ void initializePartitionTable(long long *d_partition_begin, long long *d_partition_last)
{
	int index=blockIdx.x*blockDim.x+threadIdx.x;
	if(index<numKernels)
	{
		d_partition_begin[index]=-1;
		d_partition_last[index]=-1;
	}
}


__device__ binarySearch(int *d_x, int *d_y, int *d_z,long long begin, long long offset,int hashvalue)
{
	long long mid;
	if(begin<=offset)
	{
		mid=(begin+offset)/2;
		if(((d_x[mid]+d_y[mid]+d_z[mid])%hashTableWidth)==hashvalue)
		{
			for(i=mid+1;i<=offset;i++)
			{
				if(((d_x[i]+d_y[i]+d_z[i])%hashTableWidth)>hashvalue)
					break;
			}
			return i;
		}
		else if(((d_x[mid]+d_y[mid]+d_z[mid])%hashTableWidth)<hashvalue) && (((mid<offset)&&(d_x[mid+1]+d_y[mid+1]+d_z[mid+1])%hashTableWidth)>hashvalue)||(mid==offset))
		{
			return mid;
		}
		else if ((d_x[mid]+d_y[mid]+d_z[mid])%hashTableWidth)>hashvalue)
			return binarySearch(d_x,d_y,d_z,begin,mid-1,hashvalue);
		else
			return binarySearch(d_x,d_y,d_z,mid+1, offset,hashvalue);
	}



}

__device__ shiftPos(int *d_x, int *d_y, int *d_z,long long *d_indexOrder, long long offset)
{
	int index=blockIdx.x*blockDim.x+threadIdx.x;
	int x,y,z;
	long long indexOrder;
	x=d_x[offset-index];
	y=d_y[offset-index];
	z=d_z[offset-index];
	indexOrder=d_indexOrder[offset-index];
	_syncthreads();
	d_x[offset-index+1]=x;
	d_y[offset-index+1]=y;
	d_z[offset-index+1]=z;
	d_indexOrder[offset-index+1]=indexOrder;
}


__global__ setIndexOrder(long long *d_indexOrder, long long count)
{
	long long index=blockIdx.x*blockDim.x+threadIdx.x;
	for(i=index*(count/numKernels);i<(index+1)*(count/numKernels);i++)
		d_indexOrder[i]=i;
}


__global__ sortAndBFS(int *d_x,int *d_y, int *d_z,long long *d_indexOrder, long long *d_partition_begin, long long *d_partition_last,long long numPartitions, long long count, int *d_labels,long long *d_queue,long long *d_front,long long *d_rear,long long *d_numGroups, long long *d_neighbours)
{
	long long begin, last,i,j,sortPos,x,y,z,numberOfLabels, indexOrder,index, front, rear;
	
	int hashvalue,dx[6],dy[6],dz[6];

	partitionIndex=blockIdx.x*blockDim.x+threadIdx.x;
	if(partitionIndex>=numPartitions)
		return;
	__shared__ long long hashGlobalMemory[hashTableWidth];
	__shared__ int hashedElements[hashTableWidth];
	begin=d_partition_begin[partitionIndex];
	last=d_partition_last[partitionIndex];
	for(i=begin;i<=last;i++)
	{
		hashvalue=(d_x[i]+d_y[i]+d_z[i])%hashTableWidth;
		sortPos=binarySearch<<<1,1>>>(d_x,d_y, d_z,begin,i-1,hashvalue);
		x=d_x[i];
		y=d_y[i];
		z=d_z[i];
		indexzOrder=d_indexOrder[i];
		for(j=i-1;j-numThreadsPerBlock>=sortPos;j=j-numThreadsPerBlock)
		{
			shiftPos<<<1,numThreadsPerBlock>>>(d_x,d_y,d_z,j);
		}
		for(;j>=sortPos;j--)
		{
			shiftPos<<<1,1>>>(d_x,d_y,d_z,j);
		}
		d_x[sortPos]=x;
		d_y[sortPos]=y;
		d_z[sortPos]=z;
		d_indexOrder[sortPos]=indexOrder;

	}
	hashvalue=-1;
	for(i=0;i<hashTableWidth;i++)
	{
		hashGlobalMemory[i]=-1;
		hashedElements[i]=0;
	}

	for(i=begin;i<=last;i++)
	{
		d_labels[i]=-1;
		tempHashValue=(d_x[i]+d_y[i]+d_z[i])%hashTableWidth;
		if(tempHashValue!=hashvalue)
		{
			hashGlobalMemory[tempHashValue]=i;
			hashvalue=i;
		}
		hashedElements[tempHashValue]++;

	}




	if(partitionIndex==0)
	{
		d_front[0]=-1;
		d_rear[0]=-1;
		d_numGroups[0]=0;
	}



	_syncthreads();

	for(i=0;i<count;i++)
	{
		if(partitionIndex==0)
		{
			if(labels[d_indexOrder[i]]==-1)
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
					if(hashGlobalMemory[((dx[j]+dy[j]+dz[j])%hashTableWidth)]!=-1)
					{
						for(k=hashGlobalMemory[((dx[j]+dy[j]+dz[j])%hashTableWidth)];k<(hashGlobalMemory[((dx[j]+dy[j]+dz[j])%hashTableWidth)]+hashedElements[((dx[j]+dy[j]+dz[j])%hashTableWidth)]);k++)
						{
							if(compare(d_x[k],d_y[k],d_z[k],dx[j],dy[j],dz[j])==0)
							{
								if(d_labels[d_indexOrder[k]]==-1)
								{
									d_neighbours[j]=k;
									d_labels[d_indexOrder[k]]=numGroups[0];
								}
								break;
							}

						}
					}
				}
			}
			_syncthreads();
			if(partitionIndex==0)
			{
				for(k=0;k<6;k++)
				{
					if(d_neighbours[k]!=0)
					{
						d_rear[0]++;
						d_queue[d_rear[0]]=d_neighbours[k];
					}
				}
				d_front[0]++;
			}
			_syncthreads();
		}
	}

}


int main()
{
	long long count,i,*d_front,*d_rear,h_front,h_rear,*h_labels,*d_labels,*h_partition_begin,*d_partition_begin, *orderedGroups;
	long long *h_partition_last,*d_partition_last,*d_queue, *h_numGroups, *d_numGroups,*d_indexOrder,*d_neighbours;
	int nx,ny,nz,x,y,z,*h_x,*h_y,*h_z,*d_x,*d_y,*d_z,h_numPartitions,*d_numPartition;
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
	printf("%lld\n",count);
	fclose(fp);
	h_x=(int*)malloc(sizeof(int)*count);
	h_y=(int*)malloc(sizeof(int)*count);
	h_z=(int*)malloc(sizeof(int)*count);
	h_partition_begin=(long long*)malloc(sizeof(long long));
	h_partition_last=(long long*)malloc(sizeof(long long));
	h_labels=(long long*)malloc(sizeof(long long)*count);

	cudaMalloc(&d_x,count*sizeof(int));
	cudaMalloc(&d_y,count*sizeof(int));
	cudaMalloc(&d_z,count*sizeof(int));
	cudaMalloc(&d_partition_begin,numKernels*sizeof(int));
	cudaMalloc(&d_partition_last,numKernels*sizeof(int));
	cudaMalloc(&d_queue,sizeof(long long)*count);
	cudaMalloc(&d_labels,sizeof(int)*count);
	cudaMalloc(&d_front,sizeof(long long));
	cudaMalloc(&d_rear,sizeof(long long));
	cudaMalloc(&d_labels,sizeof(long long)*count);
	cudaMalloc(&d_numPartitions,sizeof(int));
	cudaMalloc(&d_numGroups,sizeof(long long));
	cudaMalloc(&d_indexOrder, sizeof(long long)*count);
	cudaMalloc(&d_neighbours,sizeof(long long)*6); /*0:left, 1: behind, 2:right,3: front, 4: top, 5:bottom*/

	h_front=-1;
	h_rear=-1

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
	cudaMemcpy(d_front,h_front,sizeof(long long),cudaMemcpyHostToDevice);
	cudaMemcpy(d_rear,h_rear,sizeof(long long),cudaMemcpyHostToDevice);
	cudaMemcpy(d_rear,h_rear,sizeof(long long),cudaMemcpyHostToDevice);
	cudaMemcpy(d_partition_begin,h_partition_begin,sizeof(long long),cudaMemcpyHostToDevice);
	cudaMemcpy(d_partition_last,h_partition_last,sizeof(long long),cudaMemcpyHostToDevice);

	setIndexOrder<<<numKernels,1>>>(d_indexOrder,count);

	initializePartitionTable<<<numKernels,1>>>(d_partition_begin,d_partition_last);

	for(i=1;i<=(numKernels/2);i=i*2)
	{

		createPartition<<<i,1>>>(count,d_x,d_y,d_z,d_partition_begin,d_partition_last);

	}

	countPartitions<<<1,1>>>(d_partition_begin,d_partition_last,d_numPartitions);

	cudaMemcpy(h_numPartitions,d_numPartitions,sizeof(int),cudaMemcpyDeviceToHost);

	sortAndBFS<<<numKernels,1>>>(d_x,d_y,d_z,d_labels,d_queue,d_front,d_rear,d_partition_front,d_partition_rear,count);

	cudaMemcpy(h_labels,d_labels,count*sizeof(long long),cudaMemcpyDeviceToHost);
	cudaMemcpy(h_numGroups,d_groups,sizeof(long long),cudaMemcpyDeviceToHost);
	orderedGroups=(long long*)malloc(sizeof(long long)*(h_numGroups+1));
	for(i=1;i<=h_numGroups;i++)
		orderedGroups[i]=-1;
	groupLabelJockey=0;
	for(i=0;i<count;i++)
	if(orderedGroups[h_labels[i]]==-1)
	{
		groupLabelJockey++;
		orderedGroups[h_labels[i]]=groupLabelJockey;
	}
	ofp=fopen("result.txt","w");
	fprintf(ofp,"Nx=%d Ny=%d  Nz=%d Cluster=%lld\nX Y Z id\n",nx,ny,nz,h_numGroups);
	for(i=0;i<count;i++)
		fprintf(ofp,"%d %d %d %lld\n",h_x[i],h_y[i],h_z[i],orderedGroups[h_labels[i]]);
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

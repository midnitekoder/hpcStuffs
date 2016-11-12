#include<stdio.h>
#include<math.h>
#include<limits.h>
#define numKernels 512
#define kernelBatchSize 45
#define numThreads 256
#define hashTableWidth 2048
#define searchKernels 1024
#define numThreadsPerBatch 32

int maxnum(int num1,int num2)
{
if(num1>num2)
return num1;
else
return num2;
}

int maxdimension(int num1,int num2,int num3)
{

return maxnum(num1,maxnum(num2,num3));
}


__device__ void swapTriplets(int *d_x,int *d_y,int *d_z, long long index1,long long index2)
{
int temp1,temp2,temp3;
temp1=d_x[index1];
temp2=d_y[index1];
temp3=d_z[index1];
d_x[index1]=d_x[index2];
d_y[index1]=d_y[index2];
d_z[index1]=d_z[index2];
d_x[index2]=temp1;
d_y[index2]=temp2;
d_z[index2]=temp3;
}


int compare(int x1,int y1, int z1, int x2, int y2, int z2)
{
if((x1==x2) && (y1==y2) && (z1==z2))
return 0;
if((x1>x2) ||((x1==x2) && (y1>y2)) || ((x1==x2) && (y1==y2) && (z1>z2)))
return 1;
return -1;
}

__global__ void createPartition(long long count,int *d_x,int *d_y,int *d_z,int *d_partition_begin,int *d_partition_last)
{
int partitionIndex=blockIdx.x*blockDim.x+threadIdx.x;
long long begin,last,i,j;
int x,y,z;
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
i = (begin - 1);

for ( j = begin; j <= last- 1; j++)
{
if (compare(d_x[j],d_y[j],d_z[j],x,y,z)>=0)
{
    i++;
    swapTriplets (&arr[i], &arr[j]);
}
}
swapTriplets (d_x,d_y,d_z,(i+1),(last));

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


__device__ binarySearch(int *d_x, int *d_y, int *d_z,long long begin, long long offset)
{
long long mid, index=blockIdx.x*blockDim.x+threadIdx.x;
}



__device__ insertionSortParallel(int *d_x,int *d_y, int *d_z, long long *d_partition_begin, long long *d_partition_last, long long count, int partitionIndex)
{
long long begin, last,i,j,sortPos;

begin=d_partition_begin[partitionIndex];
last=d_partition_last[partitionIndex];
for(i=begin;i<=last;i++)
{
sortPos=binarySearch(d_x,d_y, d_z,begin,i-1);
for(j=i-1;j-numThradsPerBlock>=sortPos;j=j-numThreadsPerBlock)
{
shiftPos<<1,numThreadsPerBlock>>(j,d_x,d_y,d_z);
}
for(:j>=sortPos;j--)
{
shiftPos<<1,1>>(j,d_x,d_y,d_z);
}


}


__global__ sortandBfs(int *d_x,int *d_y,int *d_z,int *d_labels,long long *d_queue,long long *d_front,long long *d_rear,long long *d_partition_front,long long *d_partition_rear,long long count)
{
}



int main()
{
long long count,i,*d_front,*d_rear,h_front,h_rear,*h_labels,*d_labels,h_numLabels,*d_numLabels,*h_partition_begin,*d_partition_begin,*h_partition_last,*d_partition_last,*d_queue;
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
cudaMalloc(&d_numLabels,sizeof(long long));
cudaMalloc(&d_numPartitions,sizeof(int));
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

initializePartitionTable<<<numKernels,1>>>(d_partition_begin,d_partition_last);
for(i=1;i<=(numKernels/2);i=i*2)
{
createPartition<<<i,1>>>(count,d_x,d_y,d_z,d_partition_begin,d_partition_last);

}
countPartitions<<<1,1>>>(d_partition_begin,d_partition_last,d_numPartitions);
cudaMemcpy(h_numPartitions,d_numPartitions,sizeof(int),cudaMemcpyDeviceToHost);
sortandBfs<<<numKernels,1>>>(d_x,d_y,d_z,d_labels,d_queue,d_front,d_rear,d_partition_front,d_partition_rear,count);

cudaMemcpy(h_labels,d_labels,count*sizeof(long long),cudaMemcpyDeviceToHost);
cudaMemcpy(h_numLabels,d_numLabels,sizeof(long long),cudaMemcpyDeviceToHost);

ofp=fopen("result.txt","w");
fprintf(ofp,"Nx=%d Ny=%d  Nz=%d Cluster=%lld\nX Y Z id\n",nx,ny,nz,h_numLabels);
for(i=0;i<count;i++)
fprintf(ofp,"%d %d %d %lld\n",h_x[i],h_y[i],h_z[i],h_labels[i]);
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
cudaFree(d_numLabels);
cudaFree(d_partition_begin);
cudaFree(d_partition_last);
cudaFree(d_x);
cudaFree(d_y);
cudaFree(*d_z);
cudaFree(d_queue);

return 0;
}

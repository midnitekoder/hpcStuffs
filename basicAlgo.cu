#include<stdio.h>
#include<math.h>
#include<limits.h>
#define numKernels 512
#define kernelBatchSize 45
#define numThreads 256
#define hashTableWidth 2048
__global__ 
void getQuadrant(int count,int *d_x,int *d_y,int *d_z,int *d_quad,int windowlen,int numDivisions,int nx,int ny,int nz)
{
//printf("yo bitch\n");
int index=blockIdx.x*blockDim.x+threadIdx.x;
int x,y,z;

//printf("yo\n");
if (index<count)
{
x=d_x[index];
y=d_y[index];
z=d_z[index];

printf("count= %d d_x[%d]=%d d_y=%d d_z=%d windowlen=%d numDivisions=%d\n",count,x,index,y,z,windowlen,numDivisions); 
//d_quad[index]=(4*nx*ny*(nz+z)+2*nx*(ny+y)+(nx+x)+(boxvolume-1))/(boxvolume);
d_quad[index]=numDivisions*numDivisions*((nz+z+windowlen-1)/windowlen)+numDivisions*((ny+y+windowlen-1)/windowlen)+((nx+x+windowlen-1)/windowlen);
printf("%d %d\n",index,d_quad[index]);
}



}


__device__ void buildSearchHash(int batchIndex,int count, int *d_x,int *d_y,int *d_z,int *d_quad, int numDivisions,int offsetSize)
{
extern __shared__ int quadCoordinateMap[];
int index, newEntryColumn,i,j;
for(i=0;i<numDivisions;i++)
for(j=1;j<hashTableWidth/4;j++)
quadCoordinateMap[i*(hashTableWidth/4)+j]=INT_MIN;
for(i=0;i<numDivsions;i++)
quadCoordianteMap[i*(hastTableWidth/4)]=0;
index=batchIndex*offsetSize+blockIdx.x*blockDim.x;
for(i=0;i<offsetSize;i++)
{
newEntryColumn=quadCoordinateMap[d_quad[index+i]*(hashTableWidth/4)];
quadCoordinateMap[d_quad[index+i]*(hashTableWidth/4)+newEntryColumn*4+1]=d_x[index+i];
quadCoordinateMap[d_quad[index+i]*(hashTableWidth/4)+newEntryColumn*4+2]=d_y[index+i];
quadCoordinateMap[d_quad[index+i]*(hashTableWidth/4)+newEntryColumn*4+3]=d_z[index+i];
quadCoordinateMap[d_quad[index+i]*(hashTableWidth/4)+newEntryColumn*4+4]=index+i;
quadCoordinateMap[d_quad[index+i]*(hashTableWidth/4)]++;
}

_syncthreads();
__shared__ int findCoordinates[42];
for(i=0;i<count;i++)
{
for(k=0;k<42;k++)
findCoordinates[k]=0;
searchNeighbours<<<1,42>>>(quadCoordinateMap,count,x,y,z,quad,numDivisions,findCoordinates);
for(k=0;k<6;k++)
{
for(l=0;l<7;l++)
if(findCoordinates[k*7+l]!=0 && findCoordinates[k*7+l]!=d_quad[i])
{




/*for(j=0;j<quadCoordinateMap[d_quad[i]*(hashTableWidth/4)];j++)
{
if((quadCoordinateMap[d_quad[i]*(hashTableWidth/4)+j*3+1]==d_x[i]) && (quadCoordinateMap[d_quad[i]*(hashTableWidth/4)+j*3+2]==d_y[i]) && (quadCoordinateMap[d_quad[i]*(hashTableWidth/4)+j*3+3]==d_z[i]))
{
*/





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

int main()
{
long long count,i;
int nx,ny,nz,x,y,z,*h_x,*h_y,*h_z,*h_quad,*d_quad,*d_x,*d_y,*d_z,numDivisions,windowlen,temp;
double cubeRoot;
//long long sharedMemSize=1000;
FILE *fp;
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
h_quad=(int*)malloc(sizeof(int)*count);
cudaMalloc(&d_x,count*sizeof(int));
cudaMalloc(&d_y,count*sizeof(int));
cudaMalloc(&d_z,count*sizeof(int));
cudaMalloc(&d_quad,count*sizeof(int));


fp=fopen("data.txt","r");
fscanf(fp,"Nx=%d Ny=%d Nz=%d",&nx,&ny,&nz);
printf("%d %d %d\n",nx,ny,nz);
for(i=0;i<count;i++)
fscanf(fp,"%d %d %d\n",&h_x[i],&h_y[i],&h_z[i]);
/*printf("yoo");
for(i=0;i<count;i++)
printf("yo %d %d %d\n",h_x[i],h_y[i],h_z[i]);
*/
/*boxx=(2*nx)/numKernels;
boxy=(2*ny)/numKernels;
boxz=(2*nz)/numKernels;
*/
cubeRoot=cbrt(double(numKernels));
numDivisions=cubeRoot;
printf("cubeRoot %lf numKernels %d\n",cubeRoot,numKernels);
windowlen=(2*maxdimension(nx,ny,nz)+cubeRoot-1)/cubeRoot;
printf("numDivisions %d windowlen %d\n",numDivisions,windowlen);
cudaMemcpy(d_x,h_x,count*sizeof(int),cudaMemcpyHostToDevice);
cudaMemcpy(d_y,h_y,count*sizeof(int),cudaMemcpyHostToDevice);
cudaMemcpy(d_z,h_z,count*sizeof(int),cudaMemcpyHostToDevice);
getQuadrant<<<(count+numThreads-1)/numThreads,numThreads>>>(count,d_x,d_y,d_z,d_quad,windowlen,numDivisions,nx,ny,nz);

cudaMemcpy(h_quad,d_quad,count*sizeof(int),cudaMemcpyDeviceToHost);
for(i=0;i<5;i++)
printf("yo4 %d\n",h_quad[i]);




}

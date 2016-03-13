#include <windows.h>
#include <stdio.h>
#include <math.h>
#include <GL/gl.h>
#include <GL/glut.h>
#include <GL/glu.h>
//#include "Vector3D.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
    
void CreateTable(float* arr,int min, int max, int rgb);
void ScanWithCuda();
void CudaUpdate();
__global__ void Scanning_Kernal(unsigned char* vol, unsigned char* val);
__global__ void Scanning_Kernal2(unsigned char* vol, unsigned char* val, float* dev_Dir, float* dev_Up, float* dev_Cross, float* dev_Eye, int Width, int Height, int Depth);


GLubyte map[256*256*4];
float Eye[3] = {0, 0, 0};
bool EyeHasChanged = false;
int Min, Max;

const int Width = 256;
const int Height = 256;
const int Depth = 225;

bool Interpolation=false;

__constant__ int cWidth;
__constant__ int cHeight;
__constant__ int cDepth;
__constant__ bool cInterpolation;

__constant__ float cEye[3];
__constant__ float cU[3];
__constant__ float cCross[3];
__constant__ float cDir[3];
__constant__ float cColor_r[256];
__constant__ float cColor_g[256];
__constant__ float cColor_b[256];
__constant__ float cAlpha[256];


    const int MapSize = Width*Height*Depth*sizeof(unsigned char);
	cudaEvent_t start, stop;//타임이벤트
	float Time;
    unsigned char *vol = new unsigned char[Width * Height * Depth];
    unsigned char *val = new unsigned char[Width * Height];
   
    float Look[3] = {128,128,112.5};
	float Up[3] = {0, 1, 0};
	float Cross[3] = {0, 0, 0};
	float Dir[3] = {0, 0, 0};
	float u[3] = {0, 0, 0};
    
    float Color_r[256];
    float Color_g[256];
    float Color_b[256];
    float Alpha[256];
	float *dev_Eye;
	float *dev_Dir;
	float *dev_U;
	float *dev_Cross;
	unsigned char* dev_vol;
	unsigned char* dev_val;


void myreshape(int w, int h)
{
 
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);
	glEnable(GL_TEXTURE_2D);
 
	GLuint texId;
	glGenTextures(1, &texId);
	glBindTexture(GL_TEXTURE_2D, texId);
 
	glTexImage2D(
		GL_TEXTURE_2D, 0, GL_RGBA,
		Width, Height, 0,
		GL_RGBA, GL_UNSIGNED_BYTE,
		map
		);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
 
 
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glOrtho(0, Width, 0, Height, 0, 1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	glBindTexture(GL_TEXTURE_2D, texId);
}

void mykeyboard(unsigned char keyPressed, int x, int y)
{ 
 switch (keyPressed)//키를 누루면 위치 이동
 {
 case 'a'://left
  Eye[0]-=5;
  EyeHasChanged = true;
  break;
 case 'd'://right
  Eye[0]+=5;
  EyeHasChanged = true;
  break;
 case 'w'://up
  Eye[1]+=5;
  EyeHasChanged = true;
  break;
 case 's'://down
  Eye[1]-=5;
  EyeHasChanged = true;
  break;
 case 'q'://down
  Eye[2]-=5;
  EyeHasChanged = true;
  break;
 case 'e'://down
  Eye[2]+=5;
  EyeHasChanged = true;
  break;
 case 'z'://down
  Min-=5;
  Max-=5;
  break;
 case 'x'://down
  Min+=5;
  Max+=5;
  break;
 case 'A'://left
  Eye[0]-=10;
  EyeHasChanged = true;
  break;
 case 'D'://right
  Eye[0]+=10;
  EyeHasChanged = true;
  break;
 case 'W'://up
  Eye[1]+=10;
  EyeHasChanged = true;
  break;
 case 'S'://down
  Eye[1]-=10;
  EyeHasChanged = true;
  break;
 case 'Q'://down
  Eye[2]-=10;
  EyeHasChanged = true;
  break;
 case 'E'://down
  Eye[2]+=10;
  EyeHasChanged = true;
  break;
 case 'Z'://down
  Min-=10;
  Max-=10;
  break;
 case 'X'://down
  Min+=10;
  Max+=10;
  break;
case '['://ZoomIn
  Cross[0]*=0.8;
  Cross[1]*=0.8;
  Cross[2]*=0.8;
  u[0]*=0.8;
  u[1]*=0.8;
  u[2]*=0.8;
  EyeHasChanged = false;
  break;
 case ']'://ZoomOut
  Cross[0]*=1.2;
  Cross[1]*=1.2;
  Cross[2]*=1.2;
  u[0]*=1.2;
  u[1]*=1.2;
  u[2]*=1.2;
  EyeHasChanged = false;
  break;
 }
 CudaUpdate();
 glutPostRedisplay();
}

void MyMouseWheelFunc(int wheel, int direction, int x, int y)
{
    if(direction > 0)
    {
        Min-=10;
        Max-=10;
    }
    else   
    {
        Min+=10;
        Max+=10;
    }

 CudaUpdate();
 glutPostRedisplay();

}

void MyDisplay() {

    
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);
	glEnable(GL_TEXTURE_2D);
 
	GLuint texId;
	glGenTextures(1, &texId);
	glBindTexture(GL_TEXTURE_2D, texId);
 
	glTexImage2D(
		GL_TEXTURE_2D, 0, GL_RGBA,
		Width, Height, 0,
		GL_RGBA, GL_UNSIGNED_BYTE,
		map
		);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
 
 
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glOrtho(0, Width, 0, Height, 0, 1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	glBindTexture(GL_TEXTURE_2D, texId);
 
	glBegin(GL_QUADS);
	glTexCoord2f(0.0, 0.0); glVertex2f(0, Height);
	glTexCoord2f(1.0, 0.0); glVertex2f(Width, Height);
	glTexCoord2f(1.0, 1.0); glVertex2f(Width, 0);
	glTexCoord2f(0.0, 1.0); glVertex2f(0, 0);
	glEnd();
	
	glFlush();
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(600, 600);
	glutInitWindowPosition(200, 200);
	glutCreateWindow("openGL Sample Drawing");
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
    ScanWithCuda();
	glutDisplayFunc(MyDisplay);
    //glutReshapeFunc(myreshape);
    glutKeyboardFunc(mykeyboard);
    //glutMouseWheelFunc(MyMouseWheelFunc);
	glutMainLoop();
	return 0;
}

void ScanWithCuda(){
    
    
    int Choice;
    
    cudaError_t cudaStatus;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
    
    
	//빅헤드 불러오기
	FILE *fp = fopen("C:/Bighead.den", "rb");
	fread(vol, Width*Height*Depth, 1, fp);
	fclose(fp);
	

    printf("Input eye position and min, max\nx y z min max\n");
	scanf("%f %f %f %d %d", &Eye[0], &Eye[1], &Eye[2], &Min, &Max);
	
    CreateTable(Color_r,Min,Max,1);
    CreateTable(Color_g,Min,Max,2);
    CreateTable(Color_b,Min,Max,3);
    CreateTable(Alpha,Min,Max,0);
    
    printf("What you want to try?\n1. __constant__\n2. __global__\n");
	scanf("%d", &Choice);
	
    int inter;
    printf("Do you want to Interpolate?\n1. yes \n2. no\n");
    scanf("%d",&inter);
    if(inter == 1) Interpolation = true;
	Dir[0] = Look[0] - Eye[0];
	Dir[1] = Look[1] - Eye[1];
	Dir[2] = Look[2] - Eye[2];
    //Normalization
	float scalar = sqrt( (Dir[0] * Dir[0]) + (Dir[1] * Dir[1]) + (Dir[2] * Dir[2]) );
	for(int i=0; i<3; i++){
		Dir[i] /= scalar;
	}
	//Cross
	Cross[0] = Up[1]*Dir[2] - Up[2]*Dir[1];
	Cross[1] = Up[2]*Dir[0] - Up[0]*Dir[2];
	Cross[2] = Up[0]*Dir[1] - Up[1]*Dir[0];
    //Normalization
    scalar = sqrt( (Cross[0] * Cross[0]) + (Cross[1] * Cross[1]) + (Cross[2] * Cross[2]) );
	for(int i=0; i<3; i++){
		Cross[i] /= scalar;
	}
	//U
	u[0] = Dir[1]*Cross[2] - Dir[2]*Cross[1];
	u[1] = Dir[2]*Cross[0] - Dir[0]*Cross[2];
	u[2] = Dir[0]*Cross[1] - Dir[1]*Cross[0];

	cudaMalloc((void**)&dev_vol, MapSize);
	cudaMalloc((void**)&dev_val, Width*Height*4*sizeof(unsigned char));
	cudaMalloc((void**)&dev_Eye, 3*sizeof(float));
	cudaMalloc((void**)&dev_Dir, 3*sizeof(float));
	cudaMalloc((void**)&dev_U, 3*sizeof(float));
	cudaMalloc((void**)&dev_Cross, 3*sizeof(float));


	cudaMemcpy(dev_vol, vol, MapSize, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_Eye, Eye, 3*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_Dir, Dir, 3*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_U, u, 3*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_Cross, Cross, 3*sizeof(float), cudaMemcpyHostToDevice);

    //상수 데이터 전송
    cudaStatus = cudaMemcpyToSymbol(cWidth,&Width, sizeof(int));
    cudaStatus = cudaMemcpyToSymbol(cHeight,&Height, sizeof(int));
    cudaStatus = cudaMemcpyToSymbol(cDepth,&Depth, sizeof(int));

    cudaStatus = cudaMemcpyToSymbol(cInterpolation,&Interpolation, sizeof(bool));

    cudaStatus = cudaMemcpyToSymbol(cEye,Eye, 3*sizeof(float));
    cudaStatus = cudaMemcpyToSymbol(cU,u, 3*sizeof(float));
    cudaStatus = cudaMemcpyToSymbol(cCross,Cross, 3*sizeof(float));
    cudaStatus = cudaMemcpyToSymbol(cDir,Dir, 3*sizeof(float));

    cudaStatus = cudaMemcpyToSymbol(cColor_r,Color_r, 256*sizeof(float));
    cudaStatus = cudaMemcpyToSymbol(cColor_g,Color_g, 256*sizeof(float));
    cudaStatus = cudaMemcpyToSymbol(cColor_b,Color_b, 256*sizeof(float));
    cudaStatus = cudaMemcpyToSymbol(cAlpha,Alpha, 256*sizeof(float));
    if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpyToSymbol failed!\n");
            //goto Error;
        }

	cudaEventRecord(start,0);//시간측정 스타트
    (Choice == 1) ? Scanning_Kernal<<<256,256>>>(dev_vol, dev_val) : Scanning_Kernal2<<<256,256>>>(dev_vol, dev_val, dev_Dir, dev_U, dev_Cross, dev_Eye, Width, Height, Depth);
	cudaEventRecord(stop,0);//시간측정 종료
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&Time,start,stop);

	cudaMemcpy(map, dev_val, Width*Height*4*sizeof(unsigned char), cudaMemcpyDeviceToHost);
    
	

	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	printf("%d번 커널함수 소요시간 : %f msec\n",Choice, Time);

}


void CreateTable(float* arr,int min, int max, int rgb){
    int i=0, bone_min=50, bone_max=100;
    float j=1/(float)(max-min);
    if(rgb == 0){
        for( ;i<min; i++)
             arr[i] = 0;
        for( ;i<max; i++)
             arr[i] = arr[i-1]+j;
         for( ;i<256; i++)
             arr[i] = 1;
    }
    if(rgb == 1){
        for( ;i<bone_min; i++)
             arr[i] = 0;
        for( ;i<bone_max; i++)
             arr[i] = 0.9;
         for( ;i<256; i++)
             arr[i] = 1.0;
    }
    if(rgb == 2){
        for( ;i<bone_min; i++)
             arr[i] = 0;
        for( ;i<bone_max; i++)
             arr[i] = 0.3;
         for( ;i<256; i++)
             arr[i] = 1.0;
    }
    if(rgb == 3){
        for( ;i<bone_min; i++)
             arr[i] = 0;
        for( ;i<bone_max; i++)
             arr[i] = 0.1;
         for( ;i<256; i++)
             arr[i] = 0.2;
    }
}

void CudaUpdate()
{

    cudaError_t cudaStatus;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
    
    CreateTable(Color_r,Min,Max,1);
    CreateTable(Color_g,Min,Max,2);
    CreateTable(Color_b,Min,Max,3);
    CreateTable(Alpha,Min,Max,0);

    //눈위치 변경에 따른 업데이트
    if(EyeHasChanged){
        Dir[0] = Look[0] - Eye[0];
	    Dir[1] = Look[1] - Eye[1];
	    Dir[2] = Look[2] - Eye[2];

        //Normalization
	    float scalar = sqrt( (Dir[0] * Dir[0]) + (Dir[1] * Dir[1]) + (Dir[2] * Dir[2]) );
	    for(int i=0; i<3; i++){
		    Dir[i] /= scalar;
	    }
	    //Cross
	    Cross[0] = Up[1]*Dir[2] - Up[2]*Dir[1];
	    Cross[1] = Up[2]*Dir[0] - Up[0]*Dir[2];
	    Cross[2] = Up[0]*Dir[1] - Up[1]*Dir[0];
        //Normalization
        scalar = sqrt( (Cross[0] * Cross[0]) + (Cross[1] * Cross[1]) + (Cross[2] * Cross[2]) );
	    for(int i=0; i<3; i++){
		    Cross[i] /= scalar;
	    }
	    //U
	    u[0] = Dir[1]*Cross[2] - Dir[2]*Cross[1];
	    u[1] = Dir[2]*Cross[0] - Dir[0]*Cross[2];
	    u[2] = Dir[0]*Cross[1] - Dir[1]*Cross[0];

        EyeHasChanged == false;
    }

	cudaMemcpy(dev_Eye, Eye, 3*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_Dir, Dir, 3*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_U, u, 3*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_Cross, Cross, 3*sizeof(float), cudaMemcpyHostToDevice);

    //상수 데이터 전송
    cudaStatus = cudaMemcpyToSymbol(cEye,Eye, 3*sizeof(float));
    cudaStatus = cudaMemcpyToSymbol(cU,u, 3*sizeof(float));
    cudaStatus = cudaMemcpyToSymbol(cCross,Cross, 3*sizeof(float));
    cudaStatus = cudaMemcpyToSymbol(cDir,Dir, 3*sizeof(float));

    cudaStatus = cudaMemcpyToSymbol(cColor_r,Color_r, 256*sizeof(float));
    cudaStatus = cudaMemcpyToSymbol(cColor_g,Color_g, 256*sizeof(float));
    cudaStatus = cudaMemcpyToSymbol(cColor_b,Color_b, 256*sizeof(float));
    cudaStatus = cudaMemcpyToSymbol(cAlpha,Alpha, 256*sizeof(float));
    if (cudaStatus != cudaSuccess) {
            fprintf(stderr, "cudaMemcpyToSymbol failed!\n");
            //goto Error;
    }

	cudaEventRecord(start,0);//시간측정 스타트
    Scanning_Kernal<<<256,256>>>(dev_vol, dev_val);
	cudaEventRecord(stop,0);//시간측정 종료
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&Time,start,stop);

	cudaMemcpy(map, dev_val, Width*Height*4*sizeof(unsigned char), cudaMemcpyDeviceToHost);
    
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	printf("업데이트 커널함수 소요시간 : %f msec\n", Time);

}

__device__ unsigned char Interpolation_Kernal(float* posf, unsigned char* vol) {
    
    if(posf[2] > 220) return 0;
    float back_a = vol[(cWidth*cHeight)*(int)posf[2] + cWidth*((int)posf[1]+1) + (int)posf[0]];
    float back_b = vol[(cWidth*cHeight)*(int)posf[2] + cWidth*((int)posf[1]+1) + (int)posf[0]+1];
    float back = back_a * (1.0 - (posf[0] - (int)posf[0])) + back_b * (posf[0] - (int)posf[0]);  //위 뒤쪽

    float front_a = vol[(cWidth*cHeight)*((int)posf[2]+1) + cWidth*((int)posf[1]+1) + (int)posf[0]];
    float front_b = vol[(cWidth*cHeight)*((int)posf[2]+1) + cWidth*((int)posf[1]+1) + (int)posf[0]+1];
    float front = front_a * (1.0 - (posf[0] - (int)posf[0])) + front_b * (posf[0] - (int)posf[0]);  // 위 앞쪽

    float UpSide = back * (1.0 - (posf[1] - (int)posf[1])) + front * (posf[1] - (int)posf[1]);  //위
    

    back_a = vol[(cWidth*cHeight)*(int)posf[2] + cWidth*((int)posf[1]) + (int)posf[0]];
    back_b = vol[(cWidth*cHeight)*(int)posf[2] + cWidth*((int)posf[1]) + (int)posf[0]+1];
    back = back_a * (1.0 - (posf[0] - (int)posf[0])) + back_b * (posf[0] - (int)posf[0]); //아래 뒤쪽

    front_a = vol[(cWidth*cHeight)*((int)posf[2]+1) + cWidth*((int)posf[1]) + (int)posf[0]];
    front_b = vol[(cWidth*cHeight)*((int)posf[2]+1) + cWidth*((int)posf[1]) + (int)posf[0]+1];
    front = front_a * (1.0 - (posf[0] - (int)posf[0])) + front_b * (posf[0] - (int)posf[0]); //아래 앞쪽

    float BottomSide = back * (1.0 - (posf[1] - (int)posf[1])) + front * (posf[1] - (int)posf[1]); //아래
    
    unsigned char output = UpSide * (1.0 - (posf[2] - (int)posf[2])) + BottomSide * (posf[2] - (int)posf[2]);//가운데
    return output;
}


__global__ void Scanning_Kernal(unsigned char* vol, unsigned char* val) {

    float a_new=0,a_old=0, r_new=0,r_old=0, g_new=0,g_old=0, b_new=0,b_old=0,posf[3];
	int pos[3],tpos[3];
    float start[3];
	tpos[0] = threadIdx.x-128;
	tpos[1] = blockIdx.x-128;
	tpos[2] = cWidth*(tpos[1]+128)+tpos[0]+128;
	unsigned char d,a_sum=0,c_sum=0;
    start[0] = cEye[0] + cCross[0]*tpos[0] + cU[0]*tpos[1];
    start[1] = cEye[1] + cCross[1]*tpos[0] + cU[1]*tpos[1];
    start[2] = cEye[2] + cCross[2]*tpos[0] + cU[2]*tpos[1];

	for(int k=0;k<500; k++){
		    posf[0] = start[0] + cDir[0]*k;
		    posf[1] = start[1] + cDir[1]*k;
		    posf[2] = start[2] + cDir[2]*k;
		if(posf[0]>=0 && posf[0]<256 && posf[1]>=0 && posf[1]<256 && posf[2]>=0 && posf[2]<225){
            if(cInterpolation) d = Interpolation_Kernal(posf, vol);
			else d = vol[(cWidth*cHeight)*(int)posf[2] + cWidth*(int)posf[1] + (int)posf[0]];
            a_new = a_old + (1-a_old) * cAlpha[d];
            r_new = r_old + (1-a_old) * cColor_r[d] * cAlpha[d];
            g_new = g_old + (1-a_old) * cColor_g[d] * cAlpha[d];
            b_new = b_old + (1-a_old) * cColor_b[d] * cAlpha[d];
            r_old = r_new;
            g_old = g_new;
            b_old = b_new;
            a_old = a_new;
        }
	}
	val[tpos[2]*4+0] = r_new*255;
	val[tpos[2]*4+1] = g_new*255;
	val[tpos[2]*4+2] = b_new*255;
	val[tpos[2]*4+3] = 0xff;
}


__global__ void Scanning_Kernal2(unsigned char* vol, unsigned char* val, float* dev_Dir, float* dev_U, float* dev_Cross, float* dev_Eye, int Width, int Height, int Depth) {

	int pos[3],tpos[3];
	tpos[0] = threadIdx.x-128;
	tpos[1] = blockIdx.x-128;
	tpos[2] = Width*(tpos[1]+128)+tpos[0]+128;
	unsigned char Found_max = 0;

	for(int k=0;k<500; k++){
		pos[0] = dev_Eye[0] + dev_Cross[0]*tpos[0] + dev_U[0]*tpos[1] + dev_Dir[0]*k;
		pos[1] = dev_Eye[1] + dev_Cross[1]*tpos[0] + dev_U[1]*tpos[1] + dev_Dir[1]*k;
		pos[2] = dev_Eye[2] + dev_Cross[2]*tpos[0] + dev_U[2]*tpos[1] + dev_Dir[2]*k;

		if(pos[0]>=0 && pos[0]<256 && pos[1]>=0 && pos[1]<256 && pos[2]>=0 && pos[2]<225)
			(Found_max > vol[(Width*Height)*pos[2] + Width*pos[1] + pos[0]]) ? 1 : Found_max = vol[(Width*Height)*pos[2] + Width*pos[1] + pos[0]];
	}
	val[tpos[2]] = Found_max;
}

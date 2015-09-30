//menyelesaikan persamaan eliptik 2D untuk persegi (X = Y = N)
#include <stdio.h>
#include <time.h>
#include <math.h>

void writeToFile(int X, int Y, int length, char* filename1, float* arr) {
	//print ke file
	FILE * pFile = fopen(filename1,"w");
	fprintf(pFile, "%d\n%d\n",X,Y);
	for (int i=0;i<Y;i++) {
		for(int j=0;j<X;j++) {
			fprintf(pFile,"%.2f ",arr[(i*X)+j]);
		}
		fprintf(pFile, "\n");
	}
	fclose(pFile);
}

float absol(float input) {
	if(input < 0) {
		return (-input);
	}
	else {
		return input;
	}
}

void cpuJacobi (float * input, float * output, float * error, char ad_n, char ad_e, char ad_s, char ad_w, int N, int jumlahElemen) {
	for (int i = 0; i < jumlahElemen; i++) {
		//float A = input[i], B = input[i], C = input[i], D = input[i];
		//tetap pada batas:
		//y = 0 --> i % N = 0
		//x = 0 --> i / N = 0
		//y = y --> i % N = N-1
		//x = x --> i / N = N-1
		if(i == 0 || i == N-1 || ((i / N == 0) && (i % N == N-1)) || ((i % N == 0) && (i / N == N-1))) {
			output[i] = input[i];
		}
		else if((i%N == 0) || (i/N == 0) || (i%N == N-1) || (i/N == N-1)) {
			if(i/N == 0 && ad_n) {
				output[i] = (input[i+1] + input[i+N] + input[i-1])/3;
			}
			else if(i%N == N-1 && ad_e) {
				output[i] = (input[i-1] + input[i+N] + input[i-N])/3;
			}
			else if(i/N == N-1 && ad_s) {
				output[i] = (input[i+1] + input[i-1] + input[i-N])/3;
			}
			else if(i%N == 0 && ad_w){
				output[i] = (input[i+1] + input[i+N] + input[i-N])/3;
			}
			else{
				output[i] = input[i];
			}
		}
		else {	//tidak pada batas
			//lakukan perhitungan Jacobi
			output[i] = (input[i-1] + input[i+1] + input[i-N] + input[i+N])/(4);
		}
		//hitung error
		error[i] = output[i]-input[i];

		//input[i] = output[i];
	}
	for(int i = 0; i< jumlahElemen; i++) {
		input[i] = output[i];
	}
}

void cpuPGS(float * input, float * output, float * error, char ad_n, char ad_e, char ad_s, char ad_w, int N, int jumlahElemen) {
	//float A = 0, B = 0, C = 0, D = 0;
	for (int i = 0; i < jumlahElemen; i++) {
		//tetap pada batas:
		//y = 0 --> i % N = 0
		//x = 0 --> i / N = 0
		//y = y --> i % N = N-1
		//x = x --> i / N = N-1
		if(i == 0 || i == N-1 || ((i / N == 0) && (i % N == N-1)) || ((i % N == 0) && (i / N == N-1))) {
			output[i] = input[i];
		}
		else if((i%N == 0) || (i/N == 0) || (i%N == N-1) || (i/N == N-1)) {
			if(i/N == 0 && ad_n) {
				output[i] = (input[i+1] + input[i+N] + input[i-1])/3;
			}
			else if(i%N == N-1 && ad_e) {
				output[i] = (input[i-1] + input[i+N] + input[i-N])/3;
			}
			else if(i/N == N-1 && ad_s) {
				output[i] = (input[i+1] + input[i-1] + input[i-N])/3;
			}
			else if(i%N == 0 && ad_w){
				output[i] = (input[i+1] + input[i+N] + input[i-N])/3;
			}
			else{
				output[i] = input[i];
			}
		}
		else {	//tidak pada batas
			//lakukan perhitungan PGS
			output[i] = (output[i-1] + input[i+1] + output[i-N] + input[i+N])/(4);
		}
		//hitung error
		error[i] = output[i]-input[i];

		//input[i] = output[i];
	}
	for(int i = 0; i< jumlahElemen; i++) {
		input[i] = output[i];
	}
}

void cpuPSOR(float * input, float * output, float * error, char ad_n, char ad_e, char ad_s, char ad_w, int N, int jumlahElemen, float omega) {
	//float A = 0, B = 0, C = 0, D = 0;
	for (int i = 0; i < jumlahElemen; i++) {
		//tetap pada batas:
		//y = 0 --> i % N = 0
		//x = 0 --> i / N = 0
		//y = y --> i % N = N-1
		//x = x --> i / N = N-1
		if(i == 0 || i == N-1 || ((i / N == 0) && (i % N == N-1)) || ((i % N == 0) && (i / N == N-1))) {
			output[i] = input[i];
		}
		else if((i%N == 0) || (i/N == 0) || (i%N == N-1) || (i/N == N-1)) {
			if(i/N == 0 && ad_n) {
				output[i] = (input[i+1] + input[i+N] + input[i-1])/3;
			}
			else if(i%N == N-1 && ad_e) {
				output[i] = (input[i-1] + input[i+N] + input[i-N])/3;
			}
			else if(i/N == N-1 && ad_s) {
				output[i] = (input[i+1] + input[i-1] + input[i-N])/3;
			}
			else if(i%N == 0 && ad_w){
				output[i] = (input[i+1] + input[i+N] + input[i-N])/3;
			}
			else{
				output[i] = input[i];
			}
		}
		else {	//tidak pada batas
			//lakukan perhitungan PSOR
			output[i] = input[i] + ((output[i-1] + input[i+1] + output[i-N] + input[i+N] - (input[i]*4))*omega/(4));
		}
		//hitung error
		error[i] = output[i]-input[i];

		//input[i] = output[i];
	}
	for(int i = 0; i< jumlahElemen; i++) {
		input[i] = output[i];
	}
}

__global__ void gpuJacobi (float * input, float * error, char ad_n, char ad_e, char ad_s, char ad_w, int N) {
	//float A = 0, B = 0, C = 0, D = 0;
	//thread ID
	int tidx = (blockIdx.x * blockDim.x) + threadIdx.x;
	int tidy = (blockIdx.y * blockDim.y) + threadIdx.y;
	int tid = (tidy * N) + tidx;	// This gives every thread a unique ID.
	
	if(tid<(N*N)) {
		float hasil = 0;
		//tetap pada batas:
		//y = 0 --> i % N = 0
		//x = 0 --> i / N = 0
		//y = y --> i % N = N-1
		//x = x --> i / N = N-1
		if(tid == 0 || tid == N-1 || ((tid / N == 0) && (tid % N == N-1)) || ((tid % N == 0) && (tid / N == N-1))) {
			hasil = input[tid];
		}
		else if((tid%N == 0) || (tid/N == 0) || (tid%N == N-1) || (tid/N == N-1)) {
			if(tid/N == 0 && ad_n) {
				hasil = (input[tid+1] + input[tid+N] + input[tid-1])/3;
			}
			else if(tid%N == N-1 && ad_e) {
				hasil = (input[tid-1] + input[tid+N] + input[tid-N])/3;
			}
			else if(tid/N == N-1 && ad_s) {
				hasil = (input[tid+1] + input[tid-1] + input[tid-N])/3;
			}
			else if(tid%N == 0 && ad_w){
				hasil = (input[tid+1] + input[tid+N] + input[tid-N])/3;
			}
			else{
				hasil = input[tid];
			}
		}
		else {	//tidak pada batas
			//lakukan perhitungan Jacobi
			hasil = (input[tid-1] + input[tid+1] + input[tid-N] + input[tid+N])/(4);
			error[tid] = hasil-input[tid];
		}
		input[tid] = hasil;
		//hitung error
	}
}

__global__ void gpuPGSRed (float * input, float * error, char ad_n, char ad_e, char ad_s, char ad_w, int N) {
	//float A = 0, B = 0, C = 0, D = 0
	float hasil = 0;
	//thread ID
	int tidx = (blockIdx.x * blockDim.x) + threadIdx.x;
	int tidy = (blockIdx.y * blockDim.y) + threadIdx.y;
	int tid = (tidy * N) + tidx;	// This gives every thread a unique ID.
	//int blid = (blockIdx.y * blockDim.x) + blockIdx.x;
	
	if(tid<(N*N)) {
		//tetap pada batas:
		//y = 0 --> i % N = 0
		//x = 0 --> i / N = 0
		//y = y --> i % N = N-1
		//x = x --> i / N = N-1
		if (tid%2 != (tid/N)%2) {
			if(tid == 0 || tid == N-1 || ((tid / N == 0) && (tid % N == N-1)) || ((tid % N == 0) && (tid / N == N-1))) {
				hasil = input[tid];
			}
			else if((tid%N == 0) || (tid/N == 0) || (tid%N == N-1) || (tid/N == N-1)) {
				if(tid/N == 0 && ad_n) {
					hasil = (input[tid+1] + input[tid+N] + input[tid-1])/3;
				}
				else if(tid%N == N-1 && ad_e) {
					hasil = (input[tid-1] + input[tid+N] + input[tid-N])/3;
				}
				else if(tid/N == N-1 && ad_s) {
					hasil = (input[tid+1] + input[tid-1] + input[tid-N])/3;
				}
				else if(tid%N == 0 && ad_w){
					hasil = (input[tid+1] + input[tid+N] + input[tid-N])/3;
				}
				else{
					hasil = input[tid];
				}
			}
			else {	//tidak pada batas
				//lakukan perhitungan PGS
				//if (tid%2 != (tid/N)%2) {
					hasil = (input[tid-1] + input[tid+1] + input[tid-N] + input[tid+N])/(4);
					error[tid] = hasil-input[tid];
				//}
			}
		input[tid] = hasil;
		}
		
		//hitung error
	}
}

__global__ void gpuPGSBlack (float * input, float * error, char ad_n, char ad_e, char ad_s, char ad_w, int N) {
	//float A = 0, B = 0, C = 0, D = 0;
	float hasil = 0;
	//thread ID
	int tidx = (blockIdx.x * blockDim.x) + threadIdx.x;
	int tidy = (blockIdx.y * blockDim.y) + threadIdx.y;
	int tid = (tidy * N) + tidx;	// This gives every thread a unique ID.
	//int blid = (blockIdx.y * blockDim.x) + blockIdx.x;
	
	if(tid<(N*N)) {
		//tetap pada batas:
		//y = 0 --> i % N = 0
		//x = 0 --> i / N = 0
		//y = y --> i % N = N-1
		//x = x --> i / N = N-1
		if (tid%2 == (tid/N)%2) {
			if(tid == 0 || tid == N-1 || ((tid / N == 0) && (tid % N == N-1)) || ((tid % N == 0) && (tid / N == N-1))) {
				hasil = input[tid];
			}
			else if((tid%N == 0) || (tid/N == 0) || (tid%N == N-1) || (tid/N == N-1)) {
				if(tid/N == 0 && ad_n) {
					hasil = (input[tid+1] + input[tid+N] + input[tid-1])/3;
				}
				else if(tid%N == N-1 && ad_e) {
					hasil = (input[tid-1] + input[tid+N] + input[tid-N])/3;
				}
				else if(tid/N == N-1 && ad_s) {
					hasil = (input[tid+1] + input[tid-1] + input[tid-N])/3;
				}
				else if(tid%N == 0 && ad_w){
					hasil = (input[tid+1] + input[tid+N] + input[tid-N])/3;
				}
				else{
					hasil = input[tid];
				}
			}
			else {	//tidak pada batas
				//lakukan perhitungan PGS
				//if (tid%2 != (tid/N)%2) {
					hasil = (input[tid-1] + input[tid+1] + input[tid-N] + input[tid+N])/(4);
					error[tid] = hasil-input[tid];
				//}
			}
		input[tid] = hasil;
		}
	}
}

__global__ void gpuPSORRed (float * input, float * error, char ad_n, char ad_e, char ad_s, char ad_w, int N, float omega) {
	//float A = 0, B = 0, C = 0, D = 0;
	float hasil = 0;
	//thread ID
	int tidx = (blockIdx.x * blockDim.x) + threadIdx.x;
	int tidy = (blockIdx.y * blockDim.y) + threadIdx.y;
	int tid = (tidy * N) + tidx;	// This gives every thread a unique ID.
	//int blid = (blockIdx.y * blockDim.x) + blockIdx.x;
	
	if(tid<(N*N)) {
		//tetap pada batas:
		//y = 0 --> i % N = 0
		//x = 0 --> i / N = 0
		//y = y --> i % N = N-1
		//x = x --> i / N = N-1
		if (tid%2 != (tid/N)%2) {
			if(tid == 0 || tid == N-1 || ((tid / N == 0) && (tid % N == N-1)) || ((tid % N == 0) && (tid / N == N-1))) {
				hasil = input[tid];
			}
			else if((tid%N == 0) || (tid/N == 0) || (tid%N == N-1) || (tid/N == N-1)) {
				if(tid/N == 0 && ad_n) {
					hasil = (input[tid+1] + input[tid+N] + input[tid-1])/3;
				}
				else if(tid%N == N-1 && ad_e) {
					hasil = (input[tid-1] + input[tid+N] + input[tid-N])/3;
				}
				else if(tid/N == N-1 && ad_s) {
					hasil = (input[tid+1] + input[tid-1] + input[tid-N])/3;
				}
				else if(tid%N == 0 && ad_w){
					hasil = (input[tid+1] + input[tid+N] + input[tid-N])/3;
				}
				else{
					hasil = input[tid];
				}
			}
			else {	//tidak pada batas
				//lakukan perhitungan PSOR
				//if (tid%2 != (tid/N)%2) {
					hasil = input[tid] + ((input[tid-1] + input[tid+1] + input[tid-N] + input[tid+N] - (input[tid]*4))*(omega/4));
					error[tid] = hasil-input[tid];		
				//}
			}
			//hitung error
			input[tid] = hasil;
		}
	}
}

__global__ void gpuPSORBlack (float * input, float * error, char ad_n, char ad_e, char ad_s, char ad_w, int N, float omega) {
	//float A = 0, B = 0, C = 0, D = 0;
	float hasil = 0;
	//thread ID
	int tidx = (blockIdx.x * blockDim.x) + threadIdx.x;
	int tidy = (blockIdx.y * blockDim.y) + threadIdx.y;
	int tid = (tidy * N) + tidx;	// This gives every thread a unique ID.
	//int blid = (blockIdx.y * blockDim.x) + blockIdx.x;
	
	if(tid<(N*N)) {
		//tetap pada batas:
		//y = 0 --> i % N = 0
		//x = 0 --> i / N = 0
		//y = y --> i % N = N-1
		//x = x --> i / N = N-1
		if (tid%2 == (tid/N)%2) {
			if(tid == 0 || tid == N-1 || ((tid / N == 0) && (tid % N == N-1)) || ((tid % N == 0) && (tid / N == N-1))) {
				hasil = input[tid];
			}
			else if((tid%N == 0) || (tid/N == 0) || (tid%N == N-1) || (tid/N == N-1)) {
				if(tid/N == 0 && ad_n) {
					hasil = (input[tid+1] + input[tid+N] + input[tid-1])/3;
				}
				else if(tid%N == N-1 && ad_e) {
					hasil = (input[tid-1] + input[tid+N] + input[tid-N])/3;
				}
				else if(tid/N == N-1 && ad_s) {
					hasil = (input[tid+1] + input[tid-1] + input[tid-N])/3;
				}
				else if(tid%N == 0 && ad_w){
					hasil = (input[tid+1] + input[tid+N] + input[tid-N])/3;
				}
				else{
					hasil = input[tid];
				}
			}
			else {	//tidak pada batas
				//lakukan perhitungan PSOR
					hasil = input[tid] + ((input[tid-1] + input[tid+1] + input[tid-N] + input[tid+N] - (input[tid]*4))*(omega/4));
					error[tid] = hasil-input[tid];		
			}
			//hitung error
			input[tid] = hasil;
		}
	}
}

int main (int argc, char** argv) {
	int SISI_MATRIKS = 384;
	int totalElemen = SISI_MATRIKS*SISI_MATRIKS;
	//float pi =  3.14159265;

	char ad_n = 0, ad_e = 0, ad_s = 0, ad_w = 1;
	float omega = 1.6;
	float *awal_in, *awal_out, *awal_err;
	float *host_in, *host_out, *host_err;
	float *dev_in, *dev_err;
	float total_error = -1;
	float MAX_ERROR = 0.05;
	int CPUIter = 0, CPUIter2 = 0, CPUIter3 = 0, GPUIter = 0, GPUIter2 = 0, GPUIter3 = 0;
	float CPUTime = 0, CPUTime2 = 0, CPUTime3 = 0, GPUTime = 0, GPUTime2 = 0, GPUTime3 = 0;
	dim3 jumlahBlock, threadPerBlock;
	clock_t t1;
	clock_t t2;

	//inialisasi
	awal_in = (float *)malloc(sizeof(float) * totalElemen);
	awal_out = (float *)malloc(sizeof(float) * totalElemen);
	awal_err = (float *)malloc(sizeof(float) * totalElemen);

	//instansiasi
	printf("generating domain matrix...");
	for(int i = 0; i<totalElemen; i++) {
		//kondisi batas
		if(i<SISI_MATRIKS) {
			awal_in[i] = 100;
		}
		else {
			awal_in[i] = 0;
		}
		awal_out[i] = awal_in[i];
		awal_err[i] = 0;
	}
	printf("done\n");
	
	printf("\nMulai pengujian...\n------------------------------------------\n");
	printf("CPU\n");
	//pengujian
	printf("Metode Jacobi.....");
	//---CPU Jacobi---
	//alokasi memori domain
	host_in = (float *)malloc(sizeof(float) * totalElemen);
	host_out = (float *)malloc(sizeof(float) * totalElemen);
	host_err = (float *)malloc(sizeof(float) * totalElemen);

	//copy data dari matriks awal ke domain
	for(int i = 0; i<totalElemen; i++) {
		host_in[i] = awal_in[i];
		host_out[i] = awal_out[i];
		host_err[i] = awal_err[i];
	}
	//mulai komputasi
	t1 = clock() / (CLOCKS_PER_SEC / 1000);
	total_error = -1;
	while ((total_error == -1) || (total_error > MAX_ERROR)) {
		cpuJacobi(host_in, host_out, host_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS, totalElemen);
		
		//hitung error
		total_error = 0;
		for(int i=0; i<totalElemen; i++) {
			total_error += absol(host_err[i]);
		}
		//printf("%d, %f\n",CPUIter,total_error);
		CPUIter++;
	}
	t2 = clock() / (CLOCKS_PER_SEC / 1000);
	CPUTime = t2-t1;

	writeToFile(SISI_MATRIKS, SISI_MATRIKS, totalElemen, "CPUJacobi.hasil", host_out);

	//bebaskan memori
	free(host_in);
	free(host_out);
	free(host_err);
	printf("done\n");
	printf("Banyak iterasi = %d\nWaktu komputasi = %f\n", CPUIter, CPUTime);
	
	//---CPU PGS---
	printf("\n");
	printf("Metode Point Gauss-Seidel.....");
	//alokasi memori domain
	host_in = (float *)malloc(sizeof(float) * totalElemen);
	host_out = (float *)malloc(sizeof(float) * totalElemen);
	host_err = (float *)malloc(sizeof(float) * totalElemen);

	//copy data dari matriks awal ke domain
	for(int i = 0; i<totalElemen; i++) {
		host_in[i] = awal_in[i];
		host_out[i] = awal_out[i];
		host_err[i] = awal_err[i];
	}
	//mulai komputasi
	t1 = clock() / (CLOCKS_PER_SEC / 1000);
	total_error = -1;
	while ((total_error == -1) || (total_error > MAX_ERROR)) {
		cpuPGS(host_in, host_out, host_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS, totalElemen);
		
		//hitung error
		total_error = 0;
		for(int i=0; i<totalElemen; i++) {
			total_error += absol(host_err[i]);
		}
		//printf("%d, %f\n",CPUIter,total_error);
		CPUIter2++;
	}
	t2 = clock() / (CLOCKS_PER_SEC / 1000);
	CPUTime2 = t2-t1;

	writeToFile(SISI_MATRIKS, SISI_MATRIKS, totalElemen, "CPUPGS.hasil", host_out);

	//bebaskan memori
	free(host_in);
	free(host_out);
	free(host_err);
	printf("done\n");
	printf("Banyak iterasi = %d\nWaktu komputasi = %f\n", CPUIter2, CPUTime2);

	//---CPU PSOR---
	printf("\n");
	printf("Metode Point Successive-Over-Relaxation.....");
	//alokasi memori domain
	host_in = (float *)malloc(sizeof(float) * totalElemen);
	host_out = (float *)malloc(sizeof(float) * totalElemen);
	host_err = (float *)malloc(sizeof(float) * totalElemen);

	//copy data dari matriks awal ke domain
	for(int i = 0; i<totalElemen; i++) {
		host_in[i] = awal_in[i];
		host_out[i] = awal_out[i];
		host_err[i] = awal_err[i];
	}
	//mulai komputasi
	t1 = clock() / (CLOCKS_PER_SEC / 1000);
	total_error = -1;
	while ((total_error == -1) || (total_error > MAX_ERROR)) {
		cpuPSOR(host_in, host_out, host_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS, totalElemen, omega);
		
		//hitung error
		total_error = 0;
		for(int i=0; i<totalElemen; i++) {
			total_error += absol(host_err[i]);
		}
		//printf("%d, %f\n",CPUIter,total_error);
		CPUIter3++;
	}
	t2 = clock() / (CLOCKS_PER_SEC / 1000);
	CPUTime3 = t2-t1;

	writeToFile(SISI_MATRIKS, SISI_MATRIKS, totalElemen, "CPUPSOR.hasil", host_out);

	//bebaskan memori
	free(host_in);
	//free(host_out);
	//free(host_err);
	printf("done\n");
	printf("Banyak iterasi = %d\nWaktu komputasi = %f\n", CPUIter3, CPUTime3);
	
	printf("\n----\n");

	printf("\nGPU: 1x1");
	//sebelum GPU:inisialisasi dimensi matriks
	int threadx = 1, thready = 1;
	int blockx, blocky;
	blockx = (SISI_MATRIKS/threadx)+1;
	blocky = (SISI_MATRIKS/thready)+1;
	jumlahBlock = dim3(blockx, blocky);
	threadPerBlock = dim3(threadx,thready);


	//---GPU Jacobi---
	printf("\n");
	printf("GPU Jacobi....");
	//alokasi memori domain
	cudaMalloc( (void **)&dev_in, sizeof(float) * totalElemen) ;
	cudaMalloc( (void **)&dev_err, sizeof(float) * totalElemen);

	//copy data dari matriks awal ke domain
	cudaMemcpy(dev_in, awal_in, totalElemen*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_err, awal_err, totalElemen*sizeof(float), cudaMemcpyHostToDevice);

	//mulai komputasi
	t1 = clock() / (CLOCKS_PER_SEC / 1000);
	total_error = -1;
	while ((total_error == -1) || (total_error > MAX_ERROR)) {
		gpuJacobi<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS);
		
		//hitung error
		total_error = 0;
		cudaMemcpy(host_err, dev_err, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);

		for(int i=0; i<totalElemen; i++) {
			total_error += absol(host_err[i]);
		}
		//printf("%d, %f\n",CPUIter,total_error);
		GPUIter++;
	}
	t2 = clock() / (CLOCKS_PER_SEC / 1000);
	GPUTime = t2-t1;

	cudaMemcpy(host_out, dev_in, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);
	writeToFile(SISI_MATRIKS, SISI_MATRIKS, totalElemen, "GPUJacobi1x1.hasil", host_out);

	//bebaskan memori
	cudaFree(dev_in);
	cudaFree(dev_err);
	printf("done\n");
	printf("Banyak iterasi = %d\nWaktu komputasi = %f\n", GPUIter, GPUTime);
	
	//---GPU PGS---
	printf("\n");
	printf("GPU PGS....");
	//alokasi memori domain
	cudaMalloc( (void **)&dev_in, sizeof(float) * totalElemen) ;
	cudaMalloc( (void **)&dev_err, sizeof(float) * totalElemen);

	//copy data dari matriks awal ke domain
	cudaMemcpy(dev_in, awal_in, totalElemen*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_err, awal_err, totalElemen*sizeof(float), cudaMemcpyHostToDevice);

	//mulai komputasi
	t1 = clock() / (CLOCKS_PER_SEC / 1000);
	total_error = -1;
	while ((total_error == -1) || (total_error > MAX_ERROR)) {
		gpuPGSRed<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS);
		gpuPGSBlack<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS);
		
		//hitung error
		total_error = 0;
		cudaMemcpy(host_err, dev_err, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);

		for(int i=0; i<totalElemen; i++) {
			total_error += absol(host_err[i]);
		}
		//printf("%d, %f\n",CPUIter,total_error);
		GPUIter2++;
	}
	t2 = clock() / (CLOCKS_PER_SEC / 1000);
	GPUTime2 = t2-t1;

	cudaMemcpy(host_out, dev_in, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);
	writeToFile(SISI_MATRIKS, SISI_MATRIKS, totalElemen, "GPUPGS1x1.hasil", host_out);

	//bebaskan memori
	cudaFree(dev_in);
	cudaFree(dev_err);
	printf("done\n");
	printf("Banyak iterasi = %d\nWaktu komputasi = %f\n", GPUIter2, GPUTime2);

	//---GPU PSOR---
	printf("\n");
	printf("GPU PSOR....");
	//alokasi memori domain
	cudaMalloc( (void **)&dev_in, sizeof(float) * totalElemen) ;
	cudaMalloc( (void **)&dev_err, sizeof(float) * totalElemen);

	//copy data dari matriks awal ke domain
	cudaMemcpy(dev_in, awal_in, totalElemen*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_err, awal_err, totalElemen*sizeof(float), cudaMemcpyHostToDevice);

	//mulai komputasi
	t1 = clock() / (CLOCKS_PER_SEC / 1000);
	total_error = -1;
	while ((total_error == -1) || (total_error > MAX_ERROR)) {
		gpuPSORRed<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS, omega);
		gpuPSORBlack<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS, omega);
		
		//hitung error
		total_error = 0;
		cudaMemcpy(host_err, dev_err, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);

		for(int i=0; i<totalElemen; i++) {
			total_error += absol(host_err[i]);
		}
		//printf("%d, %f\n",GPUIter3,total_error);
		GPUIter3++;
	}
	t2 = clock() / (CLOCKS_PER_SEC / 1000);
	GPUTime3 = t2-t1;

	cudaMemcpy(host_out, dev_in, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);
	writeToFile(SISI_MATRIKS, SISI_MATRIKS, totalElemen, "GPUPSOR1x1.hasil", host_out);

	//bebaskan memori
	cudaFree(dev_in);
	cudaFree(dev_err);
	printf("done\n");
	printf("Banyak iterasi = %d\nWaktu komputasi = %f\n", GPUIter3, GPUTime3);

	printf("\n----\n");

	printf("\nGPU: 2x2");
	//sebelum GPU:inisialisasi dimensi matriks
	threadx = 2;
	thready = 2;
	blockx = (SISI_MATRIKS/threadx)+1;
	blocky = (SISI_MATRIKS/thready)+1;
	jumlahBlock = dim3(blockx, blocky);
	threadPerBlock = dim3(threadx,thready);
	GPUIter = 0; GPUIter2 = 0; GPUIter3 = 0;
	GPUTime = 0; GPUTime2 = 0; GPUTime3 = 0;


	//---GPU Jacobi---
	printf("\n");
	printf("GPU Jacobi....");
	//alokasi memori domain
	cudaMalloc( (void **)&dev_in, sizeof(float) * totalElemen) ;
	cudaMalloc( (void **)&dev_err, sizeof(float) * totalElemen);

	//copy data dari matriks awal ke domain
	cudaMemcpy(dev_in, awal_in, totalElemen*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_err, awal_err, totalElemen*sizeof(float), cudaMemcpyHostToDevice);

	//mulai komputasi
	t1 = clock() / (CLOCKS_PER_SEC / 1000);
	total_error = -1;
	while ((total_error == -1) || (total_error > MAX_ERROR)) {
		gpuJacobi<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS);
		
		//hitung error
		total_error = 0;
		cudaMemcpy(host_err, dev_err, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);

		for(int i=0; i<totalElemen; i++) {
			total_error += absol(host_err[i]);
		}
		//printf("%d, %f\n",CPUIter,total_error);
		GPUIter++;
	}
	t2 = clock() / (CLOCKS_PER_SEC / 1000);
	GPUTime = t2-t1;

	cudaMemcpy(host_out, dev_in, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);
	writeToFile(SISI_MATRIKS, SISI_MATRIKS, totalElemen, "GPUJacobi2x2.hasil", host_out);

	//bebaskan memori
	cudaFree(dev_in);
	cudaFree(dev_err);
	printf("done\n");
	printf("Banyak iterasi = %d\nWaktu komputasi = %f\n", GPUIter, GPUTime);
	
	//---GPU PGS---
	printf("\n");
	printf("GPU PGS....");
	//alokasi memori domain
	cudaMalloc( (void **)&dev_in, sizeof(float) * totalElemen) ;
	cudaMalloc( (void **)&dev_err, sizeof(float) * totalElemen);

	//copy data dari matriks awal ke domain
	cudaMemcpy(dev_in, awal_in, totalElemen*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_err, awal_err, totalElemen*sizeof(float), cudaMemcpyHostToDevice);

	//mulai komputasi
	t1 = clock() / (CLOCKS_PER_SEC / 1000);
	total_error = -1;
	while ((total_error == -1) || (total_error > MAX_ERROR)) {
		gpuPGSRed<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS);
		gpuPGSBlack<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS);
		
		//hitung error
		total_error = 0;
		cudaMemcpy(host_err, dev_err, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);

		for(int i=0; i<totalElemen; i++) {
			total_error += absol(host_err[i]);
		}
		//printf("%d, %f\n",CPUIter,total_error);
		GPUIter2++;
	}
	t2 = clock() / (CLOCKS_PER_SEC / 1000);
	GPUTime2 = t2-t1;

	cudaMemcpy(host_out, dev_in, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);
	writeToFile(SISI_MATRIKS, SISI_MATRIKS, totalElemen, "GPUPGS2x2.hasil", host_out);

	//bebaskan memori
	cudaFree(dev_in);
	cudaFree(dev_err);
	printf("done\n");
	printf("Banyak iterasi = %d\nWaktu komputasi = %f\n", GPUIter2, GPUTime2);

	//---GPU PSOR---
	printf("\n");
	printf("GPU PSOR....");
	//alokasi memori domain
	cudaMalloc( (void **)&dev_in, sizeof(float) * totalElemen) ;
	cudaMalloc( (void **)&dev_err, sizeof(float) * totalElemen);

	//copy data dari matriks awal ke domain
	cudaMemcpy(dev_in, awal_in, totalElemen*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_err, awal_err, totalElemen*sizeof(float), cudaMemcpyHostToDevice);

	//mulai komputasi
	t1 = clock() / (CLOCKS_PER_SEC / 1000);
	total_error = -1;
	while ((total_error == -1) || (total_error > MAX_ERROR)) {
		gpuPSORRed<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS, omega);
		gpuPSORBlack<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS, omega);
		
		//hitung error
		total_error = 0;
		cudaMemcpy(host_err, dev_err, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);

		for(int i=0; i<totalElemen; i++) {
			total_error += absol(host_err[i]);
		}
		//printf("%d, %f\n",CPUIter,total_error);
		GPUIter3++;
	}
	t2 = clock() / (CLOCKS_PER_SEC / 1000);
	GPUTime3 = t2-t1;

	cudaMemcpy(host_out, dev_in, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);
	writeToFile(SISI_MATRIKS, SISI_MATRIKS, totalElemen, "GPUPSOR2x2.hasil", host_out);

	//bebaskan memori
	cudaFree(dev_in);
	cudaFree(dev_err);
	printf("done\n");
	printf("Banyak iterasi = %d\nWaktu komputasi = %f\n", GPUIter3, GPUTime3);

	printf("\n----\n");

	printf("\nGPU: 4x4");
	//sebelum GPU:inisialisasi dimensi matriks
	threadx = 4;
	thready = 4;
	blockx = (SISI_MATRIKS/threadx)+1;
	blocky = (SISI_MATRIKS/thready)+1;
	jumlahBlock = dim3(blockx, blocky);
	threadPerBlock = dim3(threadx,thready);
	GPUIter = 0; GPUIter2 = 0; GPUIter3 = 0;
	GPUTime = 0; GPUTime2 = 0; GPUTime3 = 0;


	//---GPU Jacobi---
	printf("\n");
	printf("GPU Jacobi....");
	//alokasi memori domain
	cudaMalloc( (void **)&dev_in, sizeof(float) * totalElemen) ;
	cudaMalloc( (void **)&dev_err, sizeof(float) * totalElemen);

	//copy data dari matriks awal ke domain
	cudaMemcpy(dev_in, awal_in, totalElemen*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_err, awal_err, totalElemen*sizeof(float), cudaMemcpyHostToDevice);

	//mulai komputasi
	t1 = clock() / (CLOCKS_PER_SEC / 1000);
	total_error = -1;
	while ((total_error == -1) || (total_error > MAX_ERROR)) {
		gpuJacobi<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS);
		
		//hitung error
		total_error = 0;
		cudaMemcpy(host_err, dev_err, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);

		for(int i=0; i<totalElemen; i++) {
			total_error += absol(host_err[i]);
		}
		//printf("%d, %f\n",CPUIter,total_error);
		GPUIter++;
	}
	t2 = clock() / (CLOCKS_PER_SEC / 1000);
	GPUTime = t2-t1;

	cudaMemcpy(host_out, dev_in, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);
	writeToFile(SISI_MATRIKS, SISI_MATRIKS, totalElemen, "GPUJacobi4x4.hasil", host_out);

	//bebaskan memori
	cudaFree(dev_in);
	cudaFree(dev_err);
	printf("done\n");
	printf("Banyak iterasi = %d\nWaktu komputasi = %f\n", GPUIter, GPUTime);
	
	//---GPU PGS---
	printf("\n");
	printf("GPU PGS....");
	//alokasi memori domain
	cudaMalloc( (void **)&dev_in, sizeof(float) * totalElemen) ;
	cudaMalloc( (void **)&dev_err, sizeof(float) * totalElemen);

	//copy data dari matriks awal ke domain
	cudaMemcpy(dev_in, awal_in, totalElemen*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_err, awal_err, totalElemen*sizeof(float), cudaMemcpyHostToDevice);

	//mulai komputasi
	t1 = clock() / (CLOCKS_PER_SEC / 1000);
	total_error = -1;
	while ((total_error == -1) || (total_error > MAX_ERROR)) {
		gpuPGSRed<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS);
		gpuPGSBlack<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS);
		
		//hitung error
		total_error = 0;
		cudaMemcpy(host_err, dev_err, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);

		for(int i=0; i<totalElemen; i++) {
			total_error += absol(host_err[i]);
		}
		//printf("%d, %f\n",CPUIter,total_error);
		GPUIter2++;
	}
	t2 = clock() / (CLOCKS_PER_SEC / 1000);
	GPUTime2 = t2-t1;

	cudaMemcpy(host_out, dev_in, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);
	writeToFile(SISI_MATRIKS, SISI_MATRIKS, totalElemen, "GPUPGS4x4.hasil", host_out);

	//bebaskan memori
	cudaFree(dev_in);
	cudaFree(dev_err);
	printf("done\n");
	printf("Banyak iterasi = %d\nWaktu komputasi = %f\n", GPUIter2, GPUTime2);

	//---GPU PSOR---
	printf("\n");
	printf("GPU PSOR....");
	//alokasi memori domain
	cudaMalloc( (void **)&dev_in, sizeof(float) * totalElemen) ;
	cudaMalloc( (void **)&dev_err, sizeof(float) * totalElemen);

	//copy data dari matriks awal ke domain
	cudaMemcpy(dev_in, awal_in, totalElemen*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_err, awal_err, totalElemen*sizeof(float), cudaMemcpyHostToDevice);

	//mulai komputasi
	t1 = clock() / (CLOCKS_PER_SEC / 1000);
	total_error = -1;
	while ((total_error == -1) || (total_error > MAX_ERROR)) {
		gpuPSORRed<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS, omega);
		gpuPSORBlack<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS, omega);
		
		//hitung error
		total_error = 0;
		cudaMemcpy(host_err, dev_err, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);

		for(int i=0; i<totalElemen; i++) {
			total_error += absol(host_err[i]);
		}
		//printf("%d, %f\n",CPUIter,total_error);
		GPUIter3++;
	}
	t2 = clock() / (CLOCKS_PER_SEC / 1000);
	GPUTime3 = t2-t1;

	cudaMemcpy(host_out, dev_in, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);
	writeToFile(SISI_MATRIKS, SISI_MATRIKS, totalElemen, "GPUPSOR4x4.hasil", host_out);

	//bebaskan memori
	cudaFree(dev_in);
	cudaFree(dev_err);
	printf("done\n");
	printf("Banyak iterasi = %d\nWaktu komputasi = %f\n", GPUIter3, GPUTime3);
	
	printf("\n----\n");

	printf("\nGPU: 8x8");
	//sebelum GPU:inisialisasi dimensi matriks
	threadx = 8;
	thready = 8;
	blockx = (SISI_MATRIKS/threadx)+1;
	blocky = (SISI_MATRIKS/thready)+1;
	jumlahBlock = dim3(blockx, blocky);
	threadPerBlock = dim3(threadx,thready);
	GPUIter = 0; GPUIter2 = 0; GPUIter3 = 0;
	GPUTime = 0; GPUTime2 = 0; GPUTime3 = 0;


	//---GPU Jacobi---
	printf("\n");
	printf("GPU Jacobi....");
	//alokasi memori domain
	cudaMalloc( (void **)&dev_in, sizeof(float) * totalElemen) ;
	cudaMalloc( (void **)&dev_err, sizeof(float) * totalElemen);

	//copy data dari matriks awal ke domain
	cudaMemcpy(dev_in, awal_in, totalElemen*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_err, awal_err, totalElemen*sizeof(float), cudaMemcpyHostToDevice);

	//mulai komputasi
	t1 = clock() / (CLOCKS_PER_SEC / 1000);
	total_error = -1;
	while ((total_error == -1) || (total_error > MAX_ERROR)) {
		gpuJacobi<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS);
		
		//hitung error
		total_error = 0;
		cudaMemcpy(host_err, dev_err, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);

		for(int i=0; i<totalElemen; i++) {
			total_error += absol(host_err[i]);
		}
		//printf("%d, %f\n",CPUIter,total_error);
		GPUIter++;
	}
	t2 = clock() / (CLOCKS_PER_SEC / 1000);
	GPUTime = t2-t1;

	cudaMemcpy(host_out, dev_in, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);
	writeToFile(SISI_MATRIKS, SISI_MATRIKS, totalElemen, "GPUJacobi8x8.hasil", host_out);

	//bebaskan memori
	cudaFree(dev_in);
	cudaFree(dev_err);
	printf("done\n");
	printf("Banyak iterasi = %d\nWaktu komputasi = %f\n", GPUIter, GPUTime);
	
	//---GPU PGS---
	printf("\n");
	printf("GPU PGS....");
	//alokasi memori domain
	cudaMalloc( (void **)&dev_in, sizeof(float) * totalElemen) ;
	cudaMalloc( (void **)&dev_err, sizeof(float) * totalElemen);

	//copy data dari matriks awal ke domain
	cudaMemcpy(dev_in, awal_in, totalElemen*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_err, awal_err, totalElemen*sizeof(float), cudaMemcpyHostToDevice);

	//mulai komputasi
	t1 = clock() / (CLOCKS_PER_SEC / 1000);
	total_error = -1;
	while ((total_error == -1) || (total_error > MAX_ERROR)) {
		gpuPGSRed<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS);
		gpuPGSBlack<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS);
		
		//hitung error
		total_error = 0;
		cudaMemcpy(host_err, dev_err, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);

		for(int i=0; i<totalElemen; i++) {
			total_error += absol(host_err[i]);
		}
		//printf("%d, %f\n",CPUIter,total_error);
		GPUIter2++;
	}
	t2 = clock() / (CLOCKS_PER_SEC / 1000);
	GPUTime2 = t2-t1;

	cudaMemcpy(host_out, dev_in, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);
	writeToFile(SISI_MATRIKS, SISI_MATRIKS, totalElemen, "GPUPGS8x8.hasil", host_out);

	//bebaskan memori
	cudaFree(dev_in);
	cudaFree(dev_err);
	printf("done\n");
	printf("Banyak iterasi = %d\nWaktu komputasi = %f\n", GPUIter2, GPUTime2);

	//---GPU PSOR---
	printf("\n");
	printf("GPU PSOR....");
	//alokasi memori domain
	cudaMalloc( (void **)&dev_in, sizeof(float) * totalElemen) ;
	cudaMalloc( (void **)&dev_err, sizeof(float) * totalElemen);

	//copy data dari matriks awal ke domain
	cudaMemcpy(dev_in, awal_in, totalElemen*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_err, awal_err, totalElemen*sizeof(float), cudaMemcpyHostToDevice);

	//mulai komputasi
	t1 = clock() / (CLOCKS_PER_SEC / 1000);
	total_error = -1;
	while ((total_error == -1) || (total_error > MAX_ERROR)) {
		gpuPSORRed<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS, omega);
		gpuPSORBlack<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS, omega);
		
		//hitung error
		total_error = 0;
		cudaMemcpy(host_err, dev_err, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);

		for(int i=0; i<totalElemen; i++) {
			total_error += absol(host_err[i]);
		}
		//printf("%d, %f\n",CPUIter,total_error);
		GPUIter3++;
	}
	t2 = clock() / (CLOCKS_PER_SEC / 1000);
	GPUTime3 = t2-t1;
	
	cudaMemcpy(host_out, dev_in, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);
	writeToFile(SISI_MATRIKS, SISI_MATRIKS, totalElemen, "GPUPSOR8x8.hasil", host_out);

	//bebaskan memori
	cudaFree(dev_in);
	cudaFree(dev_err);
	printf("done\n");
	printf("Banyak iterasi = %d\nWaktu komputasi = %f\n", GPUIter3, GPUTime3);

	printf("\n----\n");

	printf("\nGPU: 16x16");
	//sebelum GPU:inisialisasi dimensi matriks
	threadx = 16;
	thready = 16;
	blockx = (SISI_MATRIKS/threadx)+1;
	blocky = (SISI_MATRIKS/thready)+1;
	jumlahBlock = dim3(blockx, blocky);
	threadPerBlock = dim3(threadx,thready);
	GPUIter = 0; GPUIter2 = 0; GPUIter3 = 0;
	GPUTime = 0; GPUTime2 = 0; GPUTime3 = 0;


	//---GPU Jacobi---
	printf("\n");
	printf("GPU Jacobi....");
	//alokasi memori domain
	cudaMalloc( (void **)&dev_in, sizeof(float) * totalElemen) ;
	cudaMalloc( (void **)&dev_err, sizeof(float) * totalElemen);

	//copy data dari matriks awal ke domain
	cudaMemcpy(dev_in, awal_in, totalElemen*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_err, awal_err, totalElemen*sizeof(float), cudaMemcpyHostToDevice);

	//mulai komputasi
	t1 = clock() / (CLOCKS_PER_SEC / 1000);
	total_error = -1;
	while ((total_error == -1) || (total_error > MAX_ERROR)) {
		gpuJacobi<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS);
		
		//hitung error
		total_error = 0;
		cudaMemcpy(host_err, dev_err, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);

		for(int i=0; i<totalElemen; i++) {
			total_error += absol(host_err[i]);
		}
		//printf("%d, %f\n",CPUIter,total_error);
		GPUIter++;
	}
	t2 = clock() / (CLOCKS_PER_SEC / 1000);
	GPUTime = t2-t1;

	cudaMemcpy(host_out, dev_in, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);
	writeToFile(SISI_MATRIKS, SISI_MATRIKS, totalElemen, "GPUJacobi16x16.hasil", host_out);

	//bebaskan memori
	cudaFree(dev_in);
	cudaFree(dev_err);
	printf("done\n");
	printf("Banyak iterasi = %d\nWaktu komputasi = %f\n", GPUIter, GPUTime);
	
	//---GPU PGS---
	printf("\n");
	printf("GPU PGS....");
	//alokasi memori domain
	cudaMalloc( (void **)&dev_in, sizeof(float) * totalElemen) ;
	cudaMalloc( (void **)&dev_err, sizeof(float) * totalElemen);

	//copy data dari matriks awal ke domain
	cudaMemcpy(dev_in, awal_in, totalElemen*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_err, awal_err, totalElemen*sizeof(float), cudaMemcpyHostToDevice);

	//mulai komputasi
	t1 = clock() / (CLOCKS_PER_SEC / 1000);
	total_error = -1;
	while ((total_error == -1) || (total_error > MAX_ERROR)) {
		gpuPGSRed<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS);
		gpuPGSBlack<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS);
		
		//hitung error
		total_error = 0;
		cudaMemcpy(host_err, dev_err, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);

		for(int i=0; i<totalElemen; i++) {
			total_error += absol(host_err[i]);
		}
		//printf("%d, %f\n",CPUIter,total_error);
		GPUIter2++;
	}
	t2 = clock() / (CLOCKS_PER_SEC / 1000);
	GPUTime2 = t2-t1;

	cudaMemcpy(host_out, dev_in, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);
	writeToFile(SISI_MATRIKS, SISI_MATRIKS, totalElemen, "GPUPGS16x16.hasil", host_out);

	//bebaskan memori
	cudaFree(dev_in);
	cudaFree(dev_err);
	printf("done\n");
	printf("Banyak iterasi = %d\nWaktu komputasi = %f\n", GPUIter2, GPUTime2);

	//---GPU PSOR---
	printf("\n");
	printf("GPU PSOR....");
	//alokasi memori domain
	cudaMalloc( (void **)&dev_in, sizeof(float) * totalElemen) ;
	cudaMalloc( (void **)&dev_err, sizeof(float) * totalElemen);

	//copy data dari matriks awal ke domain
	cudaMemcpy(dev_in, awal_in, totalElemen*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_err, awal_err, totalElemen*sizeof(float), cudaMemcpyHostToDevice);

	//mulai komputasi
	t1 = clock() / (CLOCKS_PER_SEC / 1000);
	total_error = -1;
	while ((total_error == -1) || (total_error > MAX_ERROR)) {
		gpuPSORRed<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS, omega);
		gpuPSORBlack<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS, omega);
		
		//hitung error
		total_error = 0;
		cudaMemcpy(host_err, dev_err, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);

		for(int i=0; i<totalElemen; i++) {
			total_error += absol(host_err[i]);
		}
		//printf("%d, %f\n",CPUIter,total_error);
		GPUIter3++;
	}
	t2 = clock() / (CLOCKS_PER_SEC / 1000);
	GPUTime3 = t2-t1;

	cudaMemcpy(host_out, dev_in, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);
	writeToFile(SISI_MATRIKS, SISI_MATRIKS, totalElemen, "GPUPSOR16x16.hasil", host_out);

	//bebaskan memori
	cudaFree(dev_in);
	cudaFree(dev_err);
	printf("done\n");
	printf("Banyak iterasi = %d\nWaktu komputasi = %f\n", GPUIter3, GPUTime3);

	printf("\n----\n");

	printf("\nGPU: 32x32");
	//sebelum GPU:inisialisasi dimensi matriks
	threadx = 32;
	thready = 32;
	blockx = (SISI_MATRIKS/threadx)+1;
	blocky = (SISI_MATRIKS/thready)+1;
	jumlahBlock = dim3(blockx, blocky);
	threadPerBlock = dim3(threadx,thready);
	GPUIter = 0; GPUIter2 = 0; GPUIter3 = 0;
	GPUTime = 0; GPUTime2 = 0; GPUTime3 = 0;


	//---GPU Jacobi---
	printf("\n");
	printf("GPU Jacobi....");
	//alokasi memori domain
	cudaMalloc( (void **)&dev_in, sizeof(float) * totalElemen) ;
	cudaMalloc( (void **)&dev_err, sizeof(float) * totalElemen);

	//copy data dari matriks awal ke domain
	cudaMemcpy(dev_in, awal_in, totalElemen*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_err, awal_err, totalElemen*sizeof(float), cudaMemcpyHostToDevice);

	//mulai komputasi
	t1 = clock() / (CLOCKS_PER_SEC / 1000);
	total_error = -1;
	while ((total_error == -1) || (total_error > MAX_ERROR)) {
		gpuJacobi<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS);
		
		//hitung error
		total_error = 0;
		cudaMemcpy(host_err, dev_err, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);

		for(int i=0; i<totalElemen; i++) {
			total_error += absol(host_err[i]);
		}
		//printf("%d, %f\n",CPUIter,total_error);
		GPUIter++;
	}
	t2 = clock() / (CLOCKS_PER_SEC / 1000);
	GPUTime = t2-t1;

	cudaMemcpy(host_out, dev_in, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);
	writeToFile(SISI_MATRIKS, SISI_MATRIKS, totalElemen, "GPUJacobi32x32.hasil", host_out);

	//bebaskan memori
	cudaFree(dev_in);
	cudaFree(dev_err);
	printf("done\n");
	printf("Banyak iterasi = %d\nWaktu komputasi = %f\n", GPUIter, GPUTime);
	
	//---GPU PGS---
	printf("\n");
	printf("GPU PGS....");
	//alokasi memori domain
	cudaMalloc( (void **)&dev_in, sizeof(float) * totalElemen) ;
	cudaMalloc( (void **)&dev_err, sizeof(float) * totalElemen);

	//copy data dari matriks awal ke domain
	cudaMemcpy(dev_in, awal_in, totalElemen*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_err, awal_err, totalElemen*sizeof(float), cudaMemcpyHostToDevice);

	//mulai komputasi
	t1 = clock() / (CLOCKS_PER_SEC / 1000);
	total_error = -1;
	while ((total_error == -1) || (total_error > MAX_ERROR)) {
		gpuPGSRed<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS);
		gpuPGSBlack<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS);
		
		//hitung error
		total_error = 0;
		cudaMemcpy(host_err, dev_err, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);

		for(int i=0; i<totalElemen; i++) {
			total_error += absol(host_err[i]);
		}
		//printf("%d, %f\n",CPUIter,total_error);
		GPUIter2++;
	}
	t2 = clock() / (CLOCKS_PER_SEC / 1000);
	GPUTime2 = t2-t1;

	cudaMemcpy(host_out, dev_in, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);
	writeToFile(SISI_MATRIKS, SISI_MATRIKS, totalElemen, "GPUPGS32x32.hasil", host_out);

	//bebaskan memori
	cudaFree(dev_in);
	cudaFree(dev_err);
	printf("done\n");
	printf("Banyak iterasi = %d\nWaktu komputasi = %f\n", GPUIter2, GPUTime2);

	//---GPU PSOR---
	printf("\n");
	printf("GPU PSOR....");
	//alokasi memori domain
	cudaMalloc( (void **)&dev_in, sizeof(float) * totalElemen) ;
	cudaMalloc( (void **)&dev_err, sizeof(float) * totalElemen);

	//copy data dari matriks awal ke domain
	cudaMemcpy(dev_in, awal_in, totalElemen*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_err, awal_err, totalElemen*sizeof(float), cudaMemcpyHostToDevice);

	//mulai komputasi
	t1 = clock() / (CLOCKS_PER_SEC / 1000);
	total_error = -1;
	while ((total_error == -1) || (total_error > MAX_ERROR)) {
		gpuPSORRed<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS, omega);
		gpuPSORBlack<<<jumlahBlock, threadPerBlock>>>(dev_in, dev_err, ad_n, ad_e, ad_s, ad_w, SISI_MATRIKS, omega);
		
		//hitung error
		total_error = 0;
		cudaMemcpy(host_err, dev_err, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);

		for(int i=0; i<totalElemen; i++) {
			total_error += absol(host_err[i]);
		}
		//printf("%d, %f\n",CPUIter,total_error);
		GPUIter3++;
	}
	t2 = clock() / (CLOCKS_PER_SEC / 1000);
	GPUTime3 = t2-t1;

	cudaMemcpy(host_out, dev_in, totalElemen*sizeof(float), cudaMemcpyDeviceToHost);
	writeToFile(SISI_MATRIKS, SISI_MATRIKS, totalElemen, "GPUPSOR32x32.hasil", host_out);

	//bebaskan memori
	cudaFree(dev_in);
	cudaFree(dev_err);
	printf("done\n");
	printf("Banyak iterasi = %d\nWaktu komputasi = %f\n", GPUIter3, GPUTime3);

	getchar();
	return 0;
}

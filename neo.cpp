#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#define MAX_STR 512
#define MAX 1.0e12


int main(void) {
	/**********************variable**********************/
	FILE* ipf;
	FILE* svf;
	int i, j, h;
	float **X; //data table
	float **center;//cluster center
	float *mean; //save cluster means 
	float **dist; //save distance 
	float **M; //cluster mean
	int **U;

	int *bestCenter;
	char **label; //save labels
	char **cLabel; //save center labels
	int it, nRow, nCol;
	int k = 3; //cluster number

	/**********************file open**********************/
	ipf = fopen("input.txt", "r");
	if (ipf == NULL) {
		printf("Input file open error\n");
		exit(1);
	}
	svf = fopen("output.txt", "w");
	if (ipf == NULL) {
		printf("Output file open error\n");
		exit(1);
	}

	if (fscanf(ipf, "rows=%d columns=%d\n", &nRow, &nCol) != 2)
	{ //write rows=x columns=y at first line of file
		printf("Format error in first line\n");
		exit(1);
	}

	/**********************make matrix**********************/
	X = (float**)calloc(nRow, sizeof(float*));
	dist = (float**)calloc(nRow, sizeof(float*));	
	label = (char**)calloc(nRow, sizeof(char*));
	bestCenter = (int*)calloc(nRow, sizeof(int));
	U = (int**)calloc(nRow, sizeof(int*));

	center = (float**)calloc(k, sizeof(float*));
	cLabel = (char**)calloc(k, sizeof(char*));

	mean = (float*)calloc(nCol, sizeof(float));
	M = (float**)calloc(nCol, sizeof(float*));


	
	for (i = 0; i < nRow; i++) {
		X[i] = (float*)calloc(nCol, sizeof(float));
		label[i] = (char*)calloc(MAX_STR, sizeof(char));
		dist[i] = (float*)calloc(k, sizeof(float));
		U[i] = (int*)calloc(k, sizeof(int));
	}

	for (i = 0; i < k; i++) {
		center[i] = (float*)calloc(nCol, sizeof(float));
		cLabel[i] = (char*)calloc(MAX_STR, sizeof(char));
	}

	for (i = 0; i < nCol; i++) {
		M[i] = (float*)calloc(k, sizeof(float));
	}

	/**********************read labels and data**********************/
	for (i = 0; i < nRow; i++) {
		fscanf(ipf, "%s", label[i]);
		for (j = 0; j < nCol; j++) {
			fscanf(ipf, "%f", &(X[i][j]));
			mean[j] += X[i][j];
		}
		fscanf(ipf, "\n");

	}

	for (i = 0; i < nRow; i++) {
		for (j = 0; j < nCol; j++) {
			printf("%f ", X[i][j]);
		}
		printf("\n");

	}

	/**********************choose centers randomly**********************/
	int* a = (int*)calloc(k, sizeof(int));
	int temp; i = 0;
	while (i < k) {
		temp = (rand() % nRow) + 1;
		for (j = i - 1; j >= 0; j--) {
			if (a[j] == temp) {
				temp = -1;
			}
		}

		if (temp > 0) {
			a[i] = temp;
			i++;
		}
		
	}
	for (i = 0; i < k; i++) {
		printf("first center%d\n", a[i]);
		cLabel[i] = label[a[i]];
		center[i] = X[a[i]];
	}
	//cluster mean

	for (i = 0; i < nCol; i++) {
		mean[i] /= nRow;
	}


	int t = 0;
	int t_max = 100;
	int J = MAX;
	int oldJ = 0;
	int epsilon = 0;


	/**********************iteration**********************/
	while (t <= t_max) {

		oldJ = J;
		J = 0;
		//compute cluster means
		float sum; int num;
		for (i = 0; i < k; i++) {
			for (j = 0; j < nCol; j++) {
				sum = 0;
				num = 0;
				for (h = 0; h < nRow; h++) {
					if (U[h][i] == 1) {
						sum += X[h][j];
						num++;
					}	
				} //i번째 클러스터에서 X의 한 col 더함
				sum = sum / (float)num;
				M[j][i] = sum;
			}
		}


		//calc distance
		float rMin;
		for (i = 0; i < nRow; i++) {
			rMin = MAX;
			for (j = 0; j < k; j++) {
				dist[i][j] = 0;
				for (int e = 0; e < nCol; e++) {
					dist[i][j] += (X[i][e] - center[j][e])*(X[i][e] - center[j][e]);
				}
				if (dist[i][j] < rMin) {
					bestCenter[i] = j;
					rMin = dist[i][j];
				}
			}
		}





		//update cluster
		for (j = 0; i < k; j++) {
			float *cent_ = (float*)calloc(nCol, sizeof(float));
			int total = 0;
			for (i = 0; i < nRow; i++) {
				if (bestCenter[i] == j) {
					for (int e = 0; e < nCol; e++) {
						cent_[e] += X[i][e];
					}
					total++;
				}
			}

			if (total > 0) {
				for (int e = 0; e < nCol; e++) {
					center[j][e] = cent_[e] / total;
				}
			}
			else {
				for (int e = 0; e < nCol; e++) {
					center[j][e] = mean[e];
				}
			}
		}






		t = t + 1;

	}


	//file print
	for (i = 0; i < k; i++) {
		for (j = 0; j < nCol; j++) {
			fprintf(svf, "%f ", center[i][j]);
		}
		fprintf(svf, "\n");
	}

	for (i = 0; i < k; i++) {
		fprintf(svf, "Cluster %d\n", k);
		for (j = 0; j < nRow; j++) {
			if (bestCenter[j] == i) {
				fprintf(svf, "%s ", label[j]);
				for (h = 0; h < nCol; h++) {
					fprintf(svf, "%f ", X[j][h]);
				}
				fprintf(svf, "\n");
			}
		}

	}

	//cmd print
	for (j = 0; j < k; j++) {
		printf("\nCluster%d\n===========\n", j);
		for (i = 0; i < nRow; i++) {
			if (bestCenter[i] == j)
				printf("%s\n", label[i]);
		}

	}



}
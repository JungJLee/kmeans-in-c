#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#define MAX_STR 512
#define MAX 1.0e12


int main(void) {

	FILE* ipf;
	FILE* svf;
	int i, j;
	int h;
	float **X; //data table
	float **center;//cluster center
	float *mean; //save cluster means 
	float **dist; //save distance 

	int *bestCenter;
	char **label; //save labels
	char **cLabel; //save center labels
	int it, nRow, nCol;
	int k = 3; //cluster number

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

	//make matrices
	X = (float**)calloc(nRow, sizeof(float*));
	mean = (float*)calloc(nCol, sizeof(float));
	dist = (float**)calloc(nRow, sizeof(float*));
	center = (float**)calloc(k, sizeof(float*));
	label = (char**)calloc(nRow, sizeof(char*));	
	cLabel = (char**)calloc(k, sizeof(char*));
	bestCenter = (int*)calloc(nRow, sizeof(int));

	for (i = 0; i < nRow; i++) {
		X[i] = (float*)calloc(nCol, sizeof(float));
		label[i] = (char*)calloc(MAX_STR, sizeof(char));
		dist[i] = (float*)calloc(k, sizeof(float));
	}

	for (i = 0; i < k; i++) {
		center[i] = (float*)calloc(nCol, sizeof(float));
		cLabel[i] = (char*)calloc(MAX_STR, sizeof(char));
	}

	h = nRow;

	//read labels and data
	for (i = 0; i < nRow; i++) {
		fscanf(ipf, "%s", label[i]);
		for ( j = 0; j < nCol; j++) {
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

	//choose centers
	int* a = (int*)calloc(k, sizeof(int));
	int temp; i = 0;
	while (i < k) {
		temp = (rand() % nRow)+1;
		for (j = i-1; j >= 0; j--) {
			if (a[j] == temp) {
				temp = -1;
			}
		}

		if (temp > 0) {
			a[i] = temp;
			i++;
		}
		//
	}
	for (i = 0; i < k; i++) {
		printf("first center%d\n", a[i]);
		cLabel[i] = label[a[i]];
		center[i] = X[a[i]];
	}


	

	//n = number of data point
	//dim = dimension
	int t = 0;
	int t_max = 100;
	//int J = MAX;
	//int oldJ = 0;
	//int epsilon = 0;


		//cluster mean
		for (i = 0; i < nCol; i++) {
			mean[i] /= nRow;
		}

	while ( t <= t_max) {

		//calc distance
		float rMin;
		int count = 0;

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
			count++;
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



	for (j = 0; j < k; j++) {
		printf("\nCluster%d\n===========\n", j);
		for (i = 0; i < nRow; i++) {
			if (bestCenter[i] == j)
				printf("%s\n", label[i]);
		}

	}



}
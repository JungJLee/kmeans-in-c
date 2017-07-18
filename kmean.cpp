#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#define MAX_STR 512
#define MAX 2147483647


int main(void) {

	FILE* fp;
	int i, j;
	float **X;
	float **center;//cluster center
	float *mean; //cluster means save
	float **dist; //distance save

	char **label;
	char *cLabel;
	int it, nRow, nCol, nCol2, nCent;


	fp = fopen("input.txt", "r");
	if (fp == NULL) {
		perror("Eror opening file");
	}

	if (fscanf(fp, "rows=%d columns=%d\n", &nRow, &nCol) != 2)
	{
		printf("Format error in first line\n");
		exit(1);
	}

	X = (float**)calloc(nRow, sizeof(float*));

	mean = (float*)calloc(nCol, sizeof(float));
	dist = (float**)calloc(nRow, sizeof(float*));
	label = (char**)calloc(nRow, sizeof(char*));

	for (j = 0; j < nRow; j++) {
		X[j] = (float*)calloc(nCol, sizeof(float));
		label[j] = (char*)calloc(MAX_STR, sizeof(char));
	}
	i = nRow;

	//read labels




	//n = number of data point
	//dim = dimension
	//k = cluster number
	int t = 0;
	int t_max = 100;

	//int J = MAX;
	//int oldJ = 0;
	//int epsilon = 0;



	while ( t <= t_max) {

		//cluster mean

		//distance

		//update cluster

		t = t + 1;

	}



}
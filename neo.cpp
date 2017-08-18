#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#define MAX_STR 512
#define MAX 1.0e12
#define CLUSTER 5
#define R 5
#define C 158


int main(void) {
	time_t startTime = 0, endTime = 0;
	float gap;
	startTime = clock();


	/**********************cluster variable**********************/
	int k = CLUSTER; //cluster number
	float alpha = 0.1, beta = 0.005;
	int nRow, nCol;

	/**********************variables**********************/
	FILE* ipf;
	FILE* svf;
	FILE* svnf;
	int i, j, h;

	float **X; //data table : n x dim
	char **label; //save labels

	float **center;//cluster centers : k x dim
	float **oldcenter;
	float *Rmean; //save row means : dim
	int *cluNum; //final cluster number : n

	float **D; //save distance  : n x k
	float **M; //cluster mean : dim x k
	int **U; //cluster inform : n x k



	int it; // nRow = n , nCol = dim

	
					 
	/**********************file open******************** **/
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
	svnf = fopen("neo_output.txt", "w");
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
	D = (float**)calloc(nRow, sizeof(float*));
	label = (char**)calloc(nRow, sizeof(char*));
	cluNum = (int*)calloc(nRow, sizeof(int));
	U = (int**)calloc(nRow, sizeof(int*));

	center = (float**)calloc(k, sizeof(float*));
	oldcenter = (float**)calloc(k, sizeof(float*));

	Rmean = (float*)calloc(nCol, sizeof(float));
	M = (float**)calloc(nCol, sizeof(float*));

	for (i = 0; i < nRow; i++) {
		X[i] = (float*)calloc(nCol, sizeof(float));
		label[i] = (char*)calloc(MAX_STR, sizeof(char));
		D[i] = (float*)calloc(k, sizeof(float));
		U[i] = (int*)calloc(k, sizeof(int));
	}

	for (i = 0; i < k; i++) {
		center[i] = (float*)calloc(nCol, sizeof(float));
		oldcenter[i] = (float*)calloc(nCol, sizeof(float));
	}

	for (i = 0; i < nCol; i++) {
		M[i] = (float*)calloc(k, sizeof(float));
	}

	/**********************read labels and data**********************/
	for (i = 0; i < nRow; i++) {
		
		//lable 개수만큼
		fscanf(ipf, "%s", label[i]);
		fscanf(ipf, "%s", label[i]);
		
		for ( j = 0; j < nCol; j++) {
			fscanf(ipf, "%f", &(X[i][j]));
			Rmean[j] += X[i][j];
		}
		fscanf(ipf, "\n");

	}

	//row mean
	for (i = 0; i < nCol; i++) {
		Rmean[i] /= nRow;
	}

	/**********************choose centers randomly**********************/
	int* a = (int*)calloc(k, sizeof(int));
	int temp; 
	i = 0;
	srand((unsigned)time(NULL));
	while (i < k) {
		temp = (rand() % nRow);
		for (j = i-1; j >= 0; j--) {
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

		for (j = 0; j < nCol; j++) {;
		center[i][j] = X[a[i]][j];
		}
	}

	int t = 0;
	int t_max = 100;

	/**********************iteration**********************/


	float rMin;

	while ( t <= t_max ) {
		for (i = 0; i < k; i++) {
			for (j = 0; j < nCol; j++) {
				oldcenter[i][j] = center[i][j];
			}
		}

		//calc distance
		for (i = 0; i < nRow; i++) {
			rMin = MAX;
			for (j = 0; j < k; j++) {
				D[i][j] = 0;
				for (h = 0; h < nCol; h++) {
					D[i][j] += (X[i][h] - center[j][h])*(X[i][h] - center[j][h]);
				}
				if (D[i][j] < rMin) {
					cluNum[i] = j;
					rMin = D[i][j];
				}
			}
		}

		//update cluster
		float *cent_ = (float*)calloc(nCol, sizeof(float));
		for (j = 0; j < k; j++) {
			int total = 0;
			for (i = 0; i < nCol; i++) {
				cent_[i] = 0;
			}
			for (i = 0; i < nRow; i++) {
				if (cluNum[i] == j) {
					for (h = 0; h < nCol; h++) {
						cent_[h] += X[i][h];
					}
					total++;
				}
			}

			if (total > 0) {
				for (h = 0 ;h < nCol; h++) {
					center[j][h] = cent_[h] / total;
				}
			}
			else {
				for (h = 0 ; h < nCol; h++) {
					center[j][h] = Rmean[h];
					
				}
			}
		}


		//test
		int flag = 0;
		for (i = 0; i < k; i++) {
			for (j = 0; j < nCol; j++) {
				if (oldcenter[i][j] != center[i][j]) {
					flag++;
				}
			}			
		}
		if ( t>9 && flag == 0) break;

		t = t + 1;
	}//end of iteration

	for (i = 0; i < nRow; i++) {
			U[i][cluNum[i]] = 1;
	}

	endTime = clock();
	gap = (float)(endTime - startTime) / (CLOCKS_PER_SEC);
	printf("\n%f %d \n\n", gap,t);

	/**********************result print**********************/
	for (i = 0; i < k; i++) {
		for (j = 0; j < nCol; j++) {
			fprintf(svf, "%f ", center[i][j]);
		}
		fprintf(svf, "\n");
	}

	for (i = 0; i < k; i++) {
		fprintf(svf, "**********************Cluster %d\n", i);
		for (j = 0; j < nRow; j++) {
			if (cluNum[j] == i) {
				fprintf(svf, "%d ", j);
				for (h = 0; h < nCol; h++) {
					fprintf(svf, "%f ", X[j][h]);
				}
				fprintf(svf, "\n");
			}
		}

	}

	free(center);
	free(Rmean);
	free(cluNum);
	free(label);








	/*************************************************************************************
	******************************Non-exhaustiveness, overlapping*************************
	**************************************************************************************/

	float J = INFINITY;
	float oldJ = 0;
	int epsilon = 0;
	
	float N = (float)nRow;
	float alphaN = round(alpha*N);
	float betaN = round(beta*N);
	int nAssign;

	float** dnk = (float**)calloc(nRow, sizeof(float*));
	for (i = 0; i < nRow; i++) {
		dnk[i] = (float*)calloc(3, sizeof(float));
	}

	t = 0;
	
	float abj = (oldJ > J) ? (oldJ - J) : (J - oldJ);
	//================iteration
	while ((abj > epsilon) && (t <= t_max)) {
		oldJ = J;
		J = 0;
		//================compute cluster means
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
				if (num > 0) {
					sum = sum / (float)num;
				}
				M[j][i] = sum;

			}
		}


		//================compute distance

		float dif;
		for (i = 0; i < k; i++) {
			for (j = 0; j < nRow; j++) {
				dif = 0;
				for (h = 0; h < nCol; h++) {
					dif += (X[j][h] - M[h][i])*(X[j][h] - M[h][i]);
				}
				D[j][i] = dif;
			}
		}


		//================make N-betaN assignments

		float min;
		int mindx;
		for (i = 0; i < nRow; i++) {
			min = INFINITY;
			mindx = 0;
			for (j = 0; j < k; j++) {
				if (D[i][j] < min) {
					min = D[i][j];
					//printf("%f\n", min);
					mindx = j;
				}
			}
			dnk[i][0] = min;
			dnk[i][1] = (float)i;
			dnk[i][2] = (float)mindx;
		}


		float* dnkt = (float*)calloc(3, sizeof(float));
		for (i = 0; i < nRow - 1; i++) { //bubble sort
			for (j = 0; j < nRow - 1 - i; j++) {
				if (dnk[j][0] > dnk[j + 1][0]) {
					dnkt = dnk[j];
					dnk[j] = dnk[j + 1];
					dnk[j + 1] = dnkt;
				}
			}
		}


		nAssign = N - betaN;
		for (i = 0; i < nAssign; i++) {
			J = J + dnk[i][0];
		}


		int** tmp = (int**)calloc(nAssign, sizeof(int*));
		for (i = 0; i < nAssign; i++) {
			tmp[i] = (int*)calloc(2, sizeof(int));
		}

		for (i = 0; i < nRow; i++) {
			for (j = 0; j < k; j++) {
				U[i][j] = 0;
			}
		}
		int r, cl;
		for (i = 0; i < nAssign; i++) {
			r = dnk[i][1];
			cl = dnk[i][2];
			U[r][cl] = 1;
			tmp[i][0] = r;
			tmp[i][1] = cl;
		}


		for (i = 0; i < nAssign; i++) {
			D[tmp[i][0]][tmp[i][1]] = INFINITY;
		}

		//================make(alphaN + betaN) assignments
		int n = 0;
		float min_d;
		int i_star, j_star;
		while (n < alphaN+betaN) {
			min_d = INFINITY;
			for (i = 0; i < nRow; i++) {
				for (j = 0; j < k; j++) {
					if (D[i][j] < min_d) {
						min_d = D[i][j];
						i_star = i;
						j_star = j;
					}
				}
			}
			J = J + min_d;
			U[i_star][j_star] = 1;
			D[i_star][j_star] = INFINITY;

			n++;
		}

		t++;
		printf("Iteratoin %d, objective : %f\n", t, J);
		abj = (oldJ > J) ? (oldJ - J) : (J - oldJ);
		
	}


	for (i = 0; i < nRow; i++) {
		for (j = 0; j < k; j++) {
			fprintf(svnf, "%d ", U[i][j]);
		}
		fprintf(svnf, "\n");
	}




	fclose(ipf);
	fclose(svf);
	fclose(svnf);
	endTime = clock();
	gap = (float)(endTime - startTime) / (CLOCKS_PER_SEC);
	printf("%f", gap);
}


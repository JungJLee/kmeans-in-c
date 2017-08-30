#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#define MAX_STR 512
#define MAX 1.0e12
#define CLUSTER 3
float alpha;
float beta;


void quickSort(int left, int right, float** data);
void estimate_alpha_beta(float** X, float** C, float alpha_delta, float beta_delta, int k, int n, int dim);


int main(void) {
	time_t startTime = 0, endTime = 0;
	float gap;
	printf("k = %d\n", CLUSTER);
	/**********************cluster variable**********************/
	int k = CLUSTER; //cluster number

	//alpha, beta :set by user
	alpha = 0.1, beta = 0.005;
	// alpha_delta: (1) between -1 and 3.5 (the first strategy)
	// beta_delta: 3 or 6
	float alpha_delta = -0.1;
	float beta_delta = 6;

	int nk = 3;
	int nRow, nCol;// nRow = n , nCol = dim
	/**********************variables**********************/
	FILE* ipf;
	FILE* svf;
	FILE* svnf;
	int i, j, h;

	float **X; //data table : n x dim
	float **D; //save distance  : n x k

	char **label; //save labels
	float **center;//cluster centers : k x dim
	float **oldcenter;
	float *Rmean; //save row means : dim
	int *cluNum; //final cluster number : n


			 
	/**********************file open******************** **/
	ipf = fopen("input.txt", "r");
	if (ipf == NULL) {
		printf("Input file open error\n");
		exit(1);
	}
	//svf = fopen("output.txt", "w");
	//if (ipf == NULL) {
	//	printf("Output file open error\n");
	//	exit(1);
	//}
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


	center = (float**)calloc(k, sizeof(float*));
	oldcenter = (float**)calloc(k, sizeof(float*));

	Rmean = (float*)calloc(nCol, sizeof(float));


	for (i = 0; i < nRow; i++) {
		X[i] = (float*)calloc(nCol, sizeof(float));
		label[i] = (char*)calloc(MAX_STR, sizeof(char));
		D[i] = (float*)calloc(k, sizeof(float));

	}

	for (i = 0; i < k; i++) {
		center[i] = (float*)calloc(nCol, sizeof(float));
		oldcenter[i] = (float*)calloc(nCol, sizeof(float));
	}

	/**********************read labels and data**********************/
	for (i = 0; i < nRow; i++) {
		
		//lable 개수만큼 scan
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
	i = 0; j = 0;
	//srand((unsigned)time(NULL));
	while (i < k) {
		//temp = (rand() % nRow);
		//for (j = i-1; j >= 0; j--) {
		//	if (a[j] == temp) {
		//		temp = -1;
		//	}
		//}
		//if (temp > 0) {
		//	a[i] = temp;
		//	i++;
		//}
		a[i] = i;
		i++;

	}
	for (i = 0; i < k; i++) {

		for (j = 0; j < nCol; j++) {;
		center[i][j] = X[a[i]][j];
		}
	}

	int t = 0;
	int t_max = 100;


	startTime = clock();
	/*********************************************iteration*********************************************/


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


		// if the cluster == old center, exit while
		int flag = 0;
		for (i = 0; i < k; i++) {
			for (j = 0; j < nCol; j++) {
				if (oldcenter[i][j] != center[i][j]) {
					flag++;
				}
			}			
		}
		if ( flag == 0) break; 

		t = t + 1;
	}//end of iteration



	endTime = clock();
	gap = (float)(endTime - startTime) / (CLOCKS_PER_SEC);
	printf("\n%f %d \n", gap,t);

	/**********************result print**********************/
	//for (i = 0; i < k; i++) {
	//	for (j = 0; j < nCol; j++) {
	//		fprintf(svf, "%f ", center[i][j]);
	//	}
	//	fprintf(svf, "\n");
	//}

	//for (i = 0; i < k; i++) {
	//	fprintf(svf, "**********************Cluster %d\n", i);
	//	for (j = 0; j < nRow; j++) {
	//		if (cluNum[j] == i) {
	//			fprintf(svf, "%d ", j);
	//			for (h = 0; h < nCol; h++) {
	//				fprintf(svf, "%f ", X[j][h]);
	//			}
	//			fprintf(svf, "\n");
	//		}
	//	}

	//}

	estimate_alpha_beta(X, center, alpha_delta, beta_delta, k, nRow, nCol);
	printf("alpha : %f, beta : %f \n", alpha, beta);

	free(Rmean);
	free(label);
	free(center);
	free(oldcenter);



	int **U; //cluster inform : n x k
	U = (int**)calloc(nRow, sizeof(int*));
	for (i = 0; i < nRow; i++) {
		U[i] = (int*)calloc(nk, sizeof(int));
	}
	for (i = 0; i < nRow; i++) {
		U[i][0] = cluNum[i]+1; //cluster number : 1 ~ k
	}

	/*************************************************************************************
	******************************Non-exhaustiveness, overlapping*************************
	**************************************************************************************/

	double J = INFINITY;
	double oldJ = 0;
	int epsilon = 0;
	
	float N = (float)nRow;
	int alphaN = round(alpha*N);
	int betaN = round(beta*N);
	int nAssign = N - betaN;

	float **M; //cluster mean : dim x k
	M = (float**)calloc(nCol, sizeof(float*));
	for (i = 0; i < nCol; i++) {
		M[i] = (float*)calloc(k, sizeof(float));
	}

	float** dnk = (float**)calloc(nRow, sizeof(float*));
	for (i = 0; i < nRow; i++) {
		dnk[i] = (float*)calloc(3, sizeof(float));
	}

	t = 0;
	double abj = (oldJ > J) ? (oldJ - J) : (J - oldJ);

	
	//================iteration
	while ((abj > epsilon) && (t <= t_max)) {
		oldJ = J;
		J = 0;
		//================compute cluster means
		float sum;
		int num;
		int* clnum = (int*)calloc(k, sizeof(int));

		for (i = 0; i < nCol; i++) {
			for (j = 0; j < k; j++) {
				M[i][j] = 0;
			}
		}
		for (i = 0; i < nk; i++) {
			clnum[i] = 0;
		}

		for (i = 0; i < nRow; i++) {
			for (h = 0; h < nk; h++) {
				if (U[i][h] > 0) {
					for (j = 0; j < nCol; j++) {
						M[j][U[i][h] - 1] += X[i][j];
					}
					clnum[U[i][h] - 1] += 1;
				}
				else
					break;


			}
		}

		for (i = 0; i < nCol; i++) {
			for (j = 0; j < k; j++) {
				if (clnum[j] > 0) {
					M[i][j] = M[i][j] / clnum[j];
				}
			}
		}

		


		//for (i = 0; i < k; i++) {
		//	for (j = 0; j < nCol; j++) {
		//		sum = 0;
		//		num = 0;
		//		for (h = 0; h < nRow; h++) {
		//			for (int e = 0; e < nk; e++) {
		//				if (U[h][e] == 0) break;
		//				if (U[h][e] - 1 == i) {
		//					sum += X[h][j];
		//					num++;
		//				}
		//				if (e+1<nRow && U[h][e + 1] == 0) break;

		//			}
		//		} //ith cluster, sum of column at X
		//		if (num > 0) {
		//			sum = sum / (float)num;
		//		}
		//		M[j][i] = sum;

		//	}
		//}



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
					mindx = j;
				}
			}
			dnk[i][0] = min;
			dnk[i][1] = (float)i;
			dnk[i][2] = (float)mindx;
		}


		quickSort(0, nRow-1, dnk);


		for (i = 0; i < nAssign; i++) {
			J = J + dnk[i][0];
		}


		int** tmp = (int**)calloc(nAssign, sizeof(int*));
		for (i = 0; i < nAssign; i++) {
			tmp[i] = (int*)calloc(2, sizeof(int));
		}

		for (i = 0; i < nRow; i++) {
			for (j = 0; j < nk; j++) {
				U[i][j] = 0;
			}
		}
		int r, cl;
		for (i = 0; i < nAssign; i++) {
			r = (int)dnk[i][1];
			cl = (int)dnk[i][2];
			U[r][0] = cl + 1;

			D[r][cl] = INFINITY;
		}

		//for (i = 0; i < nRow; i++) {
		//	for (j = 0; j < k; j++) {
		//		printf("%d ", U[i][j]);
		//	}
		//	printf("\n");
		//}


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
			
			for (int e = 0; e < nk; e++) {
				if (U[i_star][e] == 0) {
					U[i_star][e] = j_star + 1;
					J = J + min_d;
					break;
				}
			}
			D[i_star][j_star] = INFINITY;

			n++;
		}

		t++;
		printf("Iteratoin %2d, objective : %f\n", t, J);
		abj = (oldJ > J) ? (oldJ - J) : (J - oldJ);
		
	}

	printf("alphaN = %d\n", alphaN);

	for (i = 0; i < nRow; i++) {
		for (j = 0; j < nk; j++) {
			if (j > 0 && U[i][j] == 0) {
				break;
			}
			fprintf(svnf, "%d ", U[i][j]);
		}

		fprintf(svnf, "\n");
	}




	fclose(ipf);
	//fclose(svf);
	fclose(svnf);
	endTime = clock();
	gap = (float)(endTime - startTime) / (CLOCKS_PER_SEC);
	printf("%f", gap);
}


void quickSort(int left, int right, float** data) {
	int pivot = left;
	int j = pivot;
	int i = left + 1;
	float* dnkt = (float*)calloc(3, sizeof(float));

	if (left < right) {
		for (i; i <= right; i++) {
			if (data[i][0] < data[pivot][0]) {
				j++;
				dnkt = data[j];
				data[j] = data[i];
				data[i] = dnkt;

			}
		}
		dnkt = data[left];
		data[left] = data[j];
		data[j] = dnkt;
		pivot = j;

		quickSort(left, pivot - 1, data);
		quickSort(pivot + 1, right, data);
	}

}

void estimate_alpha_beta(float** X, float** C, float alpha_delta, float beta_delta, int k, int n, int dim) {
	int i, j, h;
	float** D; //n*k
			   //X : n*dim C : k*dim 
	D = (float**)calloc(n, sizeof(float*));
	for (i = 0; i < n; i++) {
		D[i] = (float*)calloc(k, sizeof(float));
	}

	//compute distance
	float dif;
	for (i = 0; i < k; i++) {
		for (j = 0; j < n; j++) {
			dif = 0;
			for (h = 0; h < dim; h++) {
				dif += (X[j][h] - C[i][h])*(X[j][h] - C[i][h]);
			}
			D[j][i] = dif;
		}
	}


	//compute mean
	float** dist = (float**)calloc(n, sizeof(float*));
	for (i = 0; i < n; i++) {
		dist[i] = (float*)calloc(2, sizeof(float));
	}

	float min;
	int mindx;
	for (i = 0; i < n; i++) {
		min = INFINITY;
		mindx = 0;
		for (j = 0; j < k; j++) {
			if (D[i][j] < min) {
				min = D[i][j];
				mindx = j;
			}
		}
		dist[i][0] = min;
		dist[i][1] = (float)mindx;
	}

	//mean, standard dev
	float mean = 0, sum_dev = 0;
	for (i = 0; i < n; i++) {
		mean += dist[i][0];
	}
	mean = mean / n;
	for (i = 0; i < n; i++) {
		sum_dev += (dist[i][0] - mean)*(dist[i][0] - mean);
	}
	sum_dev = sqrt(sum_dev / n);

	//estimate beta
	int betaN = 0;
	float fnz = mean + beta_delta*sum_dev;
	for (i = 0; i < n; i++) {
		if (dist[i][0] > fnz)
			betaN++;
	}
	beta = (float)betaN / n;

	//estimate alpha
	int overlap = 0;
	float threshold;
	int v;
	int num;
	for (j = 0; j < k; j++) {
		mean = 0; num = 0; sum_dev = 0; v = 0;
		for (i = 0; i < n; i++) {
			if (dist[i][1] == j) {
				mean += dist[i][0];
				num++;
			}
		}
		if (num > 0) {
			mean = mean / num;
		}
		for (i = 0; i < n; i++) {
			if (dist[i][1] == j) {
				sum_dev += (dist[i][0] - mean)*(dist[i][0] - mean);
			}
		}
		if (num > 0) {
			sum_dev = sqrt(sum_dev / num);
		}
		threshold = mean + alpha_delta*sum_dev;
		for (i = 0; i < n; i++) {
			if (dist[i][1] != j) {
				overlap += (D[i][j] <= threshold);
			}

		}
	}//end of for
	alpha = (float)overlap / n;

	free(D);
	free(dist);
}
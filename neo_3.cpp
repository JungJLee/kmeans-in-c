#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#define MAX_STR 100
#define MAX 1.0e12

/**********************cluster option**********************/
int k; //cluster number
int lablenum; //number of lable column

//alpha, beta :set by user
float alpha = 0.1, beta = 0.005;
// alpha_delta: (1) between -1 and 3.5 (the first strategy)
// beta_delta: 3 or 6
float alpha_delta = -0.1;
float beta_delta = 6;
int alpha_beta_flag;
int iter_flag;
int norm_flag;

int nRow, nCol;// nRow = n , nCol = dim
int t_max = 100;

void quickSort(int left, int right, float** data, int colnumber);
void quickSortint(int left, int right, int** data);
void estimate_alpha_beta(float** X, float** C, float alpha_delta, float beta_delta, int k, int n, int dim);
void kmeans(float** X, float** D, int* cluNum, float* Rmean, float** center, float**oldcenter);

int main(void) {
	//option
	int selection = 1;
	//printf("Input 0 for default, 1 for select mode : ");
	//scanf("%d", &selection);
	if (selection == 1) {
		printf("Input cluster number : ");
		scanf("%d", &k);
		if (k < 2) {
			printf("cluster error\n");
			exit(1);
		}
		printf("Input lable column number : ");
		scanf("%d", &lablenum);
		printf("Choose alpha,beta. 0 for default, 1 for estimation : ");
		scanf("%d", &alpha_beta_flag);
		printf("Choose data normalization. 0 for default, 1 for normalization : ");
		scanf("%d", &norm_flag);
		printf("Choose iteration cutting. 0 for default, 1 for soft cut, 2 for improved : ");
		scanf("%d", &iter_flag);
	}
	else {
		k = 3;
		lablenum = 2;
		alpha_beta_flag = 0;
		norm_flag = 1;
		iter_flag = 2;
		printf("Choose iteration cutting. 0 for default, 1 for soft cut, 2 for improved : ");
		scanf("%d", &iter_flag);
	}




	time_t startTime = 0, endTime = 0;
	float gap;
	printf("k = %d\n", k);

	/**********************file open**********************/
	FILE* input_ptr;
	FILE* save_ptr;
	input_ptr = fopen("input.txt", "r");
	if (input_ptr == NULL) {
		printf("Input file open error\n");
		exit(1);
	}
	if (fscanf(input_ptr, "rows=%d columns=%d\n", &nRow, &nCol) != 2)
	{ //write rows=x columns=y at first line of file
		printf("Format error in first line\n");
		exit(1);
	}

	save_ptr = fopen("neo_output.txt", "w");
	if (save_ptr == NULL) {
		printf("Output file open error\n");
		exit(1);
	}

	/**********************variables**********************/
	int i, j, h;

	float **X; //data table : n x dim
	float **D; //save distance  : n x k
	char **label; //save labels
	float *Rmean; //save row means : dim
	int *cluNum; //final cluster number : n
	float **center;//cluster centers : k x dim
	float **oldcenter;


	center = (float**)calloc(k, sizeof(float*));
	oldcenter = (float**)calloc(k, sizeof(float*));

	X = (float**)calloc(nRow, sizeof(float*));
	D = (float**)calloc(nRow, sizeof(float*));
	label = (char**)calloc(nRow, sizeof(char*));
	cluNum = (int*)calloc(nRow, sizeof(int));
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


	/*****************read labels and data*****************/
	for (i = 0; i < nRow; i++) {
		//labelscan
		for (j = 0; j < lablenum ; j++) {
			fscanf(input_ptr, "%s", label[i]);
		}
		for (j = 0; j < nCol; j++) {
			fscanf(input_ptr, "%f", &(X[i][j]));
			Rmean[j] += X[i][j];
		}
		fscanf(input_ptr, "\n");
		
	}
	
	//normalization
	//ei = (ei-min)/(max-min)
	if (norm_flag == 1) {
		float minei = INFINITY, maxei = 0;
		for (i = 0; i < nRow; i++) {
			minei = INFINITY;
			maxei = 0;
			for (j = 0; j < nCol; j++) {
				if (X[i][j] > maxei)
					maxei = X[i][j];
				if (X[i][j] < minei)
					minei = X[i][j];
			}


			for (j = 0; j < nCol; j++) {
				X[i][j] = (X[i][j] - minei) / (maxei - minei);
			}
		}
	}


	/***************************************************************
	*******************************k-means**************************
	****************************************************************/
	
	startTime = clock();
	kmeans(X, D, cluNum, Rmean, center, oldcenter);
	endTime = clock();
	gap = (float)(endTime - startTime) / (CLOCKS_PER_SEC);
	printf("normal clustering time : %f \n", gap);

	if (alpha_beta_flag == 1) {
		estimate_alpha_beta(X, center, alpha_delta, beta_delta, k, nRow, nCol);
	}
	printf("alpha : %f, beta : %f \n", alpha, beta);

	free(Rmean);
	free(label);
	for (i = 0; i < k; i++) {
		free(center[i]);
		free(oldcenter[i]);
	}
	free(center);
	free(oldcenter);


	/******************************************************************
	********************Non-exhaustiveness, overlapping****************
	*******************************************************************/

	double J = INFINITY;
	double oldJ = 0;
	double tempJ = 0;
	int epsilon = 0;
	
	float N = (float)nRow;
	int alphaN = round(alpha*N);
	int betaN = round(beta*N);
	int nAssign = N - betaN;
	int uAssign = nRow + alphaN;

	int** U; //cluster inform : n+alpha * 2 
	float** M; //cluster mean : dim x k
	float** dnk; //for N-betaN assign

	U = (int**)calloc(uAssign, sizeof(int*));
	M = (float**)calloc(nCol, sizeof(float*));
	dnk = (float**)calloc(nRow, sizeof(float*));
	int* kNumlist = (int*)calloc(k, sizeof(int));

	for (i = 0; i < uAssign; i++) {
		U[i] = (int*)calloc(2, sizeof(int));
	}
	for (i = 0; i < nRow; i++) {
		U[i][0] = i;
		U[i][1] = cluNum[i]+1;
	}
	for (i = 0; i < nCol; i++) {
		M[i] = (float*)calloc(k, sizeof(float));
	}

	for (i = 0; i < nRow; i++) {
		dnk[i] = (float*)calloc(3, sizeof(float));
	}

	float** dab;
	int modcol = 100;
	if (iter_flag == 2) {
		dab = (float**)calloc(nRow*k, sizeof(float*));
		for (i = 0; i < nRow*k; i++) {
			dab[i] = (float*)calloc(2, sizeof(float));
		}
		while (k / modcol != 0) {
			modcol *= 10;
		}
		modcol *= 10;
	}

	printf("alphaN : %d betaN : %d\n", alphaN, betaN);
	int t = 0;
	int flag = 0;
	double abj = (oldJ > J) ? (oldJ - J) : (J - oldJ);
	//========================================iteration

	while ((abj > epsilon) && (t <= t_max)) {
		oldJ = J;
		J = 0;
		int a, b;
		//초기화
		for (i = 0; i < nCol; i++) {
			for (j = 0; j < k; j++) {
				M[i][j] = 0;
			}
		}
		for (i = 0; i < k; i++) {
			kNumlist[i] = 0;
		}

		for (i = 0; i < uAssign; i++) {
			if (U[i][1] > 0) {
				a = U[i][0];
				b = U[i][1] - 1;
				for (j = 0; j < nCol; j++) {
					M[j][b] += X[a][j];
				}
				kNumlist[b] += 1;
			}
			else
				break;

		}
		for (i = 0; i < nCol; i++) {
			for (j = 0; j < k; j++) {
				if (kNumlist[j] > 0) {
					M[i][j] = M[i][j] / kNumlist[j];
				}
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
		int r, cl;

		if (iter_flag == 0) { //기존
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

			quickSort(0, nRow - 1, dnk, 3);

			flag = 0;
			for (i = 0; i < nAssign; i++) {
				r = (int)dnk[i][1]; //index
				cl = (int)dnk[i][2]; //cluster
				J = J + dnk[i][0];

				cluNum[r] = cl + 1;
				U[i][0] = r;
				U[i][1] = cl + 1;
				D[r][cl] = INFINITY;
			}

		}
		else { //soft short
			if (flag < nAssign) {
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

				quickSort(0, nRow - 1, dnk, 3);

				flag = 0;
				for (i = 0; i < nAssign; i++) {
					r = (int)dnk[i][1]; //index
					cl = (int)dnk[i][2]; //cluster
					J = J + dnk[i][0];
					if (cluNum[r] == cl + 1) {
						flag++;
					}
					cluNum[r] = cl + 1;
					U[i][0] = r;
					U[i][1] = cl + 1;
					D[r][cl] = INFINITY;
				}
				if (flag == nAssign) {
					printf("%d\n", flag);
					tempJ = J;
				}
			}
			else {
				for (i = 0; i < nAssign; i++) {
					r = U[i][0];
					cl = U[i][1];
					J += D[r][cl - 1];
					D[r][cl - 1] = INFINITY;

				}
			}
		}

		h = nAssign;
		for (i = h; i < uAssign; i++) {
			for (j = 0; j < 2; j++) {
				U[i][j] = 0;
			}
		}
		//================make(alphaN + betaN) assignments
		int n = 0;
		float min_d;
		int i_star, j_star=0;
		if (iter_flag != 2) {
			while (n < alphaN + betaN) {
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
				U[h][0] = i_star;
				U[h][1] = j_star + 1;
				h++;
				D[i_star][j_star] = INFINITY;
				n++;
			}
		}
		else {
			min_d = INFINITY;
			for (i = 0; i < nRow; i++) {
				for (j = 0; j < k; j++) {
					if (D[i][j] != INFINITY) {
						dab[n][0] = D[i][j];
						dab[n][1] = (float)i;
						n++;
					}
				}
			}
			quickSort(0, n-1, dab, 2);
			n = 0;

			while (n < alphaN + betaN) {
				J = J + dab[n][0];
				U[h][0] = i_star = (int)dab[n][1];
				for ( j = 0; j < k; j++) {
					if (dab[n][0] == D[i_star][j]) {
						j_star = j;
						D[i_star][j] = INFINITY;
						break;
					}
				}
				U[h][1] = j_star + 1 ;
				h++;
				n++;

			}
		}


		t++;
		printf("Iteratoin %2d, objective : %f\n", t, J);
		abj = (oldJ > J) ? (oldJ - J) : (J - oldJ);
	}//end of iteration


	printf("alphaN : %d, %d\n", alphaN, nAssign);
	endTime = clock();
	gap = (float)(endTime - startTime) / (CLOCKS_PER_SEC);
	printf("NEO time : %f\n", gap);


	quickSortint(0, uAssign-1, U);

	j = 0;
	for (i = 0; i < uAssign; i++) {
		fprintf(save_ptr, "%d ", U[i][0]);
		fprintf(save_ptr, "%d\n", U[i][1]);
	}


	fclose(input_ptr);
	fclose(save_ptr);
}

void quickSortint(int left, int right, int** data) {
	int pivot = left;
	int j = pivot;
	int i = left + 1;
	int* dnkt = (int*)calloc(2, sizeof(int));

	if (left < right) {
		for (i; i <= right; i++) {
			if (data[i][0] < data[pivot][0]) {
				j++;
				dnkt = data[j];
				data[j] = data[i];
				data[i] = dnkt;

			}
			else if (data[i][0] == data[pivot][0]) {
				if (data[i][1] < data[pivot][1]) {
					j++;
					dnkt = data[j];
					data[j] = data[i];
					data[i] = dnkt;
				}
			}
		}
		dnkt = data[left];
		data[left] = data[j];
		data[j] = dnkt;
		pivot = j;

		quickSortint(left, pivot - 1, data);
		quickSortint(pivot + 1, right, data);
	}



}

float* dnkt3 = (float*)calloc(3, sizeof(float));
float* dnkt2 = (float*)calloc(3, sizeof(float));

void quickSort(int left, int right, float** data,int colnumber) {
	int pivot = left;
	int j = pivot;
	int i = left + 1;
	float* dnkt;
	if (colnumber == 3) {
		dnkt = dnkt3;
	}
	else if (colnumber == 2) {
		dnkt = dnkt2;
	}
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
		//free(dnkt);
		quickSort(left, pivot - 1, data, colnumber);
		quickSort(pivot + 1, right, data, colnumber);
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

void kmeans(float** X, float** D, int* cluNum, float* Rmean, float** center, float**oldcenter) {
	//FILE* svf;
	//svf = fopen("output.txt", "w");
	//if (ipf == NULL) {
	//	printf("Output file open error\n");
	//	exit(1);
	//}

	int i, j, h;
	int t = 0;

	//row mean
	for (i = 0; i < nCol; i++) {
		Rmean[i] /= nRow;
	}

	/**********************choose centers randomly**********************/
	int* a = (int*)calloc(k, sizeof(int));
	int temp;
	i = 0; j = 0;
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
		for (j = 0; j < nCol; j++) {
			center[i][j] = X[a[i]][j];
		}
	}

	/**********************start iteration**********************/
	float rMin;
	while (t <= t_max) {
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
				for (h = 0; h < nCol; h++) {
					center[j][h] = cent_[h] / total;
				}
			}
			else {
				for (h = 0; h < nCol; h++) {
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
		if (flag == 0) break;

		t = t + 1;
	}//end of iteration

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
	//fclose(svf);



}



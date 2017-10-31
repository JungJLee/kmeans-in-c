#pragma warning(disable:4996)
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#define ROW 22237
#define COL 3831
#define LAB 3831
#define MAX 24600

using namespace std;
int a[ROW+10][40]; //각 gene의 기능 목록 1~22237
int b[MAX][2]; //클러스터링 결과 정보 0~
int c[45][2]; //selected 0~
int clusize[260]; //클러스터링 결과에서 각 클러스터 사이즈

string label[LAB+10];
int labelnum[LAB+10];

void make_dense() {
	ifstream fp("3831function.txt");
	for (int i = 0; i < ROW; i++)
	{
		for (int j = 0; j < COL; j++)
		{
			fp >> a[i][j];
		}
	}

	ofstream fop("3831function_notsparse.txt");

	for (int i = 0; i < ROW; i++)
	{
		for (int j = 0; j < COL; j++)
		{
			if (a[i][j] != 0) {
				fop << j << ' ';
			}
		}
		fop << endl;
	}

fp.close();
fop.close();


}

void count_label() {
	ifstream fdp("3831function_notsparse.txt");

	int tempi = 0;
	while (!fdp.eof()) {
		fdp >> tempi;
		labelnum[tempi]++;
	}

	ofstream fdop("labelnumber.txt");
	for (int i = 0; i < LAB; i++) {
		fdop << labelnum[i] << endl;
	}

	fdp.close();
	fdop.close();

}

void make_abc() {
	char inputstring[1000];
	char* ptoken;
	char* key;
	int i = 0;
	int index = 1;
	//각 gene의 function 정보 a에 저장
	ifstream inf("3831function_notsparse.txt");
	while (!inf.eof()) {
		inf.getline(inputstring, 1000);

		if (inputstring[0] == '\0') {
			index++;
			continue;
		}
		i = 0;
		ptoken = strtok(inputstring, " ");
		key = ptoken;
		a[index][i++] = (atoi(key)) + 1;
		while (ptoken != NULL) {
			ptoken = strtok(NULL, " ");
			if (ptoken == NULL) {
				break;
			}
			key = ptoken;
			a[index][i++] = (atoi(key)) + 1;
		}
		index++;
	}
	inf.close();

	//클러스터링 결과물 불러옴. b에 저장
	FILE *cfp;
	if (cfp = fopen("256cluster.txt", "r")) {
		index = 0;
		while (!feof(cfp)) {
			fscanf(cfp, "%d %d\n", &b[index][0], &b[index][1]);
			index++;
		}
		fclose(cfp);
	}

	//selected 불러와 c에 저장
	FILE *sfp;
	if (sfp = fopen("selected.txt", "r")) {
		index = 0;
		while (!feof(sfp)) {
			fscanf(sfp, "%d %d\n", &c[index][0], &c[index][1]);
			index++;
		}
		fclose(sfp);
	}

}



void count_clu_mem() {
	for (int i = 0; i < MAX;i++) {
		clusize[b[i][1]]++;
	}

	ofstream svc("256clustercount.txt");
	for (int i = 1; i <= 256; i++) {
		svc << i << ' ' << clusize[i] << endl;
	}

	svc.close();

	//int sum = 0;
	//for (int i = 1; i <= 256; i++) {
	//	cout << i << ' ' << clusize[i] << endl;
	//	sum += clusize[i];
	//}
	//cout << sum << endl;

}

void search_by_labelnum() {

	//function number x를 찾음
	int x;
	int sum = 0;
	//cout << "input x" << endl;
	//cin >> x;
	x = c[1][0];
	vector <int> list;
	for (int i = 1; i <= ROW; i++) {
		int j = 0;
		while (a[i][j] != 0) {
			if (a[i][j] == x) {
				list.push_back(i);
			}
			j++;
		}
	}
	cout<<list.size()<<endl;
	vector<int> clulist;
	for (int i = 0; i < list.size(); i++) {
		for (int j = 0; j < MAX; j++) {
			if (b[j][0] == list[i]) {
				clulist.push_back(b[j][1]);
			}
			if (b[j][0] > list[i]) {
				break;
			}
		}

	}

	sort(clulist.begin(), clulist.end());
	int temp = clulist[0];
	int tempn = 1;
	for (int j = 1; j < clulist.size(); j++) {
		if (clulist[j] == temp) {
			tempn++;
		}
		else {
			cout << temp << " : \t" << tempn << endl;
			sum += tempn;
			tempn = 1;
			temp = clulist[j];
		}

	}

	cout << temp << " : \t" << tempn << endl;
	sum += tempn;
	tempn = 0;
	cout << sum<<endl;


}


int result[40][260];

void search_by_selected() {

	ofstream outp("total_result.txt");

	outp << "fn/cl"<<'\t';
	for (int i = 1; i <= 256; i++) {
		outp << i << '\t';
	}
	outp << endl;

	//k max는 시도하는 selected 개수
	for (int k = 0; k < 40; k++) {
		int x = c[k][0];
		//해당 function이 존재하는 gene number를 list에 뽑아냄
		vector <int> list;
		while (list.size() > 0) {
			list.pop_back();
		}
		for (int i = 1; i <= ROW; i++) {
			int j = 0;
			while (a[i][j] != 0) {
				if (a[i][j] == x) {
					list.push_back(i);
				}
				j++;
			}
		}
		//list의 gene을 토대로 cluster 수를 뽑아냄
		int sum = 0;
		for (int i = 0; i < list.size(); i++) {
			for (int j = 0; j < MAX; j++) {
				if (b[j][0] == list[i]) {
					result[k][b[j][1]] += 1;
					sum++;
				}
				if (b[j][0] > list[i]) {
					break;
				}
			}

		}

		//출력. 클러스터 수에 따라 i max 조정
		outp << c[k][0] << '\t';
		for (int i = 1; i <= 256; i++) {
			outp << result[k][i] << '\t';
		}
		outp << sum << endl;

		outp << "ratio" << '\t';
		for (int i = 1; i <= 256; i++) {
			if (clusize[i] > 0) {
				float rat = ((float)result[k][i] / (float)clusize[i])*100;
				outp << setprecision(4)<<rat<< '\t';
			}
			else {
				outp << 0 << '\t';
			}
		}
		outp << endl;
	

	}
	
	//for (int i = 0; i < 2; i++) {
	//	for (int j = 0; j <= 256; j++) {
	//		cout << result[i][j] << '\t';
	//	}
	//	cout << endl;
	//}
	outp.close();

}


int main()
{

	//make_dense();

	//count_label();

	make_abc();
	count_clu_mem();
	//search_by_labelnum();
	//search_by_selected();






}

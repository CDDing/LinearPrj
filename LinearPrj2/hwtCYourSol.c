#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <Windows.h>
#include "MATRIX_METHODS.h"
//비트맵 헤더를 한묶음으로
typedef struct tagBITMAPHEADER {
	BITMAPFILEHEADER bf;
	BITMAPINFOHEADER bi;
	RGBQUAD hRGB[256]; //이 코드에서는 필요없음 (8bit에만 필요)
}BITMAPHEADER;

//비트맵을 읽어와서 화소정보의 포인터를 리턴
BYTE* loadBitmapFile(int bytesPerPixel, BITMAPHEADER* bitmapHeader, int* imgWidth, int* imgHeight, char* filename);

//비트맵 파일 쓰기
void writeBitmapFile(int bytesPerPixel, BITMAPHEADER outputHeader, BYTE* output, int imgSize, char* filename);
void test();
int main() {
	/*******************************************************************/
	/*************************** Read image  ***************************/
	/*******************************************************************/
	BITMAPHEADER originalHeader;	//비트맵의 헤더부분을 파일에서 읽어 저장할 구조체
	BITMAPHEADER outputHeader;		//변형을 가한 헤더부분을 저장할 구조체
	int imgSize, imgWidth, imgHeight;					//이미지의 크기를 저장할 변수
	int bytesPerPixel = 3;			//number of bytes per pixel (1 byte for R,G,B respectively)

	BYTE* image = loadBitmapFile(bytesPerPixel, &originalHeader, &imgWidth, &imgHeight, "Pikachu.bmp"); //비트맵파일을 읽어 화소정보를 저장 (불러들이는 이미지는 .c와 같은 폴더에 저장)
	if (image == NULL) return 0;

	imgSize = imgWidth * imgHeight; // total number of pixels
	BYTE* output = (BYTE*)malloc(bytesPerPixel * sizeof(BYTE) * imgSize);				//결과값을 저장할 포인터 선언 및 메모리 할당
	outputHeader = originalHeader;										//헤더정보를 출력헤더정보에 할당





	/*******************************************************************/
	/************************ Perform HWT/IHWT *************************/
	/*******************************************************************/
	//이미지 행렬 A 구성 (RGB값이 있으므로 픽셀당 값 하나씩만 읽어서 imgWidth x imgHeight 행렬 구성)
	double** A; //original image matrix
	A = (double**)malloc(sizeof(double*) * imgHeight);
	for (int i = 0; i < imgHeight; i++) {
		A[i] = (double*)malloc(sizeof(double) * imgWidth);
	}
	for (int i = 0; i < imgHeight; i++)
		for (int j = 0; j < imgWidth; j++)
			A[i][j] = image[(i * imgWidth + j) * bytesPerPixel];
	//printf("%d %lf %d %lf", A[0][0], 0.0/0.0, A[0][0]*1, A[0][0]*0.1);
	//Haar matrix H 구성 (orthonormal column을 갖도록 구성)
	int n = imgHeight; //이미지가 정사각형(Height==Width)이라고 가정; n = 2^t,t=0,1,2,...
	//...
	double** H = normalizeEachColumn(constructHaarMatrixRecursive(n), n, n);

	//HWT 수행: 행렬곱 B = H'*A*H
	//...
	double** B = multiplyTwoMatrices(transposeMatrix(H, n, n), n, n, A, n, n);
	B = multiplyTwoMatrices(B, n, n, H, n, n);
	//행렬 B 자르기: B의 upper left corner(subsquare matrix)를 잘라 Bhat에 저장
	//Bhat은 B와 사이즈가 같으며 B에서 잘라서 저장한 부분 외에는 모두 0으로 채워짐
	//...
	//...
	for (int cnt = 1; cnt < n; cnt *= 2) {
		char str[1000] = "output";
		char str2[100];
		_itoa(cnt, str2, 10);
		strcat(str, str2);
		strcat(str, ".bmp");
		//행렬 B자르기
		double** Bhat = allocateMemory(n, n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i < cnt && j < cnt) {
					Bhat[i][j] = B[i][j];
				}
				else {
					Bhat[i][j] = 0;
				}
			}
		}
		//IHWT 수행: Ahat = H*Bhat*H'
		double** Ahat = multiplyTwoMatrices(H, n, n, Bhat, n, n);
		Ahat = multiplyTwoMatrices(Ahat, n, n, transposeMatrix(H, n, n), n, n);
		/*******************************************************************/
		/******************* Write reconstructed image  ********************/
		/*******************************************************************/
		//Ahat을 이용해서 위의 image와 같은 형식이 되도록 구성 (즉, Ahat = [a b;c d]면 [a a a b b b c c c d d d]를 만들어야 함)
		BYTE* Are = (BYTE*)malloc(bytesPerPixel * sizeof(BYTE) * imgSize);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < bytesPerPixel; k++) {
					Are[(i * n + j) * bytesPerPixel + k] = Ahat[i][j];
				}
			}
		}

		//...
		writeBitmapFile(bytesPerPixel, outputHeader, Are, imgSize, str);

		//1,2번 문제 끝


		free(Are);
	}
	//3-c
	double** HT = transposeMatrix(H, n, n);
	double** HL = DivideMatrix(HT, 0, 0, n / 2, n);
	double** HH = DivideMatrix(HT, n / 2, 0, n, n);
	double** HLTHL = multiplyTwoMatrices(transposeMatrix(HL, n / 2, n), n, n / 2, HL, n / 2, n);
	double** HHTHH = multiplyTwoMatrices(transposeMatrix(HH, n / 2, n), n, n / 2, HH, n / 2, n);
	double** HBHT1 = multiplyTwoMatrices(HLTHL, n, n, A, n, n);
	HBHT1 = multiplyTwoMatrices(HBHT1, n, n, HLTHL, n, n);
	BYTE* Are_1 = (BYTE*)malloc(bytesPerPixel * sizeof(BYTE) * imgSize);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < bytesPerPixel; k++) {
				Are_1[(i * n + j) * bytesPerPixel + k] = HBHT1[i][j];
			}
		}
	}
	char HBHT_1[1000] = "HBHT1_output";
	strcat(HBHT_1, ".bmp");
	writeBitmapFile(bytesPerPixel, outputHeader, Are_1, imgSize, HBHT_1);


	double** HBHT2 = multiplyTwoMatrices(HLTHL, n, n, A, n, n);
	HBHT2 = multiplyTwoMatrices(HBHT2, n, n, HHTHH, n, n);

	BYTE* Are_2 = (BYTE*)malloc(bytesPerPixel * sizeof(BYTE) * imgSize);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < bytesPerPixel; k++) {
				Are_2[(i * n + j) * bytesPerPixel + k] = HBHT2[i][j];
			}
		}
	}
	char HBHT_2[1000] = "HBHT2_output";
	strcat(HBHT_2, ".bmp");
	writeBitmapFile(bytesPerPixel, outputHeader, Are_2, imgSize, HBHT_2);


	double** HBHT3 = multiplyTwoMatrices(HHTHH, n, n, A, n, n);
	HBHT3 = multiplyTwoMatrices(HBHT3, n, n, HLTHL, n, n);

	BYTE* Are_3 = (BYTE*)malloc(bytesPerPixel * sizeof(BYTE) * imgSize);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < bytesPerPixel; k++) {
				Are_3[(i * n + j) * bytesPerPixel + k] = HBHT3[i][j];
			}
		}
	}
	char HBHT_3[1000] = "HBHT3_output";
	strcat(HBHT_3, ".bmp");
	writeBitmapFile(bytesPerPixel, outputHeader, Are_3, imgSize, HBHT_3);


	double** HBHT4 = multiplyTwoMatrices(HHTHH, n, n, A, n, n);
	HBHT4 = multiplyTwoMatrices(HBHT4, n, n, HHTHH, n, n);
	BYTE* Are_4 = (BYTE*)malloc(bytesPerPixel * sizeof(BYTE) * imgSize);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < bytesPerPixel; k++) {
				Are_4[(i * n + j) * bytesPerPixel + k] = HBHT4[i][j];
			}
		}
	}
	char HBHT_4[1000] = "HBHT4_output";
	strcat(HBHT_4, ".bmp");
	writeBitmapFile(bytesPerPixel, outputHeader, Are_4, imgSize, HBHT_4);


	//3-d
	double** HLL = DivideMatrix(HL, 0, 0, n / 4, n);
	double** HLH = DivideMatrix(HL, n / 4, 0, n / 2, n);
	double** HLLTHLL = multiplyTwoMatrices(transposeMatrix(HLL, n / 4, n), n, n / 4, HLL, n / 4, n);
	double** HLHTHLH = multiplyTwoMatrices(transposeMatrix(HLH, n / 4, n), n, n / 4, HLH, n / 4, n);
	double** HBHT_HL1 = multiplyTwoMatrices(HLLTHLL, n, n, A, n, n);
	HBHT_HL1 = multiplyTwoMatrices(HBHT_HL1, n, n, HLLTHLL, n, n);
	BYTE* Are_HL1 = (BYTE*)malloc(bytesPerPixel * sizeof(BYTE) * imgSize);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < bytesPerPixel; k++) {
				Are_HL1[(i * n + j) * bytesPerPixel + k] = HBHT_HL1[i][j];
			}
		}
	}
	char HBHTHL1[1000] = "HBHT_HL1_output";
	strcat(HBHTHL1, ".bmp");
	writeBitmapFile(bytesPerPixel, outputHeader, Are_HL1, imgSize, HBHTHL1);


	double** HBHT_HL2 = multiplyTwoMatrices(HLLTHLL, n, n, A, n, n);
	HBHT_HL2 = multiplyTwoMatrices(HBHT_HL2, n, n, HLHTHLH, n, n);

	BYTE* Are_HL2 = (BYTE*)malloc(bytesPerPixel * sizeof(BYTE) * imgSize);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < bytesPerPixel; k++) {
				Are_HL2[(i * n + j) * bytesPerPixel + k] = HBHT_HL2[i][j];
			}
		}
	}
	char HBHTHL2[1000] = "HBHT_HL2_output";
	strcat(HBHTHL2, ".bmp");
	writeBitmapFile(bytesPerPixel, outputHeader, Are_HL2, imgSize, HBHTHL2);


	double** HBHT_HL3 = multiplyTwoMatrices(HLHTHLH, n, n, A, n, n);
	HBHT_HL3 = multiplyTwoMatrices(HBHT_HL3, n, n, HLLTHLL, n, n);

	BYTE* Are_HL3 = (BYTE*)malloc(bytesPerPixel * sizeof(BYTE) * imgSize);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < bytesPerPixel; k++) {
				Are_HL3[(i * n + j) * bytesPerPixel + k] = HBHT_HL3[i][j];
			}
		}
	}
	char HBHTHL3[1000] = "HBHT_HL3_output";
	strcat(HBHTHL3, ".bmp");
	writeBitmapFile(bytesPerPixel, outputHeader, Are_HL3, imgSize, HBHTHL3);


	double** HBHT_HL4 = multiplyTwoMatrices(HLHTHLH, n, n, A, n, n);
	HBHT_HL4 = multiplyTwoMatrices(HBHT_HL4, n, n, HLHTHLH, n, n);
	BYTE* Are_HL4 = (BYTE*)malloc(bytesPerPixel * sizeof(BYTE) * imgSize);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < bytesPerPixel; k++) {
				Are_HL4[(i * n + j) * bytesPerPixel + k] = HBHT_HL4[i][j];
			}
		}
	}
	char HBHTHL4[1000] = "HBHT_HL4_output";
	strcat(HBHTHL4, ".bmp");
	writeBitmapFile(bytesPerPixel, outputHeader, Are_HL4, imgSize, HBHTHL4);
	free(image);
	free(output);
	for (int i = 0; i < imgHeight; i++)
		free(A[i]);
	free(A);

	free(Are_1);
	free(Are_2);
	free(Are_3);
	free(Are_4);
	free(Are_HL1);
	free(Are_HL2);
	free(Are_HL3);
	free(Are_HL4);
	return 0;
}
BYTE* loadBitmapFile(int bytesPerPixel, BITMAPHEADER* bitmapHeader, int* imgWidth, int* imgHeight, char* filename)
{
	FILE* fp = fopen(filename, "rb");	//파일을 이진읽기모드로 열기
	if (fp == NULL)
	{
		printf("파일로딩에 실패했습니다.\n");	//fopen에 실패하면 NULL값을 리턴
		return NULL;
	}
	else
	{
		fread(&bitmapHeader->bf, sizeof(BITMAPFILEHEADER), 1, fp);	//비트맵파일헤더 읽기
		fread(&bitmapHeader->bi, sizeof(BITMAPINFOHEADER), 1, fp);	//비트맵인포헤더 읽기
		//fread(&bitmapHeader->hRGB, sizeof(RGBQUAD), 256, fp);	//색상팔렛트 읽기 (24bitmap 에서는 존재하지 않음)

		*imgWidth = bitmapHeader->bi.biWidth;
		*imgHeight = bitmapHeader->bi.biHeight;
		int imgSizeTemp = (*imgWidth) * (*imgHeight);	// 이미지 사이즈를 상위 변수에 할당

		printf("Size of image: Width %d   Height %d\n", bitmapHeader->bi.biWidth, bitmapHeader->bi.biHeight);
		BYTE* image = (BYTE*)malloc(bytesPerPixel * sizeof(BYTE) * imgSizeTemp);	//이미지크기만큼 메모리할당

		fread(image, bytesPerPixel * sizeof(BYTE), imgSizeTemp, fp);//이미지 크기만큼 파일에서 읽어오기

		fclose(fp);
		return image;
	}
}



void writeBitmapFile(int bytesPerPixel, BITMAPHEADER outputHeader, BYTE* output, int imgSize, char* filename)
{
	FILE* fp = fopen(filename, "wb");

	fwrite(&outputHeader.bf, sizeof(BITMAPFILEHEADER), 1, fp);
	fwrite(&outputHeader.bi, sizeof(BITMAPINFOHEADER), 1, fp);
	//fwrite(&outputHeader.hRGB, sizeof(RGBQUAD), 256, fp); //not needed for 24bitmap
	fwrite(output, bytesPerPixel * sizeof(BYTE), imgSize, fp);
	fclose(fp);
}
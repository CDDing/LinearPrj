#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <Windows.h>
#include "MATRIX_METHODS.h"
//��Ʈ�� ����� �ѹ�������
typedef struct tagBITMAPHEADER {
	BITMAPFILEHEADER bf;
	BITMAPINFOHEADER bi;
	RGBQUAD hRGB[256]; //�� �ڵ忡���� �ʿ���� (8bit���� �ʿ�)
}BITMAPHEADER;

//��Ʈ���� �о�ͼ� ȭ�������� �����͸� ����
BYTE* loadBitmapFile(int bytesPerPixel, BITMAPHEADER* bitmapHeader, int* imgWidth, int* imgHeight, char* filename);

//��Ʈ�� ���� ����
void writeBitmapFile(int bytesPerPixel, BITMAPHEADER outputHeader, BYTE* output, int imgSize, char* filename);
void test();
int main() {
	/*******************************************************************/
	/*************************** Read image  ***************************/
	/*******************************************************************/
	BITMAPHEADER originalHeader;	//��Ʈ���� ����κ��� ���Ͽ��� �о� ������ ����ü
	BITMAPHEADER outputHeader;		//������ ���� ����κ��� ������ ����ü
	int imgSize, imgWidth, imgHeight;					//�̹����� ũ�⸦ ������ ����
	int bytesPerPixel = 3;			//number of bytes per pixel (1 byte for R,G,B respectively)

	BYTE* image = loadBitmapFile(bytesPerPixel, &originalHeader, &imgWidth, &imgHeight, "Pikachu.bmp"); //��Ʈ�������� �о� ȭ�������� ���� (�ҷ����̴� �̹����� .c�� ���� ������ ����)
	if (image == NULL) return 0;

	imgSize = imgWidth * imgHeight; // total number of pixels
	BYTE* output = (BYTE*)malloc(bytesPerPixel * sizeof(BYTE) * imgSize);				//������� ������ ������ ���� �� �޸� �Ҵ�
	outputHeader = originalHeader;										//��������� ������������ �Ҵ�





	/*******************************************************************/
	/************************ Perform HWT/IHWT *************************/
	/*******************************************************************/
	//�̹��� ��� A ���� (RGB���� �����Ƿ� �ȼ��� �� �ϳ����� �о imgWidth x imgHeight ��� ����)
	double** A; //original image matrix
	A = (double**)malloc(sizeof(double*) * imgHeight);
	for (int i = 0; i < imgHeight; i++) {
		A[i] = (double*)malloc(sizeof(double) * imgWidth);
	}
	for (int i = 0; i < imgHeight; i++)
		for (int j = 0; j < imgWidth; j++)
			A[i][j] = image[(i * imgWidth + j) * bytesPerPixel];
	//printf("%d %lf %d %lf", A[0][0], 0.0/0.0, A[0][0]*1, A[0][0]*0.1);
	//Haar matrix H ���� (orthonormal column�� ������ ����)
	int n = imgHeight; //�̹����� ���簢��(Height==Width)�̶�� ����; n = 2^t,t=0,1,2,...
	//...
	double** H = normalizeEachColumn(constructHaarMatrixRecursive(n), n, n);

	//HWT ����: ��İ� B = H'*A*H
	//...
	double** B = multiplyTwoMatrices(transposeMatrix(H, n, n), n, n, A, n, n);
	B = multiplyTwoMatrices(B, n, n, H, n, n);
	//��� B �ڸ���: B�� upper left corner(subsquare matrix)�� �߶� Bhat�� ����
	//Bhat�� B�� ����� ������ B���� �߶� ������ �κ� �ܿ��� ��� 0���� ä����
	//...
	//...
	for (int cnt = 1; cnt < n; cnt *= 2) {
		char str[1000] = "output";
		char str2[100];
		_itoa(cnt, str2, 10);
		strcat(str, str2);
		strcat(str, ".bmp");
		//��� B�ڸ���
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
		//IHWT ����: Ahat = H*Bhat*H'
		double** Ahat = multiplyTwoMatrices(H, n, n, Bhat, n, n);
		Ahat = multiplyTwoMatrices(Ahat, n, n, transposeMatrix(H, n, n), n, n);
		/*******************************************************************/
		/******************* Write reconstructed image  ********************/
		/*******************************************************************/
		//Ahat�� �̿��ؼ� ���� image�� ���� ������ �ǵ��� ���� (��, Ahat = [a b;c d]�� [a a a b b b c c c d d d]�� ������ ��)
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

		//1,2�� ���� ��


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
	FILE* fp = fopen(filename, "rb");	//������ �����б���� ����
	if (fp == NULL)
	{
		printf("���Ϸε��� �����߽��ϴ�.\n");	//fopen�� �����ϸ� NULL���� ����
		return NULL;
	}
	else
	{
		fread(&bitmapHeader->bf, sizeof(BITMAPFILEHEADER), 1, fp);	//��Ʈ��������� �б�
		fread(&bitmapHeader->bi, sizeof(BITMAPINFOHEADER), 1, fp);	//��Ʈ��������� �б�
		//fread(&bitmapHeader->hRGB, sizeof(RGBQUAD), 256, fp);	//�����ȷ�Ʈ �б� (24bitmap ������ �������� ����)

		*imgWidth = bitmapHeader->bi.biWidth;
		*imgHeight = bitmapHeader->bi.biHeight;
		int imgSizeTemp = (*imgWidth) * (*imgHeight);	// �̹��� ����� ���� ������ �Ҵ�

		printf("Size of image: Width %d   Height %d\n", bitmapHeader->bi.biWidth, bitmapHeader->bi.biHeight);
		BYTE* image = (BYTE*)malloc(bytesPerPixel * sizeof(BYTE) * imgSizeTemp);	//�̹���ũ�⸸ŭ �޸��Ҵ�

		fread(image, bytesPerPixel * sizeof(BYTE), imgSizeTemp, fp);//�̹��� ũ�⸸ŭ ���Ͽ��� �о����

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
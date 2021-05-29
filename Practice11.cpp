#pragma warning(disable:4996)
#include <stdio.h>
#include <stdlib.h>
#include <Windows.h>
#include <math.h>
void InverseImage(BYTE* Img, BYTE* Out, int W, int H)
{
	int ImgSize = W * H;
	for (int i = 0; i < ImgSize; i++)
	{
		Out[i] = 255 - Img[i];
	}
}
void BrightnessAdj(BYTE* Img, BYTE* Out, int W, int H, int Val)
{
	int ImgSize = W * H;
	for (int i = 0; i < ImgSize; i++)
	{
		if (Img[i] + Val > 255)
		{
			Out[i] = 255;
		}
		else if (Img[i] + Val < 0)
		{
			Out[i] = 0;
		}
		else 	Out[i] = Img[i] + Val;
	}
}
void ContrastAdj(BYTE* Img, BYTE* Out, int W, int H, double Val)
{
	int ImgSize = W * H;
	for (int i = 0; i < ImgSize; i++)
	{
		if (Img[i] * Val > 255.0)
		{
			Out[i] = 255;
		}
		else 	Out[i] = (BYTE)(Img[i] * Val);
	}
}

void ObtainHistogram(BYTE* Img, int* Histo, int W, int H)
{
	int ImgSize = W * H;
	for (int i = 0; i < ImgSize; i++) {
		Histo[Img[i]]++;
	}
}

void ObtainAHistogram(int* Histo, int* AHisto)
{
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j <= i; j++) {
			AHisto[i] += Histo[j];
		}
	}
}

void HistogramStretching(BYTE* Img, BYTE* Out, int* Histo, int W, int H)
{
	int ImgSize = W * H;
	BYTE Low, High;
	for (int i = 0; i < 256; i++) {
		if (Histo[i] != 0) {
			Low = i;
			break;
		}
	}
	for (int i = 255; i >= 0; i--) {
		if (Histo[i] != 0) {
			High = i;
			break;
		}
	}
	for (int i = 0; i < ImgSize; i++) {
		Out[i] = (BYTE)((Img[i] - Low) / (double)(High - Low) * 255.0);
	}
}
void HistogramEqualization(BYTE* Img, BYTE* Out, int* AHisto, int W, int H)
{
	int ImgSize = W * H;
	int Nt = W * H, Gmax = 255;
	double Ratio = Gmax / (double)Nt;
	BYTE NormSum[256];
	for (int i = 0; i < 256; i++) {
		NormSum[i] = (BYTE)(Ratio * AHisto[i]);
	}
	for (int i = 0; i < ImgSize; i++)
	{
		Out[i] = NormSum[Img[i]];
	}
}

void Binarization(BYTE* Img, BYTE* Out, int W, int H, BYTE Threshold)
{
	int ImgSize = W * H;
	for (int i = 0; i < ImgSize; i++) {
		if (Img[i] < Threshold) Out[i] = 0;
		else Out[i] = 255;
	}
}

void AverageConv(BYTE* Img, BYTE* Out, int W, int H) // 박스평활화
{
	double Kernel[3][3] = { 0.11111, 0.11111, 0.11111,
										0.11111, 0.11111, 0.11111,
										0.11111, 0.11111, 0.11111 };
	double SumProduct = 0.0;
	for (int i = 1; i < H - 1; i++) { // Y좌표 (행)
		for (int j = 1; j < W - 1; j++) { // X좌표 (열)
			for (int m = -1; m <= 1; m++) { // Kernel 행
				for (int n = -1; n <= 1; n++) { // Kernel 열
					SumProduct += Img[(i + m) * W + (j + n)] * Kernel[m + 1][n + 1];
				}
			}
			Out[i * W + j] = (BYTE)SumProduct;
			SumProduct = 0.0;
		}
	}
}

void GaussAvrConv(BYTE* Img, BYTE* Out, int W, int H) // 가우시안평활화
{
	double Kernel[3][3] = { 0.0625, 0.125, 0.0625,
										0.125, 0.25, 0.125,
										0.0625, 0.125, 0.0625 };
	double SumProduct = 0.0;
	for (int i = 1; i < H - 1; i++) { // Y좌표 (행)
		for (int j = 1; j < W - 1; j++) { // X좌표 (열)
			for (int m = -1; m <= 1; m++) { // Kernel 행
				for (int n = -1; n <= 1; n++) { // Kernel 열
					SumProduct += Img[(i + m) * W + (j + n)] * Kernel[m + 1][n + 1];
				}
			}
			Out[i * W + j] = (BYTE)SumProduct;
			SumProduct = 0.0;
		}
	}
}

void Prewitt_X_Conv(BYTE* Img, BYTE* Out, int W, int H) // Prewitt 마스크 X
{
	double Kernel[3][3] = { -1.0, 0.0, 1.0,
										-1.0, 0.0, 1.0,
										-1.0, 0.0, 1.0 };
	double SumProduct = 0.0;
	for (int i = 1; i < H - 1; i++) { // Y좌표 (행)
		for (int j = 1; j < W - 1; j++) { // X좌표 (열)
			for (int m = -1; m <= 1; m++) { // Kernel 행
				for (int n = -1; n <= 1; n++) { // Kernel 열
					SumProduct += Img[(i + m) * W + (j + n)] * Kernel[m + 1][n + 1];
				}
			}
			// 0 ~ 765  =====> 0 ~ 255
			Out[i * W + j] = abs((long)SumProduct) / 3;
			SumProduct = 0.0;
		}
	}
}

void Prewitt_Y_Conv(BYTE* Img, BYTE* Out, int W, int H) // Prewitt 마스크 Y
{
	double Kernel[3][3] = { -1.0, -1.0, -1.0,
										0.0, 0.0, 0.0,
										1.0, 1.0, 1.0 };
	double SumProduct = 0.0;
	for (int i = 1; i < H - 1; i++) { // Y좌표 (행)
		for (int j = 1; j < W - 1; j++) { // X좌표 (열)
			for (int m = -1; m <= 1; m++) { // Kernel 행
				for (int n = -1; n <= 1; n++) { // Kernel 열
					SumProduct += Img[(i + m) * W + (j + n)] * Kernel[m + 1][n + 1];
				}
			}
			// 0 ~ 765  =====> 0 ~ 255
			Out[i * W + j] = abs((long)SumProduct) / 3;
			SumProduct = 0.0;
		}
	}
}

void Sobel_X_Conv(BYTE* Img, BYTE* Out, int W, int H) // Sobel 마스크 X
{
	double Kernel[3][3] = { -1.0, 0.0, 1.0,
										-2.0, 0.0, 2.0,
										-1.0, 0.0, 1.0 };
	double SumProduct = 0.0;
	for (int i = 1; i < H - 1; i++) { // Y좌표 (행)
		for (int j = 1; j < W - 1; j++) { // X좌표 (열)
			for (int m = -1; m <= 1; m++) { // Kernel 행
				for (int n = -1; n <= 1; n++) { // Kernel 열
					SumProduct += Img[(i + m) * W + (j + n)] * Kernel[m + 1][n + 1];
				}
			}
			// 0 ~ 1020  =====> 0 ~ 255
			Out[i * W + j] = abs((long)SumProduct) / 4;
			SumProduct = 0.0;
		}
	}
}

void Sobel_Y_Conv(BYTE* Img, BYTE* Out, int W, int H) // Sobel 마스크 Y
{
	double Kernel[3][3] = { -1.0, -2.0, -1.0,
										0.0, 0.0, 0.0,
										1.0, 2.0, 1.0 };
	double SumProduct = 0.0;
	for (int i = 1; i < H - 1; i++) { // Y좌표 (행)
		for (int j = 1; j < W - 1; j++) { // X좌표 (열)
			for (int m = -1; m <= 1; m++) { // Kernel 행
				for (int n = -1; n <= 1; n++) { // Kernel 열
					SumProduct += Img[(i + m) * W + (j + n)] * Kernel[m + 1][n + 1];
				}
			}
			// 0 ~ 765  =====> 0 ~ 255
			Out[i * W + j] = abs((long)SumProduct) / 4;
			SumProduct = 0.0;
		}
	}
}

void Laplace_Conv(BYTE* Img, BYTE* Out, int W, int H) // Laplacian 마스크
{
	double Kernel[3][3] = { -1.0, -1.0, -1.0,
										-1.0, 8.0, -1.0,
										-1.0, -1.0, -1.0 };
	double SumProduct = 0.0;
	for (int i = 1; i < H - 1; i++) { // Y좌표 (행)
		for (int j = 1; j < W - 1; j++) { // X좌표 (열)
			for (int m = -1; m <= 1; m++) { // Kernel 행
				for (int n = -1; n <= 1; n++) { // Kernel 열
					SumProduct += Img[(i + m) * W + (j + n)] * Kernel[m + 1][n + 1];
				}
			}
			// 0 ~ 2040  =====> 0 ~ 255
			Out[i * W + j] = abs((long)SumProduct) / 8;
			SumProduct = 0.0;
		}
	}
}

void Laplace_Conv_DC(BYTE* Img, BYTE* Out, int W, int H) // (합이 1인)Laplacian 마스크
{
	double Kernel[3][3] = { -1.0, -1.0, -1.0,
										-1.0, 9.0, -1.0,
										-1.0, -1.0, -1.0 };
	double SumProduct = 0.0;
	for (int i = 1; i < H - 1; i++) { // Y좌표 (행)
		for (int j = 1; j < W - 1; j++) { // X좌표 (열)
			for (int m = -1; m <= 1; m++) { // Kernel 행
				for (int n = -1; n <= 1; n++) { // Kernel 열
					SumProduct += Img[(i + m) * W + (j + n)] * Kernel[m + 1][n + 1];
				}
			}
			//Out[i * W + j] = abs((long)SumProduct) / 8;
			if (SumProduct > 255.0) Out[i * W + j] = 255;
			else if (SumProduct < 0.0) Out[i * W + j] = 0;
			else Out[i * W + j] = (BYTE)SumProduct;
			SumProduct = 0.0;
		}
	}
}

void SaveBMPFile(BITMAPFILEHEADER hf, BITMAPINFOHEADER hInfo,
	RGBQUAD* hRGB, BYTE* Output, int W, int H, const char* FileName)
{
	FILE* fp = fopen(FileName, "wb");
	if (hInfo.biBitCount == 24) {
		fwrite(&hf, sizeof(BYTE), sizeof(BITMAPFILEHEADER), fp);
		fwrite(&hInfo, sizeof(BYTE), sizeof(BITMAPINFOHEADER), fp);
		fwrite(Output, sizeof(BYTE), W * H * 3, fp);
	}
	else if (hInfo.biBitCount == 8) {
		fwrite(&hf, sizeof(BYTE), sizeof(BITMAPFILEHEADER), fp);
		fwrite(&hInfo, sizeof(BYTE), sizeof(BITMAPINFOHEADER), fp);
		fwrite(hRGB, sizeof(RGBQUAD), 256, fp);
		fwrite(Output, sizeof(BYTE), W * H, fp);
	}
	fclose(fp);
}

void swap(BYTE* a, BYTE* b) //swap함수
{
	BYTE temp = *a;
	*a = *b;
	*b = temp;
}

BYTE Median(BYTE* arr, int size) //배열에서 중앙값 추출하는 함수
{
	// 오름차순 정렬
	const int S = size;
	for (int i = 0; i < size - 1; i++) // pivot index
	{
		for (int j = i + 1; j < size; j++) // 비교대상 index
		{
			if (arr[i] > arr[j]) 	swap(&arr[i], &arr[j]);
		}
	}
	return arr[S / 2];
}

BYTE MaxPooling(BYTE* arr, int size) //배열에서 가장 큰 값 추출하는 함수
{
	// 오름차순 정렬
	const int S = size;
	for (int i = 0; i < size - 1; i++) // pivot index
	{
		for (int j = i + 1; j < size; j++) // 비교대상 index
		{
			if (arr[i] > arr[j]) 	swap(&arr[i], &arr[j]);
		}
	}
	return arr[S - 1];
}

BYTE MinPooling(BYTE* arr, int size) //배열에서 가장 작은 값 추출하는 함수
{
	// 오름차순 정렬
	const int S = size;
	for (int i = 0; i < size - 1; i++) // pivot index
	{
		for (int j = i + 1; j < size; j++) // 비교대상 index
		{
			if (arr[i] > arr[j]) 	swap(&arr[i], &arr[j]);
		}
	}
	return arr[0];
}

int push(short* stackx, short* stacky, int arr_size, short vx, short vy, int* top) //m_BlobColoring에서 호출되는 함수
{
	if (*top >= arr_size) return(-1);
	(*top)++;
	stackx[*top] = vx;
	stacky[*top] = vy;
	return(1);
}

int pop(short* stackx, short* stacky, short* vx, short* vy, int* top) //m_BlobColoring에서 호출되는 함수
{
	if (*top == 0) return(-1);
	*vx = stackx[*top];
	*vy = stacky[*top];
	(*top)--;
	return(1);
}


// GlassFire 알고리즘을 이용한 라벨링 함수
void m_BlobColoring(BYTE* CutImage, int height, int width)
{
	int i, j, m, n, top, area, Out_Area, index, BlobArea[1000];
	long k;
	short curColor = 0, r, c;
	//	BYTE** CutImage2;
	Out_Area = 1;


	// 스택으로 사용할 메모리 할당
	short* stackx = new short[height * width];
	short* stacky = new short[height * width];
	short* coloring = new short[height * width];

	int arr_size = height * width;

	// 라벨링된 픽셀을 저장하기 위해 메모리 할당

	for (k = 0; k < height * width; k++) coloring[k] = 0;  // 메모리 초기화

	for (i = 0; i < height; i++)
	{
		index = i * width;
		for (j = 0; j < width; j++)
		{
			// 이미 방문한 점이거나 픽셀값이 255가 아니라면 처리 안함
			if (coloring[index + j] != 0 || CutImage[index + j] != 255) continue;
			r = i; c = j; top = 0; area = 1;
			curColor++;

			while (1)
			{
			GRASSFIRE:
				for (m = r - 1; m <= r + 1; m++)
				{
					index = m * width;
					for (n = c - 1; n <= c + 1; n++)
					{
						//관심 픽셀이 영상경계를 벗어나면 처리 안함
						if (m < 0 || m >= height || n < 0 || n >= width) continue;

						if ((int)CutImage[index + n] == 255 && coloring[index + n] == 0)
						{
							coloring[index + n] = curColor; // 현재 라벨로 마크
							if (push(stackx, stacky, arr_size, (short)m, (short)n, &top) == -1) continue;
							r = m; c = n; area++;
							goto GRASSFIRE;
						}
					}
				}
				if (pop(stackx, stacky, &r, &c, &top) == -1) break;
			}
			if (curColor < 1000) BlobArea[curColor] = area;
		}
	}

	float grayGap = 255.0f / (float)curColor; //curColor는 총 영역의 개수. grayGap이 한단위씩 증가하며 영역을 채움

	// 가장 면적이 넓은 영역을 찾아내기 위함 
	for (i = 1; i <= curColor; i++)
	{
		if (BlobArea[i] >= BlobArea[Out_Area]) Out_Area = i; //BlobArea에는 각 영역의 면적 정보
	}
	// CutImage 배열 클리어~
	for (k = 0; k < width * height; k++) CutImage[k] = 255;

	// coloring에 저장된 라벨링 결과중 (Out_Area에 저장된) 영역이 가장 큰 것만 CutImage에 저장
	for (k = 0; k < width * height; k++)
	{
		if (coloring[k] == Out_Area) CutImage[k] = 0;  // 가장 큰 것만 저장(size filtering)
		//if(BlobArea[coloring[k]] > 500) CutImage[k] = 0; //면적이 500보다 큰 것만 출력
		//CutImage[k] = (unsigned char)(coloring[k] * grayGap);  //라벨링 결과 그대로 출력(각 영역 밝기를 다르게 함으로써)
	}

	delete[] coloring;
	delete[] stackx;
	delete[] stacky;
}
// 라벨링 후 가장 넓은 영역에 대해서만 뽑아내는 코드 포함

void BinaryImageEdgeDetection(BYTE* Bin, BYTE* Out, int W, int H) { //경계출력 함수
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			if (Bin[i * W + j] == 0) //전경화소라면(0이면 가장 어두운 것. 위의 코드에서 가장 큰 동전을 이진화해서 Temp에 저장했음)
			{
				if (!(Bin[(i - 1) * W + j] == 0 && Bin[(i + 1) * W + j] == 0 &&
					Bin[i * W + (j - 1)] == 0 && Bin[i * W + (j + 1)] == 0)) //4방향 화소 중 하나라도 전경이 아니라면(= 경계라면)
					Out[i * W + j] = 255; //경계를 흰색으로 
			}
		}
	}
}

BYTE DetermThGonzalez(int* H) { //Gonzalez 이진화 임계치 자동결정 함수
	BYTE ep = 3;
	BYTE Low, High;
	BYTE Th, NewTh;
	int G1 = 0, G2 = 0, cnt1 = 0, cnt2 = 0, mu1, mu2;
	for (int i = 0; i < 256; i++) {
		if (H[i] != 0) {
			Low = i;
			break;
		}
	}
	for (int i = 255; i >= 0; i--) {
		if (H[i] != 0) {
			High = i;
			break;
		}
	}

	Th = (Low + High) / 2;

	while (1) {
		for (int i = Th + 1; i <= High; i++) {
			G1 += (H[i] * i);
			cnt1 += H[i];
		}
		for (int i = Low; i <= Th; i++) {
			G2 += (H[i] * i);
			cnt2 += H[i];
		}
		if (cnt1 == 0) cnt1 = 1;
		if (cnt2 == 0) cnt2 = 1;
		mu1 = G1 / cnt1;
		mu2 = G2 / cnt2;

		NewTh = (mu1 + mu2) / 2;

		if (abs(NewTh - Th) < ep) {
			Th = NewTh;
			break;
		}
		else {
			Th = NewTh;
		}
		G1 = G2 = cnt1 = cnt2 = 0;
	}
	return Th;
}

void VerticalFlip(BYTE* Img, int W, int H) {

	for (int i = 0; i < H / 2; i++) { //y좌표
		for (int j = 0; j < W; j++) { //x좌표
			swap(&Img[i * W + j], &Img[(H - 1 - i) * W + j]);
		}
	}
}

void HorizontalFlip(BYTE* Img, int W, int H) {

	for (int i = 0; i < W / 2; i++) { //x좌표
		for (int j = 0; j < H; j++) { //y좌표
			swap(&Img[j * W + i], &Img[j * W + (W - 1 - i)]);
		}
	}
}

void Translation(BYTE* Image, BYTE* Output, int W, int H, int Tx, int Ty) { //이미지 이동시키는 함수

	Ty *= -1;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			if ((i + Ty < H && i + Ty >= 0) && (j + Tx < W && j + Tx >= 0))
				Output[(i + Ty) * W + (j + Tx)] = Image[i * W + j];
		}
	}
}

void Scaling(BYTE* Image, BYTE* Output, int W, int H, double SF_X, double SF_Y) {

	//scaling factor가 1보다 크면 크기 증가, 1보다 작으면 크기 감소
	int tmpX, tmpY;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			//tmpX = (int)(j * SF_X); //순방향 사상
			//tmpY = (int)(i * SF_Y); //순방향 사상
			tmpX = (int)(j / SF_X);
			tmpY = (int)(i / SF_Y);
			if (tmpY < H && tmpX < W)
				//Output[tmpY * W + tmpX] = Image[i * W + j]; 순방향 사상
				Output[i * W + j] = Image[tmpY * W + tmpX];
		}
	}
}

void Rotation(BYTE* Image, BYTE* Output, int W, int H, int Angle) {

	int tmpX, tmpY;
	double Radian = Angle * 3.141592 / 180.0;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			tmpX = (int)(cos(Radian) * j + sin(Radian) * i);
			tmpY = (int)(-sin(Radian) * j + cos(Radian) * i);
			//tmpX = (int)(cos(Radian) * j - sin(Radian) * i); //순방향
			//tmpY = (int)(sin(Radian) * j + cos(Radian) * i); //순방향
			if ((tmpY < H && tmpY >= 0) && (tmpX < W && tmpX >= 0))
				Output[i * W + j] = Image[tmpY * W + tmpX];
			//Output[tmpY * W + tmpX] = Image[i * W + j]; //순방향
		}
	}
}


//------------------------------------------- 중 간 고 사 ----------------------------------------------------



void MedianFiltering(BYTE* Image, BYTE* Output, int W, int H, int size) {

	int Length = size;		// 마스크의 한 변의 길이
	int Margin = Length / 2;	//마스크 한변의 길이의 절반 크기를 마진으로 정함
	int WSize = Length * Length;	//마스크 픽셀 전체 개수
	BYTE* temp = (BYTE*)malloc(sizeof(BYTE) * WSize);	//임시 배열 동적할당
	int i, j, m, n;		//for문을 위한 변수 선언

	for (i = Margin; i < H - Margin; i++) {		//마진을 제외한 영상의 세로 축
		for (j = Margin; j < W - Margin; j++) {		//마진을 제외한 영상의 가로 축
			for (m = -Margin; m <= Margin; m++) {	//마스크의 세로축
				for (n = -Margin; n <= Margin; n++) {		//마스크의 가로축
					temp[(m + Margin) * Length + (n + Margin)] = Image[(i + m) * W + j + n];	//임시 배열에 현재 마스크의 픽셀 값을 순서대로 저장
				}
			}
			Output[i * W + j] = Median(temp, WSize);	//Median함수를 통해서 임시배열의 중간값을 Output배열에 저장
		}
	}
	free(temp);		//동적으로 할당한 메모리 해제
}


// Img: 사각형을 그릴 이미지배열, W: 영상 가로사이즈, H: 영상 세로사이즈,
// LU_X: 사각형의 좌측상단 X좌표, LU_Y: 사각형의 좌측상단 Y좌표,
// RD_X: 사각형의 우측하단 X좌표, LU_Y: 사각형의 우측하단 Y좌표.
void DrawRectOutline(BYTE* Img, int W, int H, int LU_X, int LU_Y, int RD_X, int RD_Y)
{
	for (int i = LU_X; i < RD_X; i++) 	Img[LU_Y * W + i] = 255;
	for (int i = LU_X; i < RD_X; i++) 	Img[RD_Y * W + i] = 255;
	for (int i = LU_Y; i < RD_Y; i++) 	Img[i * W + LU_X] = 255;
	for (int i = LU_Y; i < RD_Y; i++) 	Img[i * W + RD_X] = 255;
}

// Img: 가로/세로 라인을 그릴 이미지배열, W: 영상 가로사이즈, H: 영상 세로사이즈,
// Cx: 가로/세로 라인이 교차되는 지점의 X좌표
// Cy: 가로/세로 라인이 교차되는 지점의 Y좌표
void DrawCrossLine(BYTE* Img, int W, int H, int Cx, int Cy)
{
	for (int i = 0; i < W - 1; i++) // horizontal line, 가로 라인 
		Img[Cy * W + i] = 255;
	for (int i = 0; i < H - 1; i++) // vertical line, 세로 라인
		Img[i * W + Cx] = 255;
}

void Obtain2DCenter(BYTE* Image, int W, int H, int* Cx, int* Cy)  //무게중심 구하는 함수
{
	int SumX = 0, SumY = 0;
	int cnt = 0;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			if (Image[i * W + j] == 0) // 동공영역이면
			{
				SumX += j;
				SumY += i;
				cnt++;
			}
		}
	}
	if (cnt == 0) 	cnt = 1;
	*Cx = SumX / cnt;
	*Cy = SumY / cnt;
	//printf("%d    %d\n", Cx, Cy);
}

void Obtain2DBoundingBox(BYTE* Image, int W, int H, int* LUX, int* LUY, int* RDX, int* RDY)
{
	int flag = 0;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			if (Image[i * W + j] == 0) {
				*LUY = i;
				flag = 1;
				break;
			}
		}
		if (flag == 1) break;
	}
	flag = 0;
	for (int i = H - 1; i >= 0; i--) {
		for (int j = 0; j < W; j++) {  //원랜 우측 하단부터 돌아야 하지만 좌측하단부터 해도 결과는 같다
			if (Image[i * W + j] == 0) {
				*RDY = i;
				flag = 1;
				break;
			}
		}
		if (flag == 1) break;
	}
	flag = 0;

	for (int j = 0; j < W; j++) {
		for (int i = 0; i < H; i++) {
			if (Image[i * W + j] == 0) {
				*LUX = j;
				flag = 1;
				break;
			}
		}
		if (flag == 1) break;
	}
	flag = 0;

	for (int j = W - 1; j >= 0; j--) {
		for (int i = 0; i < H; i++) { //하단에서 상단으로 돌아야 하지만 결과는 같으므로...
			if (Image[i * W + j] == 0) {
				*RDX = j;
				flag = 1;
				break;
			}
		}
		if (flag == 1) break;
	}
}

void FillColor(BYTE* Image, int X, int Y, int W, int H, BYTE R, BYTE G, BYTE B)
{
	Image[Y * W * 3 + X * 3] = B; // Blue 성분
	Image[Y * W * 3 + X * 3 + 1] = G; // Green 성분
	Image[Y * W * 3 + X * 3 + 2] = R; // Red 성분
}

void RGB2YCbCr(BYTE* Image, BYTE* Y, BYTE* Cb, BYTE* Cr, int W, int H)
{
	for (int i = 0; i < H; i++) { // Y좌표
		for (int j = 0; j < W; j++) { // X좌표
			Y[i * W + j] = (BYTE)(0.299 * Image[i * W * 3 + j * 3 + 2] + 0.587 * Image[i * W * 3 + j * 3 + 1] + 0.114 * Image[i * W * 3 + j * 3]);
			Cb[i * W + j] = (BYTE)(-0.16874 * Image[i * W * 3 + j * 3 + 2] - 0.3313 * Image[i * W * 3 + j * 3 + 1] + 0.5 * Image[i * W * 3 + j * 3] + 128.0);
			Cr[i * W + j] = (BYTE)(0.5 * Image[i * W * 3 + j * 3 + 2] - 0.4187 * Image[i * W * 3 + j * 3 + 1] - 0.0813 * Image[i * W * 3 + j * 3] + 128.0);
		}
	}
}

int main()
{
	BITMAPFILEHEADER hf; // 14바이트
	BITMAPINFOHEADER hInfo; // 40바이트
	RGBQUAD hRGB[256]; // 1024바이트
	FILE* fp;
	fp = fopen("fruit.bmp", "rb");
	if (fp == NULL) {
		printf("File not found!\n");
		return -1;
	}
	fread(&hf, sizeof(BITMAPFILEHEADER), 1, fp);
	fread(&hInfo, sizeof(BITMAPINFOHEADER), 1, fp);

	int ImgSize = hInfo.biWidth * hInfo.biHeight;
	int H = hInfo.biHeight, W = hInfo.biWidth;
	BYTE* Image;
	BYTE* Output;

	if (hInfo.biBitCount == 24) {					// 트루컬러
		Image = (BYTE*)malloc(ImgSize * 3);
		Output = (BYTE*)malloc(ImgSize * 3);
		fread(Image, sizeof(BYTE), ImgSize * 3, fp);
	}
	else {											// 인덱스(그레이)
		fread(hRGB, sizeof(RGBQUAD), 256, fp);
		Image = (BYTE*)malloc(ImgSize);
		Output = (BYTE*)malloc(ImgSize);
		fread(Image, sizeof(BYTE), ImgSize, fp);
	}

	/*fread(hRGB, sizeof(RGBQUAD), 256, fp);
	BYTE* Image = (BYTE*)malloc(ImgSize);
	BYTE* Temp = (BYTE*)malloc(ImgSize); // 임시배열
	BYTE* Output = (BYTE*)malloc(ImgSize);
	fread(Image, sizeof(BYTE), ImgSize, fp);*/
	fclose(fp);

	int Histo[256] = { 0 };
	int AHisto[256] = { 0 };
	

	//---------------------------------------- 중간고사--------------------------------
	//기하변화 확인
	//VerticalFlip(Image, W, H);
	//HorizontalFlip(Image, W, H);
	//Translation(Image, Output, W, H, 100, 40);
	//Scaling(Image, Output, W, H, 0.7, 0.7); //x, y의 scaling factor가 같으면 유니폼 스케일링
	//Rotation(Image, Output, W, H, 60); //원점을 중심으로 회전


	//Gonzalez 이진화 임계치 자동결정 함수
	/*BYTE Th;
	ObtainHistogram(Image, Histo, W, H);
	Th = DetermThGonzalez(Histo);
	Binarization(Image, Output, W, H, Th);
	SaveBMPFile(hf, hInfo, hRGB, Output, hInfo.biWidth, hInfo.biHeight, "outpute.bmp");*/


	//경계검출 후 이미지에 표시(pupil1.bmp)
	/*Binarization(Image, Temp, W, H, 50);
	InverseImage(Temp, Temp, W, H);
	m_BlobColoring(Temp, H, W);
	for (int i = 0; i < ImgSize; i++) Output[i] = Image[i];
	BinaryImageEdgeDetection(Temp, Output, W, H);
	SaveBMPFile(hf, hInfo, hRGB, Output, hInfo.biWidth, hInfo.biHeight, "output_bin_edge.bmp");*/


	/*ObtainHistogram(Image, Histo, hInfo.biWidth, hInfo.biHeight);
	ObtainAHistogram(Histo, AHisto);
	HistogramEqualization(Image, Output, AHisto, hInfo.biWidth, hInfo.biHeight);*/

	//GaussAvrConv(Image, Output, hInfo.biWidth, hInfo.biHeight);

	/*Sobel_X_Conv(Image, Temp, hInfo.biWidth, hInfo.biHeight);
	Sobel_Y_Conv(Image, Output, hInfo.biWidth, hInfo.biHeight);
	for (int i = 0; i < ImgSize; i++) {
		if (Temp[i] > Output[i]) 	Output[i] = Temp[i];
	}
	Binarization(Output, Output, hInfo.biWidth, hInfo.biHeight, 40);*/
	/*GaussAvrConv(Image, Temp, hInfo.biWidth, hInfo.biHeight);
	Laplace_Conv_DC(Temp, Output, hInfo.biWidth, hInfo.biHeight);*/

	//HistogramStretching(Image, Output, Histo, hInfo.biWidth, hInfo.biHeight);
	//InverseImage(Image, Output, hInfo.biWidth, hInfo.biHeight);
	//BrightnessAdj(Image, Output, hInfo.biWidth, hInfo.biHeight, 70);
	//ContrastAdj(Image, Output, hInfo.biWidth, hInfo.biHeight, 0.5);

	//AverageConv(Image, Output, hInfo.biWidth, hInfo.biHeight);


	// 2주차 과제
	/*SaveBMPFile(hf, hInfo, hRGB, Output, hInfo.biWidth, hInfo.biHeight, "output1.bmp");
	for (int i = 0; i < ImgSize; i++) Output[i] = Image[i] + 50;
	SaveBMPFile(hf, hInfo, hRGB, Output, hInfo.biWidth, hInfo.biHeight, "output2.bmp");
	for (int i = 0; i < ImgSize; i++) Output[i] = 255 - Image[i];
	SaveBMPFile(hf, hInfo, hRGB, Output, hInfo.biWidth, hInfo.biHeight, "output3.bmp");*/

	//-----------------------------------중간고사------------------------------



	/* 5주차 과제(마스크 크기 늘려가며 노이즈 제거)
	MedianFiltering(Image, Output, W, H, 11);
	SaveBMPFile(hf, hInfo, hRGB, Output, hInfo.biWidth, hInfo.biHeight, "output1.bmp");*/


	//6주차 과제(동공 추출)
	/*Binarization(Image, Output, W, H, 30);  //이진화
	InverseImage(Output, Output, W, H);  //반전
	m_BlobColoring(Output, H, W);  //가장큰 영역 라벨링

	int Cx, Cy; //무게 중심 저장
	int LUX, LUY, RDX, RDY;

	Obtain2DCenter(Output, W, H, &Cx, &Cy); // 이진영상의 무게중심 구하기
	Obtain2DBoundingBox(Output, W, H, &LUX, &LUY, &RDX, &RDY); // 이진영상의 외접직사각형 좌표 추출
	DrawCrossLine(Image, W, H, Cx, Cy);
	DrawRectOutline(Image, W, H, LUX, LUY, RDX, RDY);

	SaveBMPFile(hf, hInfo, hRGB, Image, hInfo.biWidth, hInfo.biHeight, "output.bmp");*/


	// 가로 라인을 특정 색상으로 
	/*for (int i = 0; i < W; i++) {
		FillColor(Image, i, 200, W, H, 0, 255, 255);
	}*/

	// (50, 100) ~ (300, 400) 박스 채우기
	/*for (int i = 100; i <= 400; i++) {
		for (int j = 50; j <= 300; j++) {
			FillColor(Image, j, i, W, H, 255, 0, 255);
		}
	}*/

	// 가로 띠 만들기
	/*
	// 초기화
	for (int i = 0; i <  H; i++) {
		for (int j = 0; j < W; j++) {
			Image[i * W * 3 + j * 3] = 0;
			Image[i * W * 3 + j * 3 + 1] = 0;
			Image[i * W * 3 + j * 3 + 2] = 0;
		}
	}
	// y좌표 기준 0~239 (Red)
	for (int i = 0; i < 240; i++) {
		for (int j = 0; j < W; j++) {
			Image[i*W*3 + j*3 + 2] = 255;
		}
	}
	// y좌표 기준 120 ~ 359 (Green)
	for (int i = 120; i < 360; i++) {
		for (int j = 0; j < W; j++) {
			Image[i * W * 3 + j * 3 + 1] = 255;
		}
	}
	// y좌표 기준 240 ~ 479 (Blue)
	for (int i = 240; i < 480; i++) {
		for (int j = 0; j < W; j++) {
			Image[i * W * 3 + j * 3] = 255;
		}
	}
	*/

	// 그라데이션 만들기 (B ~ R)
	/*double wt;
	for (int a = 0; a < 120; a++) {
		for (int i = 0; i < W; i++) {
			wt = i / (double)(W - 1);
			Image[a * W * 3 + i * 3] = (BYTE)(255 * (1.0 - wt));  // Blue
			Image[a * W * 3 + i * 3 + 1] = (BYTE)(255 * (1.0 - wt)); // Green
			Image[a * W * 3 + i * 3 + 2] = (BYTE)(255 * wt); // Red
		}
	}*/


	 //Red값이 큰 화소만 masking (R, G, B 모델 기준)
	/*for (int i = 0; i < H; i++) { // Y좌표
		for (int j = 0; j < W; j++) { // X좌표
			if (Image[i * W * 3 + j * 3 + 2] > 130 && 
				Image[i * W * 3 + j * 3 + 1] < 50 &&
				Image[i * W * 3 + j * 3 + 0] < 100) {
				Output[i * W * 3 + j * 3] = Image[i * W * 3 + j * 3];
				Output[i * W * 3 + j * 3 + 1] = Image[i * W * 3 + j * 3 + 1];
				Output[i * W * 3 + j * 3 + 2] = Image[i * W * 3 + j * 3 + 2];
			}
			else
				Output[i * W * 3 + j * 3] = Output[i * W * 3 + j * 3 +1] = Output[i * W * 3 + j * 3 +2] = 0;
		}
	}*/

	BYTE* Y = (BYTE *)malloc(ImgSize);
	BYTE* Cb = (BYTE*)malloc(ImgSize);
	BYTE* Cr = (BYTE*)malloc(ImgSize);

	RGB2YCbCr(Image, Y, Cb, Cr, W, H);

	// 빨간색 딸기영역만 masking (Y, Cb, Cr 모델 기준)
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			if (Cb[i * W + j] < 140 && Cr[i * W + j] > 190) {
				Output[i * W * 3 + j * 3] = Image[i * W * 3 + j * 3];
				Output[i * W * 3 + j * 3 + 1] = Image[i * W * 3 + j * 3 + 1];
				Output[i * W * 3 + j * 3 + 2] = Image[i * W * 3 + j * 3 + 2];
			}
			else 
				Output[i * W * 3 + j * 3] = Output[i * W * 3 + j * 3 + 1] = Output[i * W * 3 + j * 3 + 2] = 0;
		}
	}

	free(Y);
	free(Cb);
	free(Cr);

	// Y, Cb, Cr 분리 후 출력
	/*fp = fopen("Y.raw", "wb");
	fwrite(Y, sizeof(BYTE), W * H, fp);
	fclose(fp);
	fp = fopen("Cb.raw", "wb");
	fwrite(Cb, sizeof(BYTE), W* H, fp);
	fclose(fp);
	fp = fopen("Cr.raw", "wb");
	fwrite(Cr, sizeof(BYTE), W* H, fp);
	fclose(fp);*/

	SaveBMPFile(hf, hInfo, hRGB, Output, hInfo.biWidth, hInfo.biHeight, "output.bmp");


	free(Image);
	free(Output);
	//free(Temp);
	return 0;
}
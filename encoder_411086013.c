#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#pragma pack(2)
double PI = 3.14159265359;
// 定義DCT dimension
#define DCT_H 8
#define DCT_W 8

// 建立BMP header結構
typedef struct Bmpheader
{
    char identifier[2];      // 0x0000
    unsigned int filesize;   // 0x0002
    unsigned short reserved; // 0x0006
    unsigned short reserved2;
    unsigned int bitmap_dataoffset; // 0x000A
    unsigned int bitmap_headersize; // 0x000E
    unsigned int width;             // 0x0012
    unsigned int height;            // 0x0016
    unsigned short planes;          // 0x001A
    unsigned short bits_perpixel;   // 0x001C
    unsigned int compression;       // 0x001E
    unsigned int bitmap_datasize;   // 0x0022
    unsigned int hresolution;       // 0x0026
    unsigned int vresolution;       // 0x002A
    unsigned int usedcolors;        // 0x002E
    unsigned int importantcolors;   // 0x0032
} Bitmap;

// 建立RGB結構
typedef struct RGB
{
    unsigned char R;
    unsigned char G;
    unsigned char B;
} ImgRGB;

// 建立YCbCr結構
typedef struct YCbCr
{
    double Y;
    double Cb;
    double Cr;
} ImgYCbCr;

// 建立2DCT結構
typedef struct DCTArray
{
    double ***Y;
    double ***Cb;
    double ***Cr;
} DCTArray;

// 建立2DCT 基底結構
typedef struct BasisVector
{
    double Y;
    double Cb;
    double Cr;
} BasisVector;

// 建立複製2DCT資料結構
typedef struct DCTCoefficients
{
    double ***Y;
    double ***Cb;
    double ***Cr;
} DCTCoefficients;

// 定義用於存儲 RLE 數據的結構
typedef struct
{
    int size;
    int *data;
} RLEData;

// 定義 Y, Cb, Cr 的量化矩陣
int Q_Y[8][8] = {
    {16, 11, 10, 16, 24, 40, 51, 61},
    {12, 12, 14, 19, 26, 58, 60, 55},
    {14, 13, 16, 24, 40, 57, 69, 56},
    {14, 17, 22, 29, 51, 87, 80, 62},
    {18, 22, 37, 56, 68, 109, 103, 77},
    {24, 35, 55, 64, 81, 104, 113, 92},
    {49, 64, 78, 87, 103, 121, 120, 101},
    {72, 92, 95, 98, 112, 100, 103, 99}};

int Q_Cb[8][8] = {
    {17, 18, 24, 47, 99, 99, 99, 99},
    {18, 21, 26, 66, 99, 99, 99, 99},
    {24, 26, 56, 99, 99, 99, 99, 99},
    {47, 66, 99, 99, 99, 99, 99, 99},
    {99, 99, 99, 99, 99, 99, 99, 99},
    {99, 99, 99, 99, 99, 99, 99, 99},
    {99, 99, 99, 99, 99, 99, 99, 99},
    {99, 99, 99, 99, 99, 99, 99, 99}};

int Q_Cr[8][8] = {
    {17, 18, 24, 47, 99, 99, 99, 99},
    {18, 21, 26, 66, 99, 99, 99, 99},
    {24, 26, 56, 99, 99, 99, 99, 99},
    {47, 66, 99, 99, 99, 99, 99, 99},
    {99, 99, 99, 99, 99, 99, 99, 99},
    {99, 99, 99, 99, 99, 99, 99, 99},
    {99, 99, 99, 99, 99, 99, 99, 99},
    {99, 99, 99, 99, 99, 99, 99, 99}};

BasisVector basis_vector[DCT_H][DCT_W][DCT_H][DCT_W];

// read header
void readheader(FILE *fp, Bitmap *x)
{
    fread(&x->identifier, sizeof(x->identifier), 1, fp);
    fread(&x->filesize, sizeof(x->filesize), 1, fp);
    fread(&x->reserved, sizeof(x->reserved), 1, fp);
    fread(&x->reserved2, sizeof(x->reserved2), 1, fp);
    fread(&x->bitmap_dataoffset, sizeof(x->bitmap_dataoffset), 1, fp);
    fread(&x->bitmap_headersize, sizeof(x->bitmap_headersize), 1, fp);
    fread(&x->width, sizeof(x->width), 1, fp);
    fread(&x->height, sizeof(x->height), 1, fp);
    fread(&x->planes, sizeof(x->planes), 1, fp);
    fread(&x->bits_perpixel, sizeof(x->bits_perpixel), 1, fp);
    fread(&x->compression, sizeof(x->compression), 1, fp);
    fread(&x->bitmap_datasize, sizeof(x->bitmap_datasize), 1, fp);
    fread(&x->hresolution, sizeof(x->hresolution), 1, fp);
    fread(&x->vresolution, sizeof(x->vresolution), 1, fp);
    fread(&x->usedcolors, sizeof(x->usedcolors), 1, fp);
    fread(&x->importantcolors, sizeof(x->importantcolors), 1, fp);
}

// 建立RGB 2維陣列記憶體
ImgRGB **malloc_2D(int row, int col)
{
    ImgRGB **Array, *Data;
    int i;

    // 為指標陣列分配記憶體空間
    Array = (ImgRGB **)malloc(row * sizeof(ImgRGB *));

    // 為實際儲存圖像資料的連續記憶體空間分配空間
    Data = (ImgRGB *)malloc(row * col * sizeof(ImgRGB));

    // 初始化指標陣列，使其每個元素指向對應的行起始位置
    for (i = 0; i < row; i++, Data += col)
    {
        Array[i] = Data;
    }

    // 返回指向2維陣列首地址的指標
    return Array;
}

// 儲存RGB的值
void saveChannelToFile(const char *filename, ImgRGB **array, int H, int W, char channel)
{
    // 打開文件以便寫入
    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
        printf("Error opening file %s for writing.\n", filename);
        return;
    }

    // 遍歷圖像的每一行和每一列
    for (int i = 0; i < H; i++)
    {
        for (int j = 0; j < W; j++)
        {
            unsigned char value;

            // 根據指定的通道（R、G或B），從陣列中提取相應的值
            switch (channel)
            {
            case 'R':
                value = array[i][j].R;
                break;
            case 'G':
                value = array[i][j].G;
                break;
            case 'B':
                value = array[i][j].B;
                break;
            default:
                fclose(file);
                printf("Invalid channel specified.\n");
                return;
            }
            // 將值寫入文件
            fprintf(file, "%d ", value);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

// 建立YCbCr 2維陣列記憶體
ImgYCbCr **malloc_2D_YCbCr(int row, int col)
{
    ImgYCbCr **Array, *Data;
    int i;
    Array = (ImgYCbCr **)malloc(row * sizeof(ImgYCbCr *));
    Data = (ImgYCbCr *)malloc(row * col * sizeof(ImgYCbCr));
    for (i = 0; i < row; i++, Data += col)
    {
        Array[i] = Data;
    }
    return Array;
}

// 讀入RGB資料
void InputData(FILE *fp, ImgRGB **array, int H, int W, int skip)
{
    int temp;
    char skip_buf[3];
    int i, j;

    // 遍歷每一行和每一列
    for (i = 0; i < H; i++)
    {
        for (j = 0; j < W; j++)
        {
            // 讀取一個字節的藍色分量
            temp = fgetc(fp);
            array[i][j].B = temp;
            

            // 讀取一個字節的綠色分量
            temp = fgetc(fp);
            array[i][j].G = temp;
            

            // 讀取一個字節的紅色分量
            temp = fgetc(fp);
            array[i][j].R = temp;
            
        }

        // 跳過一些額外的資料
        if (skip != 0)
            fread(skip_buf, skip, 1, fp);

        
    }
}

// 儲存圖片大小
void saveImageDimensions(const char *filename, int width, int height)
{
    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
        printf("Error opening file %s for writing.\n", filename);
        return;
    }

    // 將圖片的寬度寫入文件
    fprintf(file, "%d\n", width);

    // 將圖片的寬度寫入文件
    fprintf(file, "%d\n", height);

    fclose(file);
}

// RGB to YCbCr轉換
void convertRGBtoYCbCr(ImgRGB **rgbData, ImgYCbCr **ycbcrData, int H, int W, int addH, int addW)
{   
    

    // 遍歷圖像的每一個像素
    for (int i = 0; i < addH; i++)
    {
        for (int j = 0; j < addW; j++)
        {
            if(i>=H||j>=W){

                // 多餘部分設定為0
                ycbcrData[i][j].Y = 0;
                ycbcrData[i][j].Cb = 0;
                ycbcrData[i][j].Cr = 0;
            }else{
                 // 從RGB數據中獲取紅色、綠色和藍色分量
                double R = rgbData[i][j].R;
                double G = rgbData[i][j].G;
                double B = rgbData[i][j].B;

                // 根據YCbCr的定義，計算Y（亮度）、Cb（藍色色度分量）和Cr（紅色色度分量）
                ycbcrData[i][j].Y = 0.299 * R + 0.587 * G + 0.114 * B;
                ycbcrData[i][j].Cb = 128 - 0.168736 * R - 0.331264 * G + 0.5 * B;
                ycbcrData[i][j].Cr = 128 + 0.5 * R - 0.418688 * G - 0.081312 * B;
            }
        }
    }
}

// 分配DCTArray動態記憶體
void allocateDCTArray(DCTArray *F, int M, int N)
{
    // 為三個顏色通道（紅、綠、藍）分配第一層指標
    F->Y = (double ***)malloc(M * sizeof(double **));
    F->Cb = (double ***)malloc(M * sizeof(double **));
    F->Cr = (double ***)malloc(M * sizeof(double **));
    for (int i = 0; i < M; i++)
    {
        // 針對每個顏色通道，進一步分配第二層指標
        F->Y[i] = (double **)malloc(N * sizeof(double *));
        F->Cb[i] = (double **)malloc(N * sizeof(double *));
        F->Cr[i] = (double **)malloc(N * sizeof(double *));
        
        // 對於每個元素，分配足夠的空間來存儲DCT（離散餘弦變換）轉換後的數據
        for (int j = 0; j < N; j++)
        {
            F->Y[i][j] = (double *)malloc(DCT_H * DCT_W * sizeof(double));
            F->Cb[i][j] = (double *)malloc(DCT_H * DCT_W * sizeof(double));
            F->Cr[i][j] = (double *)malloc(DCT_H * DCT_W * sizeof(double));
        }
    }
}

// 釋放DCTArray記憶體
void freeDCTArray(DCTArray *F, int M, int N)
{
    // Free the individual double arrays
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            free(F->Y[i][j]);
            free(F->Cb[i][j]);
            free(F->Cr[i][j]);
        }
        // Free the double pointer arrays
        free(F->Y[i]);
        free(F->Cb[i]);
        free(F->Cr[i]);
    }
    // free the triple pointer arrays
    free(F->Y);
    free(F->Cb);
    free(F->Cr);
}

// 分配DCTCoefficients動態記憶體
void allocateDCTCoefficients(DCTCoefficients *coeffs, int M, int N)
{
    coeffs->Y = (double ***)malloc(M * sizeof(double **));
    coeffs->Cb = (double ***)malloc(M * sizeof(double **));
    coeffs->Cr = (double ***)malloc(M * sizeof(double **));
    for (int i = 0; i < M; i++)
    {
        coeffs->Y[i] = (double **)malloc(N * sizeof(double *));
        coeffs->Cb[i] = (double **)malloc(N * sizeof(double *));
        coeffs->Cr[i] = (double **)malloc(N * sizeof(double *));
         // 對於每個元素，分配足夠的空間來存儲DCT轉換後複製的數據
        for (int j = 0; j < N; j++)
        {
            coeffs->Y[i][j] = (double *)malloc(DCT_H * DCT_W * sizeof(double));
            coeffs->Cb[i][j] = (double *)malloc(DCT_H * DCT_W * sizeof(double));
            coeffs->Cr[i][j] = (double *)malloc(DCT_H * DCT_W * sizeof(double));
        }
    }
}

// 釋放DCTCoefficients動態記憶體
void freeDCTCoefficients(DCTCoefficients *coeffs, int M, int N)
{
    // Free the individual double arrays
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            free(coeffs->Y[i][j]);
            free(coeffs->Cb[i][j]);
            free(coeffs->Cr[i][j]);
        }
        // Free the double pointer arrays
        free(coeffs->Y[i]);
        free(coeffs->Cb[i]);
        free(coeffs->Cr[i]);
    }
    // free the triple pointer arrays
    free(coeffs->Y);
    free(coeffs->Cb);
    free(coeffs->Cr);
}

// 產生basis
void generateBasisVectors()
{
    for (int u = 0; u < DCT_H; u++)
    {
        for (int v = 0; v < DCT_W; v++)
        {
            // 對於每個基向量，計算DCT轉換矩陣的每個元素
            for (int r = 0; r < DCT_H; r++)
            {
                for (int c = 0; c < DCT_W; c++)
                {
                    // 計算DCT的Basis值
                    double cosR = cos(PI * u * (2.0 * r + 1) / (2.0 * DCT_H));
                    double cosC = cos(PI * v * (2.0 * c + 1) / (2.0 * DCT_W));
                    // 將計算出的值賦予基向量陣列的對應元素
                    basis_vector[r][c][u][v].Y = cosR * cosC;
                    basis_vector[r][c][u][v].Cb = cosR * cosC;
                    basis_vector[r][c][u][v].Cr = cosR * cosC;
                }
            }
        }
    }
}

// 執行2D-DCT
void perform2DDCT(ImgYCbCr **image, DCTArray *F, int DCT_M, int DCT_N)
{
    // 初始化比例因子beta，用於DCT的正規化
    double beta[DCT_H];
    for (int i = 0; i < DCT_H; i++)
    {
        beta[i] = (i == 0) ? 1.0 / sqrt(2) : 1.0;
    }

    for (int m = 0; m < DCT_M; m++)
    {
        // 對每個DCT塊進行轉換
        for (int n = 0; n < DCT_N; n++)
        {
            for (int u = 0; u < DCT_H; u++)
            {
                for (int v = 0; v < DCT_W; v++)
                {
                    double sumY = 0.0, sumCb = 0.0, sumCr = 0.0;
                    // 計算DCT的總和
                    for (int x = 0; x < DCT_H; x++)
                    {
                        for (int y = 0; y < DCT_W; y++)
                        {
                            ImgYCbCr pixel = image[m * DCT_H + x][n * DCT_W + y];
                            // 將圖像數據與基向量相乘並累加
                            sumY += (pixel.Y - 128.0) * basis_vector[x][y][u][v].Y;
                            sumCb += (pixel.Cb - 128.0) * basis_vector[x][y][u][v].Cb;
                            sumCr += (pixel.Cr - 128.0) * basis_vector[x][y][u][v].Cr;
                        }
                    }
                    // 根據DCT公式計算並存儲轉換後的數據
                    F->Y[m][n][u * DCT_W + v] = 2.0 / sqrt(DCT_H * DCT_W) * beta[u] * beta[v] * sumY;
                    F->Cb[m][n][u * DCT_W + v] = 2.0 / sqrt(DCT_H * DCT_W) * beta[u] * beta[v] * sumCb;
                    F->Cr[m][n][u * DCT_W + v] = 2.0 / sqrt(DCT_H * DCT_W) * beta[u] * beta[v] * sumCr;
                }
            }
        }
    }
}    


// 複製2D-DCT後的數值
void copyDCTValues(DCTArray *DCTValues, DCTCoefficients *unquantizedValues, int DCT_M, int DCT_N)
{
    for (int m = 0; m < DCT_M; m++)
    {
        for (int n = 0; n < DCT_N; n++)
        {
            for (int i = 0; i < DCT_H; i++)
            {
                for (int j = 0; j < DCT_W; j++)
                {
                    // 將數值存入以方便計算SQNR
                    unquantizedValues->Y[m][n][i * DCT_W + j] = DCTValues->Y[m][n][i * DCT_W + j];
                    unquantizedValues->Cb[m][n][i * DCT_W + j] = DCTValues->Cb[m][n][i * DCT_W + j];
                    unquantizedValues->Cr[m][n][i * DCT_W + j] = DCTValues->Cr[m][n][i * DCT_W + j];
                }
            }
        }
    }
}

// 執行量化的函數
void quantizeDCT(DCTArray *F, int DCT_M, int DCT_N)
{   

    for (int m = 0; m < DCT_M; m++)
    {
        for (int n = 0; n < DCT_N; n++)
        {
            for (int u = 0; u < DCT_H; u++)
            {
                for (int v = 0; v < DCT_W; v++)
                {
                    F->Y[m][n][u * DCT_W + v] = round(F->Y[m][n][u * DCT_W + v] / Q_Y[u][v]);
                    F->Cb[m][n][u * DCT_W + v] = round(F->Cb[m][n][u * DCT_W + v] / Q_Cb[u][v]);
                    F->Cr[m][n][u * DCT_W + v] = round(F->Cr[m][n][u * DCT_W + v] / Q_Cr[u][v]);
                }
            }
        }
    }
}

// 儲存量化表
void saveQuantizationTableAsASCII(int Q[8][8], const char *filename)
{
    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
        printf("无法打开文件 %s\n", filename);
        return;
    }

    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            fprintf(file, "%d ", Q[i][j]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

// 儲存量化後的資料
void saveQuantizedDataAsBinary(DCTArray *F, int DCT_M, int DCT_N, const char *filename, int component)
{
    FILE *file = fopen(filename, "wb");
    if (file == NULL)
    {
        printf("无法打开文件 %s\n", filename);
        return;
    }

    for (int m = 0; m < DCT_M; m++)
    {
        for (int n = 0; n < DCT_N; n++)
        {
            for (int i = 0; i < DCT_H; i++)
            {
                for (int j = 0; j < DCT_W; j++)
                {
                    double data;
                    // 根據組件類型（Y, Cb, Cr）選擇相應的數據
                    switch (component)
                    {
                    case 0: // Y component
                        data = (double)F->Y[m][n][i * DCT_W + j];
                        break;
                    case 1: // Cb component
                        data = (double)F->Cb[m][n][i * DCT_W + j];
                        break;
                    case 2: // Cr component
                        data = (double)F->Cr[m][n][i * DCT_W + j];
                        break;
                    }
                    // 將數據以二進位形式寫入文件
                    fwrite(&data, sizeof(double), 1, file);
                }
            }
        }
    }
    fclose(file);
}

// 計算量化誤差
void calculateDCTError(DCTCoefficients *unquantizedCoeffs, DCTArray *quantizedCoeffs, int Q_Y[8][8], int Q_Cb[8][8], int Q_Cr[8][8], int DCT_M, int DCT_N)
{   
    

    for (int m = 0; m < DCT_M; m++)
    {
        for (int n = 0; n < DCT_N; n++)
        {
            for (int i = 0; i < DCT_H; i++)
            {
                for (int j = 0; j < DCT_W; j++)
                {
                    int u = i * DCT_W + j; // DCT coefficient index

                    // Y channel error
                    unquantizedCoeffs->Y[m][n][u] -= quantizedCoeffs->Y[m][n][u] * Q_Y[i][j];
                    // Cb channel error
                    unquantizedCoeffs->Cb[m][n][u] -= quantizedCoeffs->Cb[m][n][u] * Q_Cb[i][j];
                    // Cr channel error
                    unquantizedCoeffs->Cr[m][n][u] -= quantizedCoeffs->Cr[m][n][u] * Q_Cr[i][j];
                }
            }
        }
    }
}

// 保存 DCT 誤差值到二進位文件
void saveDCTErrorsToBinary(DCTCoefficients *errors, const char *filename_Y, const char *filename_Cb, const char *filename_Cr, int DCT_M, int DCT_N)
{
    //  Y 通道
    FILE *file_Y = fopen(filename_Y, "wb");
    if (file_Y == NULL)
    {
        printf("open failed %s\n", filename_Y);
        return;
    }

    // Cb 通道
    FILE *file_Cb = fopen(filename_Cb, "wb");
    if (file_Cb == NULL)
    {
        printf("open failed %s\n", filename_Cb);
        fclose(file_Y);
        return;
    }

    // Cr 通道
    FILE *file_Cr = fopen(filename_Cr, "wb");
    if (file_Cr == NULL)
    {
        printf("open failed %s\n", filename_Cr);
        fclose(file_Y);
        fclose(file_Cb);
        return;
    }

    // 遍歷並寫入
    for (int m = 0; m < DCT_M; m++)
    {
        for (int n = 0; n < DCT_N; n++)
        {
            for (int i = 0; i < DCT_H; i++)
            {
                for (int j = 0; j < DCT_W; j++)
                {
                    int u = i * DCT_W + j; // DCT coefficient index
                    double value_Y = (double)errors->Y[m][n][u];
                    double value_Cb = (double)errors->Cb[m][n][u];
                    double value_Cr = (double)errors->Cr[m][n][u];

                    // 將誤差值寫入對應的文件
                    fwrite(&value_Y, sizeof(double), 1, file_Y);
                    fwrite(&value_Cb, sizeof(double), 1, file_Cb);
                    fwrite(&value_Cr, sizeof(double), 1, file_Cr);
                }
            }
        }
    }

    // 關閉文件
    fclose(file_Y);
    fclose(file_Cb);
    fclose(file_Cr);
}

// 計算SQNR
void printSQNR(DCTCoefficients *original, DCTArray *quantized, int Q_Y[8][8], int Q_Cb[8][8], int Q_Cr[8][8], int DCT_M, int DCT_N)
{
    double sqnr_Y, sqnr_Cb, sqnr_Cr;
    int total_samples = DCT_M * DCT_N * DCT_H * DCT_W; // 總樣本數

    // 對每個 8x8 DCT 塊計算 SQNR
    for (int i = 0; i < DCT_H; i++)
    {
        for (int j = 0; j < DCT_W; j++)
        {
            double sum_signal_Y = 0.0, sum_noise_Y = 0.0;
            double sum_signal_Cb = 0.0, sum_noise_Cb = 0.0;
            double sum_signal_Cr = 0.0, sum_noise_Cr = 0.0;

            // 計算信號和噪訊的平方和
            for (int m = 0; m < DCT_M; m++)
            {
                for (int n = 0; n < DCT_N; n++)
                {
                    int u = i * DCT_W + j; // DCT coefficient index
                    double signal_Y = original->Y[m][n][u];
                    double noise_Y = signal_Y - (quantized->Y[m][n][u] * Q_Y[i][j]);

                    double signal_Cb = original->Cb[m][n][u];
                    double noise_Cb = signal_Cb - (quantized->Cb[m][n][u] * Q_Cb[i][j]);

                    double signal_Cr = original->Cr[m][n][u];
                    double noise_Cr = signal_Cr - (quantized->Cr[m][n][u] * Q_Cr[i][j]);

                    sum_signal_Y += signal_Y * signal_Y;
                    sum_noise_Y += noise_Y * noise_Y;

                    sum_signal_Cb += signal_Cb * signal_Cb;
                    sum_noise_Cb += noise_Cb * noise_Cb;

                    sum_signal_Cr += signal_Cr * signal_Cr;
                    sum_noise_Cr += noise_Cr * noise_Cr;
                }
            }

            // 計算功率
            double power_signal_Y = sum_signal_Y / total_samples;
            double power_noise_Y = sum_noise_Y / total_samples;

            double power_signal_Cb = sum_signal_Cb / total_samples;
            double power_noise_Cb = sum_noise_Cb / total_samples;

            double power_signal_Cr = sum_signal_Cr / total_samples;
            double power_noise_Cr = sum_noise_Cr / total_samples;

            // 計算 SQNR
            sqnr_Y = power_noise_Y != 0 ? 10 * log10(power_signal_Y / power_noise_Y) : DBL_MAX;
            sqnr_Cb = power_noise_Cb != 0 ? 10 * log10(power_signal_Cb / power_noise_Cb) : DBL_MAX;
            sqnr_Cr = power_noise_Cr != 0 ? 10 * log10(power_signal_Cr / power_noise_Cr) : DBL_MAX;

            // 打印 SQNR
            printf("Y[%d][%d]: %.2f dB, Cb[%d][%d]: %.2f dB, Cr[%d][%d]: %.2f dB\n", i, j, sqnr_Y, i, j, sqnr_Cb, i, j, sqnr_Cr);
        }
    }
}


// 實施 DPCM 於每個 DCT 塊的第一個係數
void performDPCM(DCTArray *F_q, int DCT_M, int DCT_N)
{
    // 遍歷每個顏色成分（Y, Cb, Cr）
    for (int color = 0; color < 3; color++)
    {
        double ***component;
        switch (color)
        {
        case 0:
            component = F_q->Y;
            break;
        case 1:
            component = F_q->Cb;
            break;
        case 2:
            component = F_q->Cr;
            break;
        }

        // 橫向處理每個塊（由右至左）
        for (int m = 0; m < DCT_M; m++)
        {
            for (int n = DCT_N - 1; n > 0; n--)
            {
                component[m][n][0] -= component[m][n - 1][0];
            }
        }

        // 縱向處理每個塊（由下至上）
        for (int m = DCT_M - 1; m > 0; m--)
        {
            component[m][0][0] -= component[m - 1][0][0];
        }
    }
}

// 定義 Zigzag 排序的順序
int zz_order[64] = {
    0, 1, 5, 6, 14, 15, 27, 28,
    2, 4, 7, 13, 16, 26, 29, 42,
    3, 8, 12, 17, 25, 30, 41, 43,
    9, 11, 18, 24, 31, 40, 44, 53,
    10, 19, 23, 32, 39, 45, 52, 54,
    20, 22, 33, 38, 46, 51, 55, 60,
    21, 34, 37, 47, 50, 56, 59, 61,
    35, 36, 48, 49, 57, 58, 62, 63};

// 執行 Zigzag 排序
void performZigzag(double ***F_q, double ***F_q_zigzag, int DCT_M, int DCT_N)
{
    for (int m = 0; m < DCT_M; m++)
    {
        for (int n = 0; n < DCT_N; n++)
        {
            for (int i = 0; i < 64; i++)
            {
                int r = i / 8;            // 計算行數
                int c = i / 8;            // 計算列數
                int zz_pos = zz_order[i]; // 獲取 Zigzag 順序中的位置
                F_q_zigzag[m][n][zz_pos] = F_q[m][n][r * 8 + c];
            }
        }
    }
}

// 修改後的 performRLE 函數
RLEData performRLE(double *block, int H, int W)
{
    int nzero = 0;     // 連續零的數量
    int capacity = 10; // 動態數組分配的初始容量
    int size = 0;      // RLE 數據的當前大小
    int *rle_data = (int *)malloc(capacity * sizeof(int));

    for (int i = 0; i < H * W; i++)
    {
        if (block[i] != 0)
        {
            if (size + 2 > capacity)
            {
                capacity *= 2;
                rle_data = (int *)realloc(rle_data, capacity * sizeof(int));
            }
            rle_data[size++] = nzero;
            rle_data[size++] = (int)block[i];
            nzero = 0;
        }
        else
        {
            nzero++;
        }
    }

    if (size + 2 > capacity)
    {
        capacity += 2;
        rle_data = (int *)realloc(rle_data, capacity * sizeof(int));
    }
    rle_data[size++] = 0;
    rle_data[size++] = 0;

    RLEData result;
    result.size = size;
    result.data = rle_data;
    return result;
}

// 修改後的 encodeImageWithRLE 函數調用
RLEData *encodeImageWithRLE(DCTArray *F_q_zigzag, int DCT_M, int DCT_N)
{
    RLEData *encodedData = (RLEData *)malloc(DCT_M * DCT_N * 3 * sizeof(RLEData));

    int index = 0;
    for (int m = 0; m < DCT_M; m++)
    {
        for (int n = 0; n < DCT_N; n++)
        {
            encodedData[index++] = performRLE(F_q_zigzag->Y[m][n], DCT_H, DCT_W);
            encodedData[index++] = performRLE(F_q_zigzag->Cb[m][n], DCT_H, DCT_W);
            encodedData[index++] = performRLE(F_q_zigzag->Cr[m][n], DCT_H, DCT_W);
        }
    }

    return encodedData;
}

// 釋放 RLE 編碼數據的函數
void freeRLEData(RLEData *data, int M, int N)
{
    int totalBlocks = M * N * 3;
    for (int i = 0; i < totalBlocks; i++)
    {
        free(data[i].data);
    }
    free(data);
}

// 函數來寫入 RLE 數據到文件
void writeRLEDataToFile(const char *filename, RLEData *rleData, int DCT_M, int DCT_N, int imgHeight, int imgWidth)
{
    FILE *file = fopen(filename, "w");
    if (file == NULL)
    {
        printf("Error opening file %s.\n", filename);
        return;
    }

    // 寫入圖片大小
    fprintf(file, "%d %d\n", imgHeight, imgWidth);

    // 每個 8x8 塊的 RLE 數據
    int index = 0;
    for (int m = 0; m < DCT_M; m++)
    {
        for (int n = 0; n < DCT_N; n++)
        {
            // Y, Cb, Cr 順序輸出
            for (int channel = 0; channel < 3; channel++)
            {
                fprintf(file, "(%d,%d, %s) ", m, n, channel == 0 ? "Y" : (channel == 1 ? "Cb" : "Cr"));
                for (int i = 0; i < rleData[index].size; i += 2)
                {
                    fprintf(file, "%d %d ", rleData[index].data[i], rleData[index].data[i + 1]);
                }
                fprintf(file, "\n");
                index++;
            }
        }
    }

    fclose(file);
}

int main(int argc, char **argv)
{
    int mode = atoi(argv[1]);

    if (mode == 0)
    {

        // 處理BMP檔案
        char *fn_in = argv[2];
        char *fn_out_R = argv[3];
        char *fn_out_G = argv[4];
        char *fn_out_B = argv[5];
        char *fn_out_dim = argv[6];

        FILE *fp_in = fopen(fn_in, "rb");
        if (!fp_in)
        {
            printf("Type input\n");
            return 1;
        }

        // 讀取BMP檔頭
        Bitmap bmpheader;
        readheader(fp_in, &bmpheader);

        // 圖像尺寸
        int H = bmpheader.height;
        int W = bmpheader.width;

        // 計算跳過的字節
        int skip = (4 - (bmpheader.width * 3) % 4) % 4;

        // 為RGB數據分配內存
        ImgRGB **Data_RGB = malloc_2D(H, W);
        InputData(fp_in, Data_RGB, H, W, skip);
        fclose(fp_in);

        // 保存R, G, B通道數據
        saveChannelToFile(fn_out_R, Data_RGB, H, W, 'R');
        saveChannelToFile(fn_out_G, Data_RGB, H, W, 'G');
        saveChannelToFile(fn_out_B, Data_RGB, H, W, 'B');

        // 保存圖像尺寸
        saveImageDimensions(fn_out_dim, W, H);

        // 釋放內存
        free(Data_RGB[0]);
        free(Data_RGB);

        printf("BMP processing completed.\n");
    }
    else if (mode == 1)
    {
        // 檢查指令
        if (mode != 1)
        {
            printf("Unsupported mode. Currently, only mode 1 is supported.\n");
            return 1;
        }

        char *input_file = argv[2];
        char *qt_y_file = argv[3];
        char *qt_cb_file = argv[4];
        char *qt_cr_file = argv[5];
        char *dim_file = argv[6];
        char *qf_y_file = argv[7];
        char *qf_cb_file = argv[8];
        char *qf_cr_file = argv[9];
        char *ef_y_file = argv[10];
        char *ef_cb_file = argv[11];
        char *ef_cr_file = argv[12];

        // 開啟輸入bmp
        FILE *fp_in = fopen(input_file, "rb");
        if (!fp_in)
        {
            printf("Cannot open input file %s.\n", input_file);
            return 1;
        }

        // 讀取bmp標頭
        Bitmap bmpheader;
        readheader(fp_in, &bmpheader);

        int H = bmpheader.height;
        int W = bmpheader.width;
        int addH=H;
        int addW=W;
        if(H%8>0){
           addH+=8-(H%8);
        }

        if(W%8>0){
           addW+=8-(W%8);
        }
        int skip = (4 - (W * 3) % 4) % 4;

        ImgRGB **Data_RGB = malloc_2D(H, W);
        ImgYCbCr **Data_YCbCr = malloc_2D_YCbCr(addH, addW);

        // 讀入bmp檔
        InputData(fp_in, Data_RGB, H, W, skip);
        fclose(fp_in);

        // 保存圖像尺寸
        saveImageDimensions(dim_file, W, H);

        convertRGBtoYCbCr(Data_RGB, Data_YCbCr, H, W, addH, addW);

        // 計算 DCT_M 和 DCT_N
        int DCT_M = addH / DCT_H;
        int DCT_N = addW / DCT_W;

        // 產生 DCT 基向量
        generateBasisVectors();

        // 分配記憶體給 DCT 係數
        DCTArray F;
        allocateDCTArray(&F, DCT_M, DCT_N);
        DCTCoefficients unquantizedDCTCoeffs;
        allocateDCTCoefficients(&unquantizedDCTCoeffs, DCT_M, DCT_N);

        // 進行 2D DCT 轉換
        perform2DDCT(Data_YCbCr, &F, DCT_M, DCT_N);

        // 複製DCT值用來計算誤差值
        copyDCTValues(&F, &unquantizedDCTCoeffs, DCT_M, DCT_N);

        // 進行量化
        quantizeDCT(&F, DCT_M, DCT_N);

        // 印出SQNR
        printSQNR(&unquantizedDCTCoeffs, &F, Q_Y, Q_Cb, Q_Cr, DCT_M, DCT_N);

        // 計算量化誤差
        calculateDCTError(&unquantizedDCTCoeffs, &F, Q_Y, Q_Cb, Q_Cr, DCT_M, DCT_N);

        // 儲存量化表
        saveQuantizationTableAsASCII(Q_Y, qt_y_file);
        saveQuantizationTableAsASCII(Q_Cb, qt_cb_file);
        saveQuantizationTableAsASCII(Q_Cr, qt_cr_file);

        // 儲存量化值
        saveQuantizedDataAsBinary(&F, DCT_M, DCT_N, qf_y_file, 0);
        saveQuantizedDataAsBinary(&F, DCT_M, DCT_N, qf_cb_file, 1);
        saveQuantizedDataAsBinary(&F, DCT_M, DCT_N, qf_cr_file, 2);

        // 儲存量化誤差
        saveDCTErrorsToBinary(&unquantizedDCTCoeffs, ef_y_file, ef_cb_file, ef_cr_file, DCT_M, DCT_N);

        free(Data_RGB[0]);
        free(Data_RGB);
        free(Data_YCbCr[0]);
        free(Data_YCbCr);
        freeDCTArray(&F, DCT_M, DCT_N);
        freeDCTCoefficients(&unquantizedDCTCoeffs, DCT_M, DCT_N);

        printf("DCT & Quantization completed\n");
    }
    else if (mode == 2)
    {
        // 檢查指令
        if (mode != 2)
        {
            printf("Unsupported mode. Currently, only mode 2 is supported.\n");
            return 1;
        }

        char *input_file = argv[2];
        char *ascii = argv[3];
        char *rlefile = argv[4];

        // 開啟輸入bmp
        FILE *fp_in = fopen(input_file, "rb");
        if (!fp_in)
        {
            printf("Cannot open input file %s.\n", input_file);
            return 1;
        }

        // 讀取bmp標頭
        Bitmap bmpheader;
        readheader(fp_in, &bmpheader);

        int H = bmpheader.height;
        int W = bmpheader.width;
        int addH = H;
        int addW = W;
        if (H % 8 > 0)
        {
            addH += 8 - (H % 8);
        }

        if (W % 8 > 0)
        {
            addW += 8 - (W % 8);
        }
        int skip = (4 - (W * 3) % 4) % 4;

        ImgRGB **Data_RGB = malloc_2D(H, W);
        ImgYCbCr **Data_YCbCr = malloc_2D_YCbCr(addH, addW);

        // 讀入bmp檔
        InputData(fp_in, Data_RGB, H, W, skip);
        fclose(fp_in);

        convertRGBtoYCbCr(Data_RGB, Data_YCbCr, H, W, addH, addW);

        // 計算 DCT_M 和 DCT_N
        int DCT_M = addH / DCT_H;
        int DCT_N = addW / DCT_W;

        // 產生 DCT 基向量
        generateBasisVectors();

        // 分配記憶體給 DCT 係數
        DCTArray F;
        allocateDCTArray(&F, DCT_M, DCT_N);

        // 進行 2D DCT 轉換
        perform2DDCT(Data_YCbCr, &F, DCT_M, DCT_N);

        // 進行量化
        quantizeDCT(&F, DCT_M, DCT_N);

        // 對量化數據執行 DPCM
        performDPCM(&F, DCT_M, DCT_N);

        // 創建並初始化一個 DCTArray 結構來存儲 Zigzag 排序後的數據
        DCTArray F_q_zigzag;
        allocateDCTArray(&F_q_zigzag, DCT_M, DCT_N);

        // 執行 Zigzag 排序
        performZigzag(F.Y, F_q_zigzag.Y, DCT_M, DCT_N);
        performZigzag(F.Cb, F_q_zigzag.Cb, DCT_M, DCT_N);
        performZigzag(F.Cr, F_q_zigzag.Cr, DCT_M, DCT_N);

        // 對整個圖像進行 RLE 編碼
        RLEData *rle_encoded_data = encodeImageWithRLE(&F_q_zigzag, DCT_M, DCT_N);
        // 將 RLE 數據寫入到文件
        writeRLEDataToFile(rlefile, rle_encoded_data, DCT_M, DCT_N, bmpheader.height, bmpheader.width);

        // 釋放 RLE 編碼數據
        freeRLEData(rle_encoded_data, DCT_M, DCT_N);

        free(Data_RGB[0]);
        free(Data_RGB);
        free(Data_YCbCr[0]);
        free(Data_YCbCr);
        freeDCTArray(&F, DCT_M, DCT_N);
        freeDCTArray(&F_q_zigzag, DCT_M, DCT_N);

        printf("RLE completed\n");
    }
    else
    {
        printf("Invalid mode. Use 0 or 1 or 2 \n");
        return 1;
    }

    return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double PI = 3.14159265359;
#define DCT_H 8
#define DCT_W 8
#pragma pack(2)

// 寫入BMP header
void writeBMPHeader(FILE *fp, int width, int height)
{
    /*construct a structure of BMP header*/
    typedef struct
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
    } Bmpheader;

    Bmpheader header;
    header.identifier[0] = 'B';
    header.identifier[1] = 'M';
    header.filesize = 54 + width * height * 3;
    header.reserved = 0;
    header.reserved2 = 0;
    header.bitmap_dataoffset = 54;
    header.bitmap_headersize = 40;
    header.width = width;
    header.height = height;
    header.planes = 1;
    header.bits_perpixel = 24;
    header.compression = 0;
    header.bitmap_datasize = width * height * 3;
    header.hresolution = 2835;
    header.vresolution = 2835;
    header.usedcolors = 0;
    header.importantcolors = 0;
    fwrite(&header, sizeof(Bmpheader), 1, fp);
}

// 定義 BMP 標頭結構
typedef struct Bmpheader
{
    char identifier[2];
    unsigned int filesize;
    unsigned short reserved;
    unsigned short reserved2;
    unsigned int bitmap_dataoffset;
    unsigned int bitmap_headersize;
    unsigned int width;
    unsigned int height;
    unsigned short planes;
    unsigned short bits_perpixel;
    unsigned int compression;
    unsigned int bitmap_datasize;
    unsigned int hresolution;
    unsigned int vresolution;
    unsigned int usedcolors;
    unsigned int importantcolors;
    unsigned int palette;
} Bitmap;

// 定義 RGB 結構
typedef struct RGB
{
    unsigned char R;
    unsigned char G;
    unsigned char B;
} ImgRGB;

// RLE的結構
typedef struct
{
    int m; // Row of the 8x8 block
    int n; // Column of the 8x8 block
    int size;
    int *data;
} RLEBlock;

typedef struct
{
    double ***Y;
    double ***Cb;
    double ***Cr;
} UnRLEData;

typedef struct
{
    RLEBlock **Y;
    RLEBlock **Cb;
    RLEBlock **Cr;
} RLEImage;

// 讀取 BMP 標頭的函數
void readheader(FILE *fp, Bitmap *bmpHeader)
{
    fread(&bmpHeader->identifier, sizeof(bmpHeader->identifier), 1, fp);
    fread(&bmpHeader->filesize, sizeof(bmpHeader->filesize), 1, fp);
    fread(&bmpHeader->reserved, sizeof(bmpHeader->reserved), 1, fp);
    fread(&bmpHeader->reserved2, sizeof(bmpHeader->reserved2), 1, fp);
    fread(&bmpHeader->bitmap_dataoffset, sizeof(bmpHeader->bitmap_dataoffset), 1, fp);
    fread(&bmpHeader->bitmap_headersize, sizeof(bmpHeader->bitmap_headersize), 1, fp);
    fread(&bmpHeader->width, sizeof(bmpHeader->width), 1, fp);
    fread(&bmpHeader->height, sizeof(bmpHeader->height), 1, fp);
    fread(&bmpHeader->planes, sizeof(bmpHeader->planes), 1, fp);
    fread(&bmpHeader->bits_perpixel, sizeof(bmpHeader->bits_perpixel), 1, fp);
    fread(&bmpHeader->compression, sizeof(bmpHeader->compression), 1, fp);
    fread(&bmpHeader->bitmap_datasize, sizeof(bmpHeader->bitmap_datasize), 1, fp);
    fread(&bmpHeader->hresolution, sizeof(bmpHeader->hresolution), 1, fp);
    fread(&bmpHeader->vresolution, sizeof(bmpHeader->vresolution), 1, fp);
    fread(&bmpHeader->usedcolors, sizeof(bmpHeader->usedcolors), 1, fp);
    fread(&bmpHeader->importantcolors, sizeof(bmpHeader->importantcolors), 1, fp);
}

// 分配二維記憶體空間
ImgRGB **malloc_2D(int row, int col)
{
    // 使用 malloc 分配 row 個指向 ImgRGB 指標的記憶體空間
    ImgRGB **array = (ImgRGB **)malloc(row * sizeof(ImgRGB *));
    for (int i = 0; i < row; i++)
    {
        array[i] = (ImgRGB *)malloc(col * sizeof(ImgRGB));
    }
    return array;
}

// 讀取 RGB 數據的函數
void InputData(FILE *fp, ImgRGB **array, int H, int W, int skip)
{
    int temp; // 用於暫時儲存讀取的數據
    char skip_buf[3];

    // 遍歷所有像素的高度
    for (int i = 0; i < H; i++)
    {
        // 遍歷所有像素的寬度
        for (int j = 0; j < W; j++)
        {
            // 讀取一個字元，並儲存到RGB
            temp = fgetc(fp);
            array[i][j].B = temp;
            temp = fgetc(fp);
            array[i][j].G = temp;
            temp = fgetc(fp);
            array[i][j].R = temp;
        }
        if (skip != 0)
        {
            fread(skip_buf, skip, 1, fp);
        }
    }
}

// 讀取量化資料
void readQuantizedData(const char *filename, double ***data, int height, int width)
{
    // 打開指定的檔案，以二進位模式讀取
    FILE *file = fopen(filename, "rb");
    if (file == NULL)
    {
        perror("Error opening file");
        exit(1);
    }

    for (int i = 0; i < height / DCT_H; i++)
    {
        for (int j = 0; j < width / DCT_W; j++)
        {
            // 從檔案中讀取一個 DCT 區塊的量化資料
            fread(data[i][j], sizeof(double), DCT_H * DCT_W, file);
        }
    }

    fclose(file);
}

// 分配三維記憶體空間(YCbCr)
double ***malloc3D(int height, int width)
{
    double ***array = (double ***)malloc(height / DCT_H * sizeof(double **));
    for (int i = 0; i < height / DCT_H; i++)
    {
        array[i] = (double **)malloc(width / DCT_W * sizeof(double *));
        for (int j = 0; j < width / DCT_W; j++)
        {
            array[i][j] = (double *)malloc(DCT_H * DCT_W * sizeof(double));
        }
    }
    return array;
}

// 釋放三維記憶體空間
void free3D(double ***array, int height, int width)
{
    for (int i = 0; i < height / DCT_H; i++)
    {
        for (int j = 0; j < width / DCT_W; j++)
        {
            free(array[i][j]);
        }
        free(array[i]);
    }
    free(array);
}

// 分配二維記憶體空間
unsigned char **malloc2D(int row, int col)
{
    unsigned char **array;
    array = (unsigned char **)malloc(row * sizeof(unsigned char *));
    for (int i = 0; i < row; i++)
    {
        array[i] = (unsigned char *)malloc(col * sizeof(unsigned char));
    }
    return array;
}

// 分配 double 型態的二維陣列
double **malloc2D_double(int row, int col)
{
    double **array = (double **)malloc(row * sizeof(double *));
    for (int i = 0; i < row; i++)
    {
        array[i] = (double *)malloc(col * sizeof(double));
    }
    return array;
}

// 釋放二維記憶體空間
void free2D(unsigned char **array, int row)
{
    for (int i = 0; i < row; i++)
    {
        free(array[i]);
    }
    free(array);
}

// 釋放二維記憶體空間
void free2D_double(double **array, int row)
{
    for (int i = 0; i < row; i++)
    {
        free(array[i]);
    }
    free(array);
}

// 讀取dim資料
void readImageDimensions(const char *filename, int *height, int *width)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        perror("Error opening file");
        exit(1);
    }

    fscanf(file, "%d", width);
    fscanf(file, "%d", height);

    fclose(file);
}

// 讀取量化表
void readQuantizationTable(const char *filename, int table[DCT_H][DCT_W])
{
    FILE *file = fopen(filename, "r");
    if (file == NULL)
    {
        perror("Error opening file");
        exit(1);
    }

    for (int i = 0; i < DCT_H; i++)
    {
        for (int j = 0; j < DCT_W; j++)
        {
            // 從文件中讀取一個整數並存儲到量化表的相應位置
            fscanf(file, "%d", &table[i][j]);
        }
    }

    fclose(file);
}

// 反量化
void unquantizeDCT(double ***quantizedData, int quantizationTable[DCT_H][DCT_W], int height, int width)
{
    // 計算圖像在垂直和水平方向分別包含多少個 DCT
    int M = height / DCT_H;
    int N = width / DCT_W;

    // 遍歷圖像中的每一個 DCT 塊
    for (int m = 0; m < M; m++)
    {
        for (int n = 0; n < N; n++)
        {
            // 遍歷 DCT 塊的每一個元素
            for (int i = 0; i < DCT_H; i++)
            {
                for (int j = 0; j < DCT_W; j++)
                {
                    // 將量化後的數據乘以相應的量化表中的值進行反量化
                    quantizedData[m][n][i * DCT_W + j] *= quantizationTable[i][j];
                }
            }
        }
    }
}

// 反2D-DCT的basis(同DCT)
void generateBasisVectors(double basis_vector[DCT_H][DCT_W][DCT_H][DCT_W])
{
    for (int u = 0; u < DCT_H; u++)
    {
        for (int v = 0; v < DCT_W; v++)
        {
            for (int r = 0; r < DCT_H; r++)
            {
                for (int c = 0; c < DCT_W; c++)
                {
                    double cosR = cos(PI * u * (2.0 * r + 1) / (2.0 * DCT_H));
                    double cosC = cos(PI * v * (2.0 * c + 1) / (2.0 * DCT_W));
                    basis_vector[u][v][r][c] = cosR * cosC;
                }
            }
        }
    }
}

// 反2D-DCT
void inverseDCT(double ***quantizedData, double **image, double basis_vector[DCT_H][DCT_W][DCT_H][DCT_W], int height, int width)
{
    // 計算圖像在垂直和水平方向分別包含多少個 DCT 塊
    int M = height / DCT_H;
    int N = width / DCT_W;

    for (int m = 0; m < M; m++)
    {
        for (int n = 0; n < N; n++)
        {
            for (int x = 0; x < DCT_H; x++)
            {
                for (int y = 0; y < DCT_W; y++)
                {
                    double sum = 0.0;
                    for (int u = 0; u < DCT_H; u++)
                    {
                        for (int v = 0; v < DCT_W; v++)
                        {
                            // 計算係數，考慮到 u=0 或 v=0 的情況
                            double coeff = (u == 0 ? 1.0 / sqrt(2) : 1.0) * (v == 0 ? 1.0 / sqrt(2) : 1.0);
                            // 進行逆轉換並累加總和
                            sum += coeff * quantizedData[m][n][u * DCT_W + v] * basis_vector[u][v][x][y];
                        }
                    }
                    // 將計算出的總和存儲
                    image[m * DCT_H + x][n * DCT_W + y] = sum * 2.0 / sqrt(DCT_H * DCT_W);
                }
            }
        }
    }
}

// 把YCbCr轉回RGB
void YCbCrToRGB(double **Y, double **Cb, double **Cr, unsigned char **R, unsigned char **G, unsigned char **B, int H, int W)
{
    // 遍歷圖像的每一個像素
    for (int i = 0; i < H; i++)
    {
        for (int j = 0; j < W; j++)
        {
            // 取出 YCbCr 各分量的值
            double y = Y[i][j] + 128.0;
            double cb = Cb[i][j] + 128.0;
            double cr = Cr[i][j] + 128.0;

            // 根據 YCbCr 到 RGB 的轉換公式計算 RGB 值
            int r = round(1.402 * (cr - 128) + y);
            int g = round(-0.34414 * (cb - 128) - 0.71414 * (cr - 128) + y);
            int b = round(1.772 * (cb - 128) + y);

            // 限制 RGB 值在 0 到 255 範圍內
            if (r < 0)
            {
                r = 0;
            }
            else if (r > 255)
            {
                r = 255;
            }

            R[i][j] = (unsigned char)r;

            if (g < 0)
            {
                g = 0;
            }
            else if (g > 255)
            {
                g = 255;
            }

            G[i][j] = (unsigned char)g;

            if (b < 0)
            {
                b = 0;
            }
            else if (b > 255)
            {
                b = 255;
            }
            B[i][j] = (unsigned char)b;
        }
    }
}

// 把RGB存成BMP
void writePixelData(FILE *fp, unsigned char **R, unsigned char **G, unsigned char **B, int width, int height)
{
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            unsigned char color[3] = {B[y][x], G[y][x], R[y][x]}; // BMP 使用 BGR 格式
            fwrite(color, sizeof(unsigned char), 3, fp);
        }
    }
}

// 計算MSE
double calculateMSE(ImgRGB **original, unsigned char **processed, int height, int width)
{
    double mseR = 0.0, mseG = 0.0, mseB = 0.0;
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            // 計算各通道的 MSE
            mseR += pow((double)original[i][j].R - (double)processed[i][j], 2);
            mseG += pow((double)original[i][j].G - (double)processed[i][j], 2);
            mseB += pow((double)original[i][j].B - (double)processed[i][j], 2);
        }
    }
    return (mseR + mseG + mseB) / (3 * height * width); // 计算所有通道的平均 MSE
}

// 計算SQNR
double calculateSQNR(ImgRGB **original, unsigned char **processed, int height, int width)
{
    double mse = calculateMSE(original, processed, height, width);
    if (mse == 0)
        return INFINITY;

    double max_signal_value = 255.0; // 对于8位元的RGB值，最大信号值为255
    double sqnr = 10 * log10(pow(max_signal_value, 2) / mse);
    return sqnr;
}

// 寫入每格的資料
void writePixelData0(FILE *fp, FILE *fR, FILE *fG, FILE *fB, int width, int height)
{
    int r, g, b;
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            fscanf(fR, "%d", &r);
            fscanf(fG, "%d", &g);
            fscanf(fB, "%d", &b);
            unsigned char color[3] = {b, g, r}; // BMP uses BGR format
            fwrite(color, sizeof(unsigned char), 3, fp);
        }
    }
}

double ***malloc3D_float(int height, int width)
{
    double ***array = (double ***)malloc(height / DCT_H * sizeof(double **));
    for (int i = 0; i < height / DCT_H; i++)
    {
        array[i] = (double **)malloc(width / DCT_W * sizeof(double *));
        for (int j = 0; j < width / DCT_W; j++)
        {
            array[i][j] = (double *)malloc(DCT_H * DCT_W * sizeof(double));
        }
    }
    return array;
}

void free3D_float(double ***array, int height, int width)
{
    for (int i = 0; i < height / DCT_H; i++)
    {
        for (int j = 0; j < width / DCT_W; j++)
        {
            free(array[i][j]);
        }
        free(array[i]);
    }
    free(array);
}

// 讀取量化誤差數據
void readQuantizationError(const char *filename, double ***errorData, int height, int width)
{
    FILE *file = fopen(filename, "rb");
    if (file == NULL)
    {
        perror("Error opening file");
        exit(1);
    }

    for (int i = 0; i < height / DCT_H; i++)
    {
        for (int j = 0; j < width / DCT_W; j++)
        {
            fread(errorData[i][j], sizeof(double), DCT_H * DCT_W, file);
        }
    }

    fclose(file);
}

// 將量化誤差加到量化數據上
void addQuantizationError(double ***quantizedData, double ***errorData, int quantizationTable[DCT_H][DCT_W], int height, int width)
{
    int M = height / DCT_H;
    int N = width / DCT_W;

    for (int m = 0; m < M; m++)
    {
        for (int n = 0; n < N; n++)
        {
            for (int i = 0; i < DCT_H; i++)
            {
                for (int j = 0; j < DCT_W; j++)
                {
                    int index = i * DCT_W + j;
                    quantizedData[m][n][index] = (double)quantizedData[m][n][index] * quantizationTable[i][j] + errorData[m][n][index];
                }
            }
        }
    }
}

int main(int argc, char *argv[])
{

    int command = atoi(argv[1]);

    if (command == 1)
    {
        if (argc == 11)
        {

            char *outputBMPFile = argv[2];
            char *inputBMPFile = argv[3];
            char *qt_Y_File = argv[4];
            char *qt_Cb_File = argv[5];
            char *qt_Cr_File = argv[6];
            char *dimFile = argv[7];
            char *qf_Y_File = argv[8];
            char *qf_Cb_File = argv[9];
            char *qf_Cr_File = argv[10];

            int height, width;
            int Q_Y[DCT_H][DCT_W], Q_Cb[DCT_H][DCT_W], Q_Cr[DCT_H][DCT_W];

            // 從 dim.txt 讀取圖像尺寸
            readImageDimensions(dimFile, &height, &width);

            int addheight = height;
            int addwidth = width;
            if (height % 8 > 0)
            {
                addheight += 8 - (height % 8);
            }

            if (width % 8 > 0)
            {
                addwidth += 8 - (width % 8);
            }

            double ***quantized_Y = malloc3D(addheight, addwidth);
            double ***quantized_Cb = malloc3D(addheight, addwidth);
            double ***quantized_Cr = malloc3D(addheight, addwidth);

            // 從讀取量化資料
            readQuantizedData(qf_Y_File, quantized_Y, addheight, addwidth);
            readQuantizedData(qf_Cb_File, quantized_Cb, addheight, addwidth);
            readQuantizedData(qf_Cr_File, quantized_Cr, addheight, addwidth);

            // 讀取量化表
            readQuantizationTable(qt_Y_File, Q_Y);
            readQuantizationTable(qt_Cb_File, Q_Cb);
            readQuantizationTable(qt_Cr_File, Q_Cr);

            // 反量化
            unquantizeDCT(quantized_Y, Q_Y, addheight, addwidth);
            unquantizeDCT(quantized_Cb, Q_Cb, addheight, addwidth);
            unquantizeDCT(quantized_Cr, Q_Cr, addheight, addwidth);

            // 生成基向量
            double basis_vector[DCT_H][DCT_W][DCT_H][DCT_W];
            generateBasisVectors(basis_vector);

            // 創建一個用於存儲反量化圖像的陣列（使用 double）
            double **image_Y = malloc2D_double(addheight, addwidth);
            double **image_Cb = malloc2D_double(addheight, addwidth);
            double **image_Cr = malloc2D_double(addheight, addwidth);

            // 反 2D-DCT
            inverseDCT(quantized_Y, image_Y, basis_vector, addheight, addwidth);
            inverseDCT(quantized_Cb, image_Cb, basis_vector, addheight, addwidth);
            inverseDCT(quantized_Cr, image_Cr, basis_vector, addheight, addwidth);

            // 生成 RGB 圖像陣列（仍然使用 unsigned char）
            unsigned char **image_R = malloc2D(height, width);
            unsigned char **image_G = malloc2D(height, width);
            unsigned char **image_B = malloc2D(height, width);

            // 將 YCbCr 轉換為 RGB
            YCbCrToRGB(image_Y, image_Cb, image_Cr, image_R, image_G, image_B, height, width);

            // 寫入新的 BMP 檔案
            FILE *fp = fopen(outputBMPFile, "wb");
            if (fp == NULL)
            {
                fprintf(stderr, "無法創建 BMP 檔案 %s。\n", outputBMPFile);
                return 1;
            }

            writeBMPHeader(fp, width, height);
            writePixelData(fp, image_R, image_G, image_B, width, height);
            fclose(fp);

            // 打開原始 BMP 檔案以計算 SQNR
            FILE *fp_bmp = fopen(inputBMPFile, "rb");
            if (!fp_bmp)
            {
                fprintf(stderr, "無法打開 BMP 檔案 %s。\n", inputBMPFile);
                return 1;
            }

            Bitmap bmpHeader;
            readheader(fp_bmp, &bmpHeader);

            int H = bmpHeader.height;
            int W = bmpHeader.width;
            int skip = (4 - (W * 3) % 4) % 4;

            ImgRGB **originalData_RGB = malloc_2D(H, W);
            InputData(fp_bmp, originalData_RGB, H, W, skip);
            fclose(fp_bmp);

            // 計算SQNR
            double sqnr_R = calculateSQNR(originalData_RGB, image_R, height, width);
            double sqnr_G = calculateSQNR(originalData_RGB, image_G, height, width);
            double sqnr_B = calculateSQNR(originalData_RGB, image_B, height, width);
            printf("SQNR for R channel: %lf dB\n", sqnr_R);
            printf("SQNR for G channel: %lf dB\n", sqnr_G);
            printf("SQNR for B channel: %lf dB\n", sqnr_B);

            // 釋放記憶體
            free2D_double(image_Y, height);
            free2D_double(image_Cb, height);
            free2D_double(image_Cr, height);
            free2D(image_R, height);
            free2D(image_G, height);
            free2D(image_B, height);

            free3D(quantized_Y, height, width);
            free3D(quantized_Cb, height, width);
            free3D(quantized_Cr, height, width);

            for (int i = 0; i < H; i++)
            {
                free(originalData_RGB[i]);
            }
            free(originalData_RGB);

            printf("Image processing completed\n");
        }
        else if (argc == 13)
        {

            char *outputBMPFile = argv[2];
            char *qt_Y_File = argv[3];
            char *qt_Cb_File = argv[4];
            char *qt_Cr_File = argv[5];
            char *dimFile = argv[6];
            char *qf_Y_File = argv[7];
            char *qf_Cb_File = argv[8];
            char *qf_Cr_File = argv[9];
            char *eF_Y_File = argv[10];
            char *eF_Cb_File = argv[11];
            char *eF_Cr_File = argv[12];

            long int height, width;

            int Q_Y[DCT_H][DCT_W], Q_Cb[DCT_H][DCT_W], Q_Cr[DCT_H][DCT_W];

            // 從 dim.txt 讀取圖像尺寸
            readImageDimensions(dimFile, (int *)&height, (int *)&width);

            int addheight = height;
            int addwidth = width;
            if (height % 8 > 0)
            {
                addheight += 8 - (height % 8);
            }

            if (width % 8 > 0)
            {
                addwidth += 8 - (width % 8);
            }

            double ***quantized_Y = malloc3D_float(addheight, addwidth);
            double ***quantized_Cb = malloc3D_float(addheight, addwidth);
            double ***quantized_Cr = malloc3D_float(addheight, addwidth);

            readQuantizedData(qf_Y_File, quantized_Y, addheight, addwidth);
            readQuantizedData(qf_Cb_File, quantized_Cb, addheight, addwidth);
            readQuantizedData(qf_Cr_File, quantized_Cr, addheight, addwidth);

            // 讀取量化表
            readQuantizationTable(qt_Y_File, Q_Y);
            readQuantizationTable(qt_Cb_File, Q_Cb);
            readQuantizationTable(qt_Cr_File, Q_Cr);

            // 讀取量化誤差
            double ***error_Y = malloc3D_float(addheight, addwidth);
            double ***error_Cb = malloc3D_float(addheight, addwidth);
            double ***error_Cr = malloc3D_float(addheight, addwidth);

            readQuantizationError(eF_Y_File, error_Y, addheight, addwidth);
            readQuantizationError(eF_Cb_File, error_Cb, addheight, addwidth);
            readQuantizationError(eF_Cr_File, error_Cr, addheight, addwidth);

            // 將量化誤差加到量化數據上
            addQuantizationError(quantized_Y, error_Y, Q_Y, addheight, addwidth);
            addQuantizationError(quantized_Cb, error_Cb, Q_Cb, addheight, addwidth);
            addQuantizationError(quantized_Cr, error_Cr, Q_Cr, addheight, addwidth);

            // 生成基向量
            double basis_vector[DCT_H][DCT_W][DCT_H][DCT_W];
            generateBasisVectors(basis_vector);

            // 創建一個用於存儲反量化圖像的陣列（使用 double）
            double **image_Y = malloc2D_double(addheight, addwidth);
            double **image_Cb = malloc2D_double(addheight, addwidth);
            double **image_Cr = malloc2D_double(addheight, addwidth);

            // 反 2D-DCT
            inverseDCT(quantized_Y, image_Y, basis_vector, addheight, addwidth);
            inverseDCT(quantized_Cb, image_Cb, basis_vector, addheight, addwidth);
            inverseDCT(quantized_Cr, image_Cr, basis_vector, addheight, addwidth);

            // 生成 RGB 圖像陣列（仍然使用 unsigned char）
            unsigned char **image_R = malloc2D(height, width);
            unsigned char **image_G = malloc2D(height, width);
            unsigned char **image_B = malloc2D(height, width);

            // 將 YCbCr 轉換為 RGB
            YCbCrToRGB(image_Y, image_Cb, image_Cr, image_R, image_G, image_B, height, width);

            // 寫入新的 BMP 檔案
            FILE *fp = fopen(outputBMPFile, "wb");
            if (fp == NULL)
            {
                fprintf(stderr, "無法創建 BMP 檔案 %s。\n", outputBMPFile);
                return 1;
            }
            writeBMPHeader(fp, width, height);
            writePixelData(fp, image_R, image_G, image_B, width, height);
            fclose(fp);

            // 釋放記憶體
            free2D_double(image_Y, height);
            free2D_double(image_Cb, height);
            free2D_double(image_Cr, height);
            free2D(image_R, height);
            free2D(image_G, height);
            free2D(image_B, height);

            free3D_float(quantized_Y, height, width);
            free3D_float(quantized_Cb, height, width);
            free3D_float(quantized_Cr, height, width);
            // 處理完畢後，釋放量化誤差數據的記憶體
            free3D_float(error_Y, height, width);
            free3D_float(error_Cb, height, width);
            free3D_float(error_Cr, height, width);
            printf("Image processing completed\n");
        }
    }
    else if (command == 0)
    {

        FILE *fR = fopen(argv[3], "r");
        FILE *fG = fopen(argv[4], "r");
        FILE *fB = fopen(argv[5], "r");
        FILE *fDim = fopen(argv[6], "r");
        FILE *fp;
        int width, height;

        if (!fR || !fG || !fB || !fDim)
        {
            printf("Error opening files\n");
            return 1;
        }

        // 從 dim.txt 讀取圖像尺寸
        fscanf(fDim, "%d %d", &width, &height);
        fclose(fDim);

        fp = fopen(argv[2], "wb");
        if (!fp)
        {
            printf("Error creating BMP file\n");
            return 1;
        }

        // RGB資料寫入bmp
        writeBMPHeader(fp, width, height);
        writePixelData0(fp, fR, fG, fB, width, height);

        fclose(fR);
        fclose(fG);
        fclose(fB);
        fclose(fp);

        printf("BMP file created successfully\n");
    }
    else
    {
        printf("Invalid command.\n");
    }

    return 0;
}
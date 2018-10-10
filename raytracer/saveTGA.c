// 32 bit image
typedef struct {
    int width, height;
    int** data;
} RBitmap;

int write_truecolor_tga( RBitmap* data ) {
    FILE *fp = fopen("out.tga", "w");
    if (fp == NULL) return 0;

// The image header
    char header[ 18 ] = { 0 }; // char = byte
    header[ 2 ] = 2; // truecolor
    header[ 12 ] = data->width & 0xFF;
    header[ 13 ] = (data->width >> 8) & 0xFF;
    header[ 14 ] = data->height & 0xFF;
    header[ 15 ] = (data->height >> 8) & 0xFF;
    header[ 16 ] = 24; // bits per pixel

    fwrite((const char*)&header, 1, sizeof(header), fp);

// The image data is stored bottom-to-top, left-to-right
    for (int y = data->height -1; y >= 0; y--)
        for (int x = 0; x < data->width; x++)
        {
            char b = (data->data[x][y] & 0x0000FF);
            char g = (data->data[x][y] & 0x00FF00) >> 8;
            char r = (data->data[x][y] & 0xFF0000) >> 16;
            putc((int)(b & 0xFF),fp);
            putc((int)(g & 0xFF),fp);
            putc((int)(r & 0xFF),fp);
        }

// The file footer
    static const char footer[ 26 ] =
            "\0\0\0\0" // no extension area
            "\0\0\0\0" // no developer directory
            "TRUEVISION-XFILE" // yep, this is a TGA file
            ".";
    fwrite((const char*)&footer, 1, sizeof(footer), fp);

    fclose(fp);
    return 1;
}
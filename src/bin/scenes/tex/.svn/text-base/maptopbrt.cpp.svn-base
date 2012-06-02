#include "imageio.h"

int main()
{
    int w, h;
    unsigned char* img  = loadImageRGBA("map.png", &w, &h);

    char* head = "Shape \"heightfield\" \"integer nu\" [%d] \"integer nv\" [%d] \"float Pz\" [\n";
    printf(head, w, h);

    for(int y=0; y<h; ++y)
    {
        for(int x=0; x<w; ++x)
            printf("%f ", (float)img[(y*w+x)*4] / 255.0f);
        printf("\n");
    }
    printf("]\n");
}
